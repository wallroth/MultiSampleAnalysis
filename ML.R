## Machine Learning functions

#make helper functions accessible
source("myFuncs.R")

decoding.SS.evaluate <- function(train, test, params, parallel=F, 
                                 permutations=10, boot.iterations=10000,
                                 CI=c(.95, .99), ...) {
  ## function to evaluate a model fit with certain parameter settings 
  ## used at the core of a parameter search
  ## decoding.SS(.parallel) is called with k=1 and evaluated via data.fit_test
  #INPUT ---
  #train, test: the training and test sets. Either corresponding to a single subject
  #             or to a list of subjects (in which case parallel must be TRUE).
  #params: list of parameters which will be fed to decoding.SS
  #        the names must correspond to valid inputs of decoding.SS
  #        corresponds to the parameter combination of a single instance of a parameter search
  #parallel: if True, data corresponds to multiple subjects and decoding.SS.parallel 
  #          will be called
  #permutations: if > 0, permutation tests will be performed
  #              the labels of the training and test set are shuffled on each permutation,
  #              i.e. both the fit and the validation are performed on random data
  #boot.iterations: if not NULL and permutations > 1, the resulting distribution
  #                 of randomly obtained AUC values will be bootstrapped with the specified
  #                 number of iterations
  #CI: the one-sided confidence interval for the non-parametric test.
  #    only relevant with permutations > 0. These correspond to the bootstrapped values
  #    or directly to the distribution of permuted AUC values (if boot.iterations = NULL)
  #    AUC values larger the CI are conisdered significant
  #...: further arguments corresponding to decoding.SS that are constant across the parameter search
  #     e.g. slice = TRUE if all instances of the search operate with sliced data
  #RETURNS ---
  #a list whose elements depend on the settigns of parallel and permutations.
  #params: parameter supplied via 'params' expanded to a df
  #average: if parallel = T (multiple subjects), an average summary is supplied,
  #         i.e. the mean of the AUC and if permutations > 0, new CI values 
  #         corresponding to the overall distribution of random AUC values
  #summary: result df with AUC, CI values and significance test (separated for sujbects
  #         if parallel is TRUE).
  #permutation.result: if permutations > 0, the randomly obtained AUC values
  #boot.result: boostrapped values if boot.iterations > 0
  #training.fit: the model fits obtained on the training set with true labels
  
  .eval_package("foreach", load=T)
  if (class(params) != "list") stop( "'params' must be a list with elements corresponding ",
                                     "to valid inputs of decoding.SS.")
  if ( any(CI > 1) || any(CI < 0) ) stop( "CI must be in the range 0|1 corresponding to 0% to 100%." )
  
  decodeFun = ifelse(parallel, "decoding.SS.parallel", "decoding.SS")
  args.in = list(...)
  args.in = modifyList( args.in, params )
  checktime = proc.time()[3]
  cat("Starting to decode with true labels . . .")
  decoding.true = do.call( decodeFun, modifyList(args.in, 
                      list( data = train, verbose = T, permute = F, k = 1) ) )
  cat(" | total time elapsed (mins):", (proc.time()[3]-checktime)/60, "\n")
  #save the true model fit for evaluation and output
  if (!parallel) {
    training.fit = decoding.true$fits.true$CV1$Fold1
  } else {
    training.fit = lapply( decoding.true$fits.true, function(sub) sub$CV1$Fold1 )
  }
  fits = list(training.fit)
  if ( permutations > 0 ) {
    training.fit.random = list()
    cat("Starting permutations . . .\n")
    for ( run in seq_len(permutations) ) {
      #shuffle labels
      if (parallel) {
        train.random = lapply(train, function(d) data.permute_labels(d, shuffle=T)$data)
      } else {
        train.random = data.permute_labels(train, shuffle=T)$data
      }
      cat("Permutation:", run)
      #get random model fit
      decoding.random = do.call( decodeFun, modifyList(args.in, 
                            list( data = train.random, verbose = T, permute = F, k = 1) ) )
      cat(" | total time elapsed (mins):", (proc.time()[3]-checktime)/60, "\n")
      #save the random fit for evaluation
      if (parallel) {
        training.fit.random[[run]] = lapply( decoding.random$fits.true, function(sub) sub$CV1$Fold1 )
      } else {
        training.fit.random[[run]] = decoding.random$fits.true$CV1$Fold1
      }
    }
    fits = c( fits, training.fit.random )
  }
  pExport = ifelse( is.null(args.in$model), "LiblineaR", 
                    ifelse( toupper(args.in$model) == "L2", "LiblineaR", "MASS" ) )
  #create evaluation function to run inside foreach
  evaluate <- function(test, fit, args.in) {
    source("myFuncs.R")
    .loadFuncs(util=T, CSP=T, ML=T)
    validation_performance = do.call(data.fit_test,
                                     modifyList( args.in, list( data = test, fit = fit ) ))
    return(validation_performance)
  }
  #save slice argument for evaluate
  args.in$slice = all(grepl( "slice", names(training.fit) ))
  if ( parallel ) {
    args.in$slice = all(grepl( "slice", names(training.fit[[1]]) ))
    subs = seq_along(train)
    if ( is.null(args.in$nCores) ) {
      #initialize all available cores but not more than needed
      args.in$nCores = min( parallel::detectCores(), length(subs) )
    }
    CPUcluster = .parallel_backend(on=T, nCores=args.in$nCores)
  }
  ## begin model evaluation
  test.performance = list()
  test.run = test
  cat("Starting model evaluation . . .\n")
  for ( fitnum in seq_along(fits) ) {
    cat("Evaluation:", fitnum)
    if ( fitnum > 1 ) { #permutation fit
      #shuffle labels
      if (parallel) {
        test.run = lapply(test, function(d) data.permute_labels(d, shuffle=T)$data)
      } else {
        test.run= data.permute_labels(test, shuffle=T)$data
      }
    }
    if (parallel) {
      test.performance[[fitnum]] =
        foreach(sub=subs, .combine=list, .multicombine=T, .packages=pExport) %dopar%
        evaluate( test.run[[sub]], fits[[fitnum]][[sub]], args.in )  
    } else {
      test.performance[[fitnum]] = evaluate( test.run, fits[[fitnum]], args.in )
    }
    cat(" | total time elapsed (mins):", (proc.time()[3]-checktime)/60, "\n")
  }
  ## extract true and random performance
  if (parallel) {
    performance.true = setNames( lapply(subs, function(sub) {
      test.performance[[1]][[sub]][[1]]
    }), paste0("subject", subs) )
    if ( permutations > 0 ) {
      performance.random = setNames( lapply(subs, function(sub) {
        perf = lapply( test.performance[2:length(fits)], "[[", sub )
        if ( args.in$slice ) { #slices
          perf = do.call( rbind, lapply(perf, "[[", "summary") )
        } else { #no slices
          perf = data.frame( AUC = do.call( c, lapply(perf, "[[", "AUC") ) )
        }
      }), paste0("subject", subs) )
    }
  } else {
    performance.true = test.performance[[1]][[1]]
    if ( permutations > 0 ) {
      performance.random = lapply( test.performance[2:length(fits)], "[[", 1 )
      if ( args.in$slice ) { #slices
        performance.random = do.call( rbind, performance.random )
      } else { #no slices
        performance.random = data.frame( AUC = do.call( c, performance.random ) )
      }
    }
  }
  ## compare random vs. true performance
  boot = permutations > 1 && !is.null(boot.iterations) && boot.iterations > 1
  if (boot) {
    cat("Bootstrapping . . .")
    if (parallel) {
      boot.result = foreach(sub=subs, .combine=list, .multicombine=T) %dopar%
        .bootstrap( performance.random[[sub]], iterations=boot.iterations )
      boot.result = setNames( boot.result, paste0("subject", subs) )
    } else {
      boot.result = list( .bootstrap( performance.random, iterations=boot.iterations ) )
    }
  }
  if (parallel) .parallel_backend(on=F, CPUcluster) #close the CPU cluster
  # create summary if permutation tests were done
  if ( permutations > 0 ) {
    if ( !parallel ) { #create list of 1 subject
      performance.random = list(performance.random)
      performance.true = list(performance.true)
    }
    summary = setNames( lapply( seq_along(performance.true), function(sub) {
      #calculate percentile thresholds
      if (boot) { #result is a matrix
        thresholds = t( apply( boot.result[[sub]], 1, quantile, prob=CI ) )
        if ( length(CI) < 2 ) thresholds = matrix(thresholds)
      } else { #result is a df
        if ( args.in$slice ) {
          thresholds = do.call( rbind, lapply( unique(performance.random[[sub]]$slice), function(slice) {
            quantile( performance.random[[sub]][ performance.random[[sub]]$slice == slice, "AUC" ], prob=CI )
          }) )
        } else {
          thresholds = t( quantile( performance.random[[sub]][, "AUC" ], prob=CI ) )
        }
      }
      if ( is.null(names( performance.true[[sub]] )) ) names( performance.true[[sub]] ) = "AUC"
      summary = setNames( data.frame( performance.true[[sub]], data.frame(thresholds), row.names=NULL ), 
                          c( names( performance.true[[sub]] ), paste0("CI", CI*100) ) )
    }), paste0("subject", seq_along(performance.true)) )
    for ( p in CI ) {
      #add 'significance' info to summary
      summary = lapply( summary, function(s) {
        s[[paste0("p",p*100)]] = s$AUC > s[[paste0("CI",p*100)]]
        s
      })
    }
  }
  cat(" Done. Total time elapsed (mins):", (proc.time()[3]-checktime)/60, "\n")
  ## compile output
  output = list( params = expand.grid(params) )
  if ( !args.in$slice ) { #add name to AUC value
    if (parallel) {
      performance.true = lapply(performance.true, setNames, nm="AUC")
    } else {
      performance.true = setNames(performance.true, nm="AUC")
    }
  }
  if (parallel) { #create an average
    average = Reduce("+", performance.true)/length(performance.true)
    if ( permutations > 0 ) {
      avg = do.call(rbind, performance.random)
      if (boot) { #redo bootstrap for average CI
        avgboot = .bootstrap( avg, iterations=boot.iterations )
        thresholds = t( apply( avgboot, 1, quantile, prob=CI ) )
      } else {
        if (args.in$slice) {
          thresholds = do.call( rbind, lapply( unique(avg$slice), function(slice) {
            quantile( avg[ avg$slice == slice, "AUC" ], prob=CI )
          }) )
        } else {
          thresholds = t( quantile( avg[, "AUC" ], prob=CI ) )
        }
      }
      average = setNames( data.frame( average, data.frame(thresholds), row.names=NULL ), 
                          c( names( average ), paste0("CI", CI*100) ) )
      for ( p in CI ) {
        #add 'significance' info
        average[[paste0("p",p*100)]] = average$AUC > average[[paste0("CI",p*100)]]
      }
    } #permutations end
    output$average = average
  }
  if ( permutations > 0 ) {
    if (!parallel) { #remove the nesting in case of only 1 subject
      summary = summary[[1]]
      performance.random = performance.random[[1]]
      if (boot) boot.result = boot.result[[1]]
    }
    output$summary = summary
    output$permutation.result = performance.random
    if (boot) output$boot.result = boot.result
  } else { #no permutation tests
    output$summary = performance.true
  }
  output$training.fit = training.fit
  return(output)
}

.bootstrap <- function(values, iterations=10000) {
  #internal function to bootstrap permutation results
  if ( "slice" %in% names(values) ) {
    boot = replicate(iterations, {
      setNames( sapply( unique(values$slice), function(slice) { #stratify the slices
        mean( sample( values[ values$slice == slice, "AUC" ], replace=T) )
      }), paste0("slice", unique(values$slice)) )
    }) #rows: slice number; columns: iteration
  } else {
    boot = t(as.matrix( replicate(iterations, {
      mean( sample( values[, "AUC" ], replace=T) )
    }) )) #convert vector to matrix with 1 row for compatibility
  }
  return( boot )
}

data.fit_test <- function(data, fit, slice=F, ...) {
  ## a function to evaluate a model fit on "new" data, i.e. a validation or test set
  ## fit is the result of a fitting procedure, e.g. decoding.SS or data.fit_model
  ## it will be used to predict the outcome that is in the 2nd column of data
  ## if fit is a list of slice fits, the data should either be sliced as well
  ## or the appropriate slicing parameters have to be provided (with slice = T)
  ## if CSP was used to create the fit, provide the arguments approximate/logtransform
  ## in case they differed from their defaults
  #INPUT ---
  #data: data frame, list of trials or slices
  #fit: either directly a model fit or the output of decoding.SS / data.fit_model
  #     if fit was "manually" created make sure scale/CSP are included if applicable
  #slice: if True, fit is a list of fits corresponding to slices and data will be 
  #       sliced according to the supplied window/overlap parameters
  #... : further arguments correspond to slicing and CSP feature settings if applicable
  #OUTPUT ---
  #a list: if no slicing, the output directly mimicks data.fit_model
  #        otherwise a list with elements "summary" and "fits"
  
  #evaluate different types of fit input
  if ( "fit" %in% names(fit) ) { #input directly from data.fit_model
    fit = fit$fit 
  } else {
    if ( any( grepl("fits", names(fit)) ) ) { #multiple fits
      temp = fit[[ which( grepl("fits", names(fit)) )[1] ]] #if true/random fits, take the first
    } else if ( class(fit) == "list" ) { 
      #perhaps a sub element was input so no "fits" name present
      temp = fit
    }
    if (slice || .is.sliced(data) ) {
      #find the slice fits
      while ( !all( grepl("slice", names(temp)) ) ) {
        if ( identical(temp, temp[[1]]) ) stop( "Could not find any fits." )
        temp = temp[[1]] #if numcv/folds, take first CV
      }
      fit = temp #list of slice fits
    } else if ( class(fit) == "list" ) {
      #find the fit
      while ( !"fit" %in% names(temp) ) {
        if ( identical(temp, temp[[1]]) ) stop( "Could not find a fit." )
        temp = temp[[1]] #if numcv/folds, take first CV
      }
      fit = temp$fit #single fit
    }
  } # "else" fit is already a model fit, e.g. of class LiblineaR or lda
  
  #evaluate function arguments
  args.slice = .eval_ellipsis("data.trials.slice", ...)
  args.slice$split = F #save RAM
  .loadFuncs(CSP=T)
  args.csp.feat = .eval_ellipsis("CSP.get_features", ...)
  #evaluate data
  if ( .is.sliced(data) ) {
    #unslice back to df
    data = .data.unslice(data, params.out=T)
    if (!slice) { #save settings if no new slice request was made on function call
      slice = T 
      args.slice = modifyList( args.slice, list(window=data$window, overlap=data$overlap) )
    }
    data = data$data
  }
  data = data.check(data, aslist=T, strip=F) #transform to list if necessary
  #evaluation end: data is a list of trials at this point 
  
  ## create prediction function ##
  get_predictions <- function(data, fit) {
    #internal helper function
    data = data.check(data, aslist=F)
    if ( !is.factor(data[,2]) ) data[,2] = as.factor(data[,2])
    if ( "CSP" %in% names( fit ) ) {
      if ( eval(args.csp.feat$approximate) ) { #1 value per trial
        nsamples = data.get_samplenum(data) 
        outcome = data[ seq(1, nrow(data), nsamples), 2]
        n = 1:length(outcome)
      } else { #dimensionality is the same as before
        n = data[,1]
        outcome = data[,2]
      }
      data = data.frame(n, outcome, do.call(CSP.get_features, 
                 modifyList( args.csp.feat, list(data=data, filters=fit$CSP$filters) )))
    }
    if ( "scale" %in% names( fit ) ) {
      data = cbind( data[, !is.datacol(data)], 
                    scale( data[, is.datacol(data)], 
                           center = fit$scale$center, 
                           scale = fit$scale$scale) )
    }
    #get predictions
    preds <- predict( fit, data[, is.datacol(data)] )
    cm = caret::confusionMatrix( preds[[1]], data[,2] )
    acc = AUC::auc( AUC::roc( preds[[1]], data[,2] ) )
    return( list(AUC=acc, CM=cm, fit=fit) )
  }
  
  ## start fit evaluation ##
  if (slice) {
    # iterate over the slices and obtain predictions
    sx = do.call(data.trials.slice, modifyList( args.slice, list(data=data) ))
    slicenums = unique(sx$splitidx) #number of slices
    fits = setNames( lapply( slicenums, function(i) { #iterate slices
      sliceidx = sx$slides[sx$splitidx == i][1:(args.slice$window)]
      slice.data = lapply( data, function(d) na.omit( d[ sliceidx, ] ) ) #imbalanced slices generate NAs
      get_predictions(slice.data, fit[[i]]$fit) # return to lapply
    }), paste0("slice",slicenums) )
    out = list( summary = data.frame(slice=1:length(fits), 
                    AUC = sapply(fits, "[[", "AUC"), row.names=NULL), fits = fits)
  } else { #single fit
    out = get_predictions(data, fit)
  }
  return(out)
}

decoding.GAT <- function(data, slice=T, numcv=1, CSP=F, verbose=F, ...) {
  ## Generalization across time: how long stays a process active, how many are there
  ## train model at one time point and test prediction accuracy on all other time points
  ## i.e. the data of one slice is the training set while the data of all other slices are 
  ## sequentially the test sets
  ## to prevent overfitting, GAT is nested within CV loop
  #INPUT ---
  #data: data frame or list of trials
  #slice: if True, slice trial data according to window/overlap,
  #       else data is expected to be sliced. Otherwise GAT is meaningless.
  #CSP: if True, apply CSP per slice and use CSP features to decode
  #     note: it is recommended to band-pass filter the data for CSP
  if ( !slice & !.is.sliced(data) ) stop( "Data must be sliced for GAT to be meaningful." )
  .eval_package("caret")
  args.slice = .eval_ellipsis("data.trials.slice", ...)
  args.slice$split = F #save RAM
  k = .eval_ellipsis("data.folds.split", ...)$k #get k for the fold split
  args.fit = list(...)
  if (CSP) {
    .loadFuncs(CSP=T)
    if ( is.null( args.fit$method) ) {
      method = "CSP.apply"
    } else {
      method = ifelse( grepl("spec", tolower(args.fit$method)), "SpecCSP.apply", "CSP.apply" )
    }    
    args.csp = modifyList( .eval_ellipsis(method, ...), list(baseline=NULL) )
    args.csp.feat = .eval_ellipsis("CSP.get_features", ...)
    #insert get_features arguments for CSP function
    args.csp[["..."]] = NULL; args.csp = modifyList(args.csp, args.csp.feat[3:4])
    #if features are to be log-transformed, turn off centering and scaling
    if ( eval( args.csp.feat$logtransform ) & eval( args.csp.feat$approximate ) ) {
      args.fit = modifyList( args.fit, list(scale=F) )
    }
  }  
  if ( .is.sliced(data) ) { 
    #unslice back to df
    temp = .data.unslice(data, params.out=T)
    data = temp$data
    if (!slice) { #save settings if no new slice request was made on function call
      slice = T 
      args.slice = modifyList( args.slice, list(window=temp$window, overlap=temp$overlap) )
    }
  }
  #get indices for slices
  sx = do.call(data.trials.slice, modifyList( args.slice, list(data=data) ))
  slicenums = unique(sx$splitidx) #number of slices
  data = data.check(data, aslist=T, strip=F) #transform to list if necessary
  outcome = as.factor( sapply(data, "[[", 1, 2) ) #1st value of every trial
  if ( CSP && length( levels(outcome) ) > 2 ) { #more than 2 classes
    print( "More than 2 classes. Switching to OVR CSP to resolve multiclass problem." )
    method = "CSP.Multiclass.apply"
  }
  #start GAT procedure
  result = setNames( replicate(numcv, {
    ## To reduce RAM demands, grab data only as needed ##
    indxFolds = caret::createFolds(y=outcome, k=k) #get k-fold indices
    foldFits = setNames( lapply(1:k, function(foldnum) { #iterate folds
      ## part 1: train slice models
      slice.fits = lapply( slicenums, function(i) {
        #grab slice data sequentially to go easy on the RAM (relevant if parallelized)
        #slice data is 100/(k-1)% of the trial data at the time point to train the model
        sliceidx = sx$slides[sx$splitidx == i][1:(args.slice$window)]
        slice = lapply( data, function(d) na.omit( d[ sliceidx, ] ) ) #imbalanced slices generate NAs
        slice.train = slice[ -indxFolds[[foldnum]] ] #train set
        slice.test = slice[ indxFolds[[foldnum]] ] #test set
        if (CSP) {
          temp = .CSP.create_sets(slice.train, slice.test, args.csp, args.csp.feat, method)
          slice.train = temp$train
          slice.test = temp$test
        }
        #get the model fit
        fit = do.call(data.fit_model, modifyList( args.fit, list(trainData=slice.train, testData=slice.test) ))$fit
        if (CSP) {
          fit$CSP = temp$CSP[!grepl("features|outcome", names(temp$CSP))]
        }
        fit #return to lapply (slice.fits)
      })
      ## part 2: omnibus prediction of slices with model fits
      GAT = setNames( lapply(slice.fits, function(fit) { #iterate the fits
        setNames( lapply( slicenums, function(i) { #iterate the slices
          sliceidx = sx$slides[sx$splitidx == i][1:(args.slice$window)]
          testData = as.data.frame( data.table::rbindlist( #to df
            lapply( data[ indxFolds[[foldnum]] ], function(d) na.omit( d[ sliceidx, ] ) )))
          outcome.test = testData[,2] #1 value per sample
          if (CSP) {
            #transform test set with CSP filters of training set
            testData = do.call(CSP.get_features, modifyList( args.csp.feat, list(data=testData, filters=fit$CSP$filters) ))
            if ( eval(args.csp.feat$approximate) ) {
              outcome.test = outcome[ indxFolds[[foldnum]] ] #1 value per trial
            }
          } 
          if ( "scale" %in% names(fit) ) {
            #scale test set with scaling parameters of training set
            testData = scale(testData[, is.datacol(testData)], fit$scale$center, fit$scale$scale)
          }
          if ( !is.factor(outcome.test) ) outcome.test = as.factor(outcome.test)
          preds = predict( fit, testData[, is.datacol(testData)] )
          cm = caret::confusionMatrix( preds[[1]], outcome.test )
          acc = AUC::auc( AUC::roc( preds[[1]], outcome.test ) )
          list(AUC=acc, CM=cm, fit=fit) #the fits are added redundantly
        }), paste0("slice",slicenums) )
      }), paste0("fit",slicenums) )
      #prepare output to foldFits
      out = list( GATmatrix = sapply(GAT, sapply, "[[", "AUC") ) #columns: fit, rows: generalization
      if (verbose) {
        #grab for every fit the outcome (AUC, CM) at each slice and take the fit only once
        out$fits = lapply(GAT, function(fit) {
          res = lapply(fit, "[", -3)
          list( fit = fit[[1]]$fit, performance = res )
          }) #grab AUC and CM, and only 1st fit entry (rest redundant)
      }
      out #return to lapply (foldFits)
    }), paste0("Fold", 1:k) )
    #prepare output to replicate
    GATs = lapply(foldFits, "[[", "GATmatrix")
    out = list( GAT = Reduce("+", GATs)/k, GATmatrices=GATs ) #averaged over folds
    if (verbose) out$fits = lapply(foldFits, "[[", "fits")
    out #return to replicate (result)
  }, simplify=F), paste0("CV", 1:numcv) )
  #preapre output
  out = list( GAT = Reduce( "+", lapply(result, "[[", "GAT") )/numcv ) #averaged for numcv
  out$GATmatrices = lapply(result, "[[", "GATmatrices")
  if (verbose) out$fits = lapply(result, "[[", "fits")
  return(out)
}

decoding.GAT.parallel <- function(data, nCores=NULL, ...) {
  ## parallelized wrapper of decoding.GAT
  #INPUT ---
  #data: list of data frames where each df corresponds to a single subject
  #      note: make sure 1st and 2nd column in these dfs is samplenum/outcome
  #other inputs as with decoding.GAT
  #RETURNS ---
  #a list with first element as the GAT matrix averaged over all subjects, 
  #and the other elements containing the individual GAT matrices
  if ( class(data) != "list" | class(data[[1]]) != "data.frame" ) {
    stop( "Please supply data of all subjects in the form of a list with data frames." )
  }
  .eval_package("foreach", load=T)
  ## initiate parallel backend
  if ( is.null(nCores) ) {
    #initialize all available cores but not more than needed
    nCores = min( parallel::detectCores(), length(data) )
  }
  CPUcluster = .parallel_backend(on=T, nCores=nCores)
  subs = seq_along(data)
  
  #TEMP while not a package: (later use .packages argument)
  wrappedGAT <- function(...) {
    source("myFuncs.R")
    .loadFuncs(util=T, ML=T)
    return( do.call(decoding.GAT, list(...)) )
  }
  args.in = list(...)
  #execute fitting procedure
  GATres = foreach(sub = subs, .combine = list, .multicombine = T) %dopar%
    do.call(wrappedGAT, utils::modifyList( args.in, list(data=data[[sub]]) ))
  
  ## close parallel backend
  .parallel_backend(on=F, CPUcluster)
  #prepare and return output
  names(GATres) = paste0("Subject", subs)
  GATavg = Reduce("+", lapply(GATres, "[[", "GAT"))/length(subs)
  return( append(GATres, list(GATAll=GATavg), after=0 ) )
}

decoding.LOSO <- function(data, slice=T, permute=T, nCores=1, verbose=F, ...) {
  ## leave one subject out decoding procedure 
  ## splits data into training and test set where 
  ##  - training set: sum subject data sets - 1
  ##  - test set: 1 subject data set (left out)
  ## i.e. no randomness in the assignment of folds
  ## can be used to detect generalizability of information content in the signals
  ## or to identify subjects that deviate strongly from the population
  #INPUT ---
  #data: list of data frames where each df corresponds to a single subject
  #      note: make sure 1st and 2nd column in these dfs is samplenum/outcome
  #slice: if True, slice trial data according to window/overlap
  #       important: if data is not sliced, one makes the strong and most likely
  #       incorrect assumption that samples are independent of time, i.e. 
  #       models will be fitted with data from beginning and end of trials simultaneously
  #       essentially, not slicing will only be appropriate if there's 1 measurement per trial
  #permute: if True, does both prediction on true and shuffled labels
  #         and evaluates significance  
  #nCores: if larger 1, uses foreach package to parallelize computation
  #verbose: if True, extensive output which includes the fitted model and CM
  #RETURNS ---
  #if verbose = F:
  #a summary df with columns: left-out subject, AUC and slice if applicable
  #if verbose = T:
  #a list with 2 elements: summary (as before) and models with elements: 
  #subjects (+ slices if applicable), each with the output of a single call to data.fit_model
  if ( class(data) != "list" | class(data[[1]]) != "data.frame" ) {
    stop( "Please supply data of all subjects in the form of a list with data frames." )
  }
  .eval_package("foreach", load=T)
  if (nCores > 1) {
    CPUcluster = .parallel_backend(on=T, nCores=nCores)
  } else {
    registerDoSEQ() #sequential backend for foreach
  }
  subs = 1:length(data)
  #slice up all subject data sets
  if (slice) { 
    args.in = .eval_ellipsis("data.trials.slice", ...)
    data = lapply(subs, function(sub) {
      do.call(data.trials.slice, modifyList( args.in, list(data=data[[sub]]) ))  
    })  #data: list of subjects, list of slices
  }
  #create internal function to run in foreach environment:
  #each function execution represents one LOSO step
  LOSO.run <- function(sub, slice, permute, verbose, ...) {
    source("myFuncs.R") # TEMP
    .loadFuncs(util=T, ML=T) # TEMP
    args.fit = list(...)
    ## NON-SLICED CASE:
    if (!slice) {
      #simply create data frames
      trainData = as.data.frame( data.table::rbindlist(data[-sub]) )
      testData = as.data.frame( data.table::rbindlist(data[sub]) )
      #fit model with true labels
      truefit = do.call(data.fit_model, 
                        utils::modifyList( args.fit, list(trainData=trainData, 
                                                          testData=testData) ))
      result = data.frame(label="true", AUC=truefit$AUC)
      if (permute) {
        #can be left at independent = F as the general case
        permTrain = data.permute_labels(trainData, shuffle=T)
        permTest = data.permute_labels(testData, shuffle=T)
        #fit model with random labels
        randfit = do.call(data.fit_model, 
                          utils::modifyList( args.fit, list(trainData=permTrain$data, 
                                                            testData=permTest$data) ))
        result = rbind( result, data.frame(label="random", AUC=randfit$AUC) )
      }
    ## SLICED CASE:
    } else {
      numslices = 1:length(data[[1]])
      trainData = lapply(numslices, function(i) {
        #pull the i-th slice out of every subj data and bind to the others
        as.data.frame( data.table::rbindlist( lapply(data[-sub], "[[", i) ) )
      })
      testData = data[[sub]]
      #fit models with true labels for every slice
      truefit = lapply(numslices, function(slice) {
        do.call(data.fit_model, 
                utils::modifyList( args.fit, list(trainData=trainData[[slice]], 
                                                  testData=testData[[slice]]) ))
      })
      result = data.frame( label="true", slice=numslices, 
                           AUC=sapply(truefit, "[[", "AUC") )
      names(truefit) = paste0("slice",numslices) #in case verbose=T
      if (permute) {
        #can be left at independent = F as the general case
        permTrain = data.permute_labels(trainData, shuffle=T)
        permTest = data.permute_labels(testData, shuffle=T)
        #fit models with random labels for every slice
        randfit = lapply(numslices, function(slice) {
          do.call(data.fit_model, 
                  utils::modifyList( args.fit, list(trainData=permTrain$data[[slice]], 
                                                    testData=permTest$data[[slice]]) ))
        })
        result = rbind(result, data.frame( label="random", slice=numslices, 
                                           AUC=sapply(randfit, "[[", "AUC") ))
      }
    }
    #create output
    out = list(summary = result)
    if (verbose) {
      out$true = truefit
      if (permute) {
        out$random = randfit
      }
    }
    return(out)
  }
  args.LOSO.run = list(slice=slice, permute=permute, verbose=verbose)
  args.LOSO.run = utils::modifyList(args.LOSO.run, list(...))
  #execute fitting procedure
  LOSOres = foreach(sub = subs, .combine = list, .multicombine = T) %dopar%
    do.call(LOSO.run, utils::modifyList(args.LOSO.run, list(sub=sub)))
    #LOSO.run(sub, slice, permute, verbose, args.fit)
  
  if (nCores > 1) {
    .parallel_backend(on=F, CPUcluster)
  }
  out = .output.prepare(LOSOres, "subject")
  #indicate that sub was left out with negative number
  out$summary$subject = -1*out$summary$subject
  #test for significance
  if (permute) {
    args.in = .eval_ellipsis("decoding.signtest", ...)
    if ( !slice ) {
      #simple dependent t.test
      args.in = modifyList( args.in, list(group="") )
    } 
    out$significance = do.call(decoding.signtest, modifyList( args.in, list(result=out$summary) ))
  }
  return(out)
}

decoding.SS.parallel <- function(data, nCores=NULL, ...) {
  ## parallelized wrapper of decoding.SS
  #INPUT ---
  #data: list of data frames where each df corresponds to a single subject
  #      note: make sure 1st and 2nd column in these dfs is samplenum/outcome
  #other inputs as with decoding.SS
  if ( class(data) != "list" ) {
    stop( "Please supply the data of your subjects as a list, ", 
          "where each element corresponds to a single subject." )
  }
  .eval_package("foreach", load=T)
  ## initiate parallel backend
  if ( is.null(nCores) ) {
    #initialize all available cores but not more than needed
    nCores = min( parallel::detectCores(), length(data) )
  }
  CPUcluster = .parallel_backend(on=T, nCores=nCores)
  subs = 1:length(data)
  
  #TEMP while not a package: (later use .packages argument)
  wrappedSS <- function(...) {
    source("myFuncs.R")
    .loadFuncs(util=T, ML=T, CSP=T)
    return( do.call(decoding.SS, list(...)) )
  }
  args.in = list(...)
  #execute fitting procedure
  SSres = foreach(sub = subs, .combine = list, .multicombine = T) %dopar%
    do.call(wrappedSS, utils::modifyList( args.in, list(data=data[[sub]]) ))

  ## close parallel backend
  .parallel_backend(on=F, CPUcluster)
  #prepare and return output
  out = .output.prepare(SSres, "subject")
  #if SpecCSP was used, average the alphas for all subjects
  if ( "fits.SpecCSP" %in% names(out) ) {
    temp = list( true = lapply( out$fits.SpecCSP, "[[", "alpha.true" ) ) #true alphas for all subs
    if ( "alpha.random" %in% names(out$fits.SpecCSP[[1]]) ) {
      temp$rand = lapply( out$fits.SpecCSP, "[[", "alpha.random" ) #random alphas for all subs
    }
    sliced = class( temp$true[[1]] ) == "list" #slice list
    specout = setNames( lapply( temp, function(alpha) {
      if (sliced) {
        slices = names( alpha[[1]] ) #sub1
        setNames( lapply(slices, function(slice) {
          Reduce("+", lapply(alpha, "[[", slice))/length(alpha)
        }), slices)
      } else {
        Reduce("+", alpha)/length(alpha)
      }
    }), paste0( "alpha.", names(temp)) )
    specout$settings = out$fits.SpecCSP[[1]]$settings
    out$fits.SpecCSP = append(out$fits.SpecCSP, list(overall=specout), after=0 )
  }
  return(out)
}

.output.prepare <- function(output, newStr = "subject") {
  ## helper to prepare output for user accessibility
  idx = 1:length(output)
  elements = names(output[[1]])
  elem_classes = sapply(output[[1]], class)
  out = list()
  for (elem in elements) {
    if ( elem_classes[elem] == "data.frame" ) {
      #merge to one df
      temp = plyr::rbind.fill( lapply(output, "[[", elem) )
      #prevent error in case of varying slice lengths:
      nrows = sapply( lapply(output, "[[", elem), nrow )
      #append newStr column, e.g. subject number
      temp = data.frame(rep( idx, nrows ), temp)
      names(temp)[1] = newStr
      out[[elem]] = temp #append to list
    } else { #list
      #grab for every newStr (e.g. subject) the corresponding element
      temp = lapply(output, "[[", elem)
      #rename it
      names(temp) = paste0(newStr, idx)
      idxname = ifelse( grepl("fits\\.", elem), "", "fits." )
      out[[paste0(idxname,elem)]] = temp #append to list
    }
  }
  return(out)
}

decoding.SS <- function(data, slice=T, permute=T, numcv=1, verbose=F, CSP=F, ...) {
  ## single subject decoding procedure with repeated k-fold CV
  ## can be easily parallelized with foreach for all subjects
  #INPUT ---
  #data: continuous df, list of trials, or sliced data (with slice=F)
  #slice: if True, data will be sliced. Set window and overlap accordingly.
  #       otherwise, data is left as is.
  #permute: if True, does both prediction on true and shuffled labels
  #         and evaluates significance
  #numcv: number of repetitions of the k-fold CV
  #verbose: if True, extensive output with all model fits is provided
  #         otherwise, only important (averaged) results are presented
  #CSP: if True, decode with logvar of CSP components per slice (must be True)
  #     note: it is recommended to band-pass filter the data for CSP
  #RETURNS ---
  #if verbose & permute = F: a summary data frame with prediction accuracy per 
  #                          repetition (numcv), fold, and slice if applicable
  #if permute = T: a list with summary df and a df with significance testing result
  #if verbose = T: a list with 2 or 4 (if permute = T) elements:
  #                1. the same summary df as with verbose = F
  #                2. a list with number of elements corresponding to numcv times
  #                   k folds, each with the model fits for true labels
  #                3. if permute = T, analogous to 2 but with random labels
  #                4. if permute = T, significance testing result df
  .eval_package("caret")
  
  ######## PREPARATION START ########
  #retrieve (...) outside of replicate
  args.slice = .eval_ellipsis("data.trials.slice", ...)
  args.slice$split = F #save RAM
  k = .eval_ellipsis("data.folds.split", ...)$k #get k for the fold split
  args.fit = list(...)
  if (CSP) {
    .loadFuncs(CSP=T)
    if ( is.null( args.fit$method) ) {
      method = "CSP.apply"
    } else {
      method = ifelse( grepl("spec", tolower(args.fit$method)), "SpecCSP.apply", "CSP.apply" )
    }    
    args.csp = .eval_ellipsis(method, ...)
    args.csp.feat = .eval_ellipsis("CSP.get_features", ...)
    #insert get_features arguments for CSP function
    args.csp[["..."]] = NULL; args.csp = modifyList(args.csp, args.csp.feat[3:4])
    #if features are to be log-transformed, turn off centering and scaling
    if ( eval( args.csp.feat$logtransform ) && eval( args.csp.feat$approximate ) ) {
      args.fit = modifyList( args.fit, list(scale=F) )
    }
  }
  #ensure legitimate values for k and numcv
  k = as.integer(k); numcv = as.integer(numcv)
  if ( numcv < 1 ) numcv = 1L
  if ( k <= 1 ) {
    k = 1L
    if ( numcv > 1 ) {
      numcv = 1L
      warning( "Repetitions are meaningless without a fold-split. Setting numcv to 1." )
    }
  }
  #evaluate data, get outcome
  Data.true = data.permute_labels(data, shuffle=F)
  if ( CSP && length( levels(Data.true$outcome ) ) > 2 ) { #more than 2 classes
    print( "More than 2 classes. Switching to OVR CSP to resolve multiclass problem." )
    method = "CSP.Multiclass.apply"
  }
  if ( .is.sliced(data) ) {
    #unslice back to df
    temp = .data.unslice(Data.true$data, params.out=T)
    Data.true$data = temp$data
    if (!slice) { #save settings if no new slice request was made on function call
      slice = T 
      args.slice = modifyList( args.slice, list(window=temp$window, overlap=temp$overlap) )
    }
  }
  Data.true$data = data.check(Data.true$data, aslist=T, strip=F) #transform to list if necessary
  #evaluation end: data is a list of trials at this point
  
  #get indices for slices if requested
  if (slice) {
    sx = do.call(data.trials.slice, modifyList( args.slice, list(data=Data.true$data) ))
    slicenums = unique(sx$splitidx) #number of slices
    #define internal slice fitting function
    slice.fitting <- function(train, test) {
      ## takes the training and test data sets (lists of trials) and fits each slice
      ## returns a list of fits with each list having the elements: AUC, CM, fit
      fits = setNames( lapply( slicenums, function(i) { #iterate slices
        #get current slices (additional memory required!) in trial list format...
        #only take the first 1:window indices to sequentially go through the trials:
        sliceidx = sx$slides[sx$splitidx == i][1:(args.slice$window)]
        slice.train = lapply( train, function(d) na.omit( d[ sliceidx, ] ) ) #imbalanced slices generate NAs
        slice.test = lapply( test, function(d) na.omit( d[ sliceidx, ] ) )
        if (CSP) {
          #transform slice data into CSP features
          csp = .CSP.create_sets(slice.train, slice.test, args.csp, args.csp.feat, method)
          slice.train = csp$train
          slice.test = csp$test            
        }
        fit = do.call(data.fit_model, modifyList( args.fit, list(trainData=slice.train, testData=slice.test) ))
        if (CSP && verbose) {
          #add CSP output to fit but without features and outcome to keep the output lean
          fit$fit$CSP = csp$CSP[!grepl("features|outcome", names(csp$CSP))]
        }
        fit # return to lapply
      }), paste0("slice",slicenums) )
      return( fits )
    }
  }
  ######## PREPARATION END ########
  
  ######## FITTING START ########
  #repeat numcv times...
  result = replicate(numcv, {
    ## To reduce RAM demands, grab data only as needed ##
    indxFolds.true = caret::createFolds(y=Data.true$outcome, k=k) #get k-fold indices (save RAM)
    #permute labels if requested:
    if (permute) { 
      #memory allocation mostly unchanged due to reference to Data.true (except for outcome column)
      Data.rand = data.permute_labels(Data.true$data, shuffle=T)
      indxFolds.rand = caret::createFolds(y=Data.rand$outcome, k=k)
    } 
    #start iteration over folds
    foldFits = setNames( lapply(1:k, function(foldnum) {
      #get training and test sets (all referenced so no additional memory required)
      train.true = Data.true$data[ -indxFolds.true[[foldnum]] ]
      test.true = Data.true$data[ indxFolds.true[[foldnum]] ]
      if (permute) {
        train.rand = Data.rand$data[ -indxFolds.rand[[foldnum]] ]
        test.rand = Data.rand$data[ indxFolds.rand[[foldnum]] ]
      }
      if ( k == 1 ) { #no CV
        train.true = test.true #training set and test set are the same
        if (permute) train.rand = test.rand
      }
      if (slice) {
        fits.true = slice.fitting(train.true, test.true)
        result.true = data.frame( label="true", slice=slicenums, AUC=sapply(fits.true, "[[", 1) )
        if (permute) {
          fits.rand = slice.fitting(train.rand, test.rand)
          result.rand = data.frame( label="random", slice=slicenums, AUC=sapply(fits.rand, "[[", 1) )
        }
      } else {
        #simple fitting procedure on trial data
        if (CSP) {
          csp.true = .CSP.create_sets(train.true, test.true, args.csp, args.csp.feat, method)
          train.true = csp.true$train; test.true = csp.true$test
          if (permute) {
            csp.rand = .CSP.create_sets(train.rand, test.rand, args.csp, args.csp.feat, method)
            train.rand = csp.rand$train; test.rand = csp.rand$test
          }
        }
        fits.true = do.call(data.fit_model, modifyList( args.fit, list(trainData=train.true, testData=test.true) ))
        result.true = data.frame(label="true", AUC=fits.true$AUC)
        if (permute) {
          fits.rand = do.call(data.fit_model, modifyList( args.fit, list(trainData=train.rand, testData=test.rand) ))
          result.rand = data.frame(label="rand", AUC=fits.rand$AUC)
        }
        if ( verbose && CSP ) {
          #add the important CSP output to fit
          fits.true$fit$CSP = csp.true$CSP[!grepl("features|outcome", names(csp.true$CSP))]
          if (permute) fits.rand$fit$CSP = csp.rand$CSP[!grepl("features|outcome", names(csp.rand$CSP))]
        }
      }
      out = list(result.true = result.true)
      if (verbose) out$fits.true = fits.true
      if (permute) {
        out$result.random = result.rand
        if (verbose) out$fits.random = fits.rand
      }
      out #return to lapply
    }), paste0("Fold", 1:k) )
    #prepare return to replicate...
    summary = plyr::ldply( lapply(foldFits, "[[", "result.true"), .id="fold" )
    if (permute) summary = rbind( summary, plyr::ldply( lapply(foldFits, "[[", "result.random"), .id="fold" ) )
    if (slice) {
      summary = summary[ order(summary$label, summary$slice), ]
      row.names(summary) = NULL #unsorted numbers due to ordering
    }
    out = list(summary = summary)
    if (verbose) {
      out$true = lapply(foldFits, "[[", "fits.true")
      if (permute) out$random = lapply(foldFits, "[[", "fits.random")
    }
    out #return to replicate
  }, simplify=F)
  ######## FITTING END ########
  
  #finalize output
  out = .output.prepare(result, "CV")
  #test for significance
  if (permute && k > 1) {
    args.in = .eval_ellipsis("decoding.signtest", ...)
    if ( !slice ) {
      #simple independent t.test
      args.in = modifyList( args.in, list(group="") )
    } 
    out$significance = do.call(decoding.signtest, modifyList( args.in, list(result=out$summary) ))
  }
  #if SpecCSP add aggregated alpha per slice 
  if (verbose && CSP && method=="SpecCSP.apply") {
    alphas = setNames( lapply( out[ grep("fits", names(out)) ], function(fits) { #true and random
      if (slice) {
        slices = names( out$fits.true[[1]][[1]] ) #cv1, fold1
        setNames( lapply(slices, function(slice) {
          temp = lapply(fits, function(cv) {
            Reduce("+", lapply(cv, function(fold) { fold[[slice]]$fit$CSP$alpha }))/k
          })
          Reduce("+", temp)/numcv
        }), slices)
      } else {
        temp = lapply(fits, function(cv) {
          Reduce("+", lapply(cv, function(fold) { fold$fit$CSP$alpha }))/k
        })
        Reduce("+", temp)/numcv
      }
    }), paste0( "alpha.", substr( names(out[ grep("fits", names(out)) ]), 6, 12) ))
    out$SpecCSP = alphas
    if (slice) {
      out$SpecCSP$settings = list(bands = out$fits.true[[1]][[1]][[1]]$fit$CSP$bands, 
                                  freqs = out$fits.true[[1]][[1]][[1]]$fit$CSP$freqs, 
                                  outcome = levels(Data.true$outcome))
    } else {
      out$SpecCSP$settings = list(bands = out$fits.true[[1]][[1]]$fit$CSP$bands, 
                                  freqs = out$fits.true[[1]][[1]]$fit$CSP$freqs, 
                                  outcome = levels(Data.true$outcome))      
    }
  }
  return(out)
}

decoding.signtest <- function(result, dv="AUC", pair="label", group="slice", 
                              adjust="none", alpha=0.01, parametric=T) {
  ## does one-sided t-tests at all slices for significant difference in decodability 
  ## between true and random label and adjusts p-values for multiple NHT
  #INPUT ---
  #result: df with columns specified by dv, pair and group
  #        or a list of such dfs corresponding to different subjects
  #dv: dependent variable, column with the classification accuracy
  #pair: pair to compare, column with the labels (true vs. random)
  #group: variable to split the data for, e.g. column with the slices
  #       if empty char, essentially simple independent t-test without split
  #adjust: correction method for multiple NHT, default is "none"
  #alpha: significance threshold
  #RETURNS ---
  #a df with slice number, mean difference (true - random), 
  #p value (with applied correction method)
  #and whether it is significant (1) or not (0), i.e. < alpha 
  #NOTE ---
  #does paired samples t-tests per slice if subject list is supplied, otherwise
  #tests for independent samples, i.e. true vs random performance at each slice
  if ( group == "" ) {
    #placeholder for simple independent t.testing (aggregate will do nothing)
    result = cbind(.id=1, result)
    group = ".id"
  }
  paired = F
  form = as.formula( paste(dv, pair, sep="~") )
  if ( class(result) == "list" | "subject" %in% names(result) ) { #subject list
    if ( class(result) == "list" & !"subject" %in% names(result[[1]]) ) {
      #list of dfs withiout subject info
      subject = rep( 1:length(result), each=nrow(result[[1]]) )
      result = cbind( subject, plyr::rbind.fill(result) )
    } else {
      #list of dfs or single df with subject info
      result = plyr::rbind.fill(result)
    }
    paired = T
    #formula to aggregate for subjects
    aggform = as.formula( paste(dv, paste(pair, "subject", sep="+"), sep="~") )
  }
  #split data
  slices = split(result, result[,group])
  #test every slice
  if (parametric) { #do t.tests
    tests = sapply(slices, function(slice) {
      if (paired) { #average for subjects
        slice = aggregate(aggform, data=slice, mean)
      }
      test = t.test(form, data=slice, paired=paired, alternative="greater")
      diff = ifelse(paired, unname( test$estimate[1] ), 
                    unname( test$estimate[1]-test$estimate[2] ))
      c(diff=diff, pval=test$p.value)
    })
    #correct for multiple NHT
    pvals = p.adjust(tests[2,], method=adjust)
    resultdf = as.data.frame( cbind( 1:length(pvals), meandiff=tests[1,], 
                                     pval=pvals, significant.p=(pvals<alpha) ) )
  } else { #do bootstrap
    #internal bootstrap function
    .bootfun <- function(x, iter=10000) { #10000 iterations
      return( replicate(iter, {
        mean( sample(x, replace=T) )
      }) )
    }
    #pair level to bootstrap
    random = ifelse( "random" %in% unique( result[[pair]] ), "random", 
                     unique( result[[pair]] )[2] )
    tests = sapply(slices, function(slice) {
      boot = .bootfun( slice[[dv]][ slice[[pair]] == random ] )
      CI = quantile(boot, probs=1-alpha)
      truemean = mean( slice[[dv]][ slice[[pair]] != random ] )
      diff = truemean - mean(boot)
      pval = sum( boot >= truemean ) / length(boot)
      c(true=truemean, diff=diff, pval=pval, ci=CI)
    })
    pvals = p.adjust(tests[3,], method=adjust)
    resultdf = as.data.frame( cbind( seq_along(slices), AUC=tests[1,],
                                     meandiff=tests[2,], 
                                     pval=pvals, ci=tests[4,],
                                     significant.p=pvals<alpha,
                                     significant.obs=tests[1,]>tests[4,] ) )
  }
  names(resultdf)[1] = group
  if (group == ".id") { resultdf = resultdf[,-1] }
  return( resultdf )
}

data.permute_labels <- function(data, shuffle=T) {
  ## evaluates data format, extracts labels from data and permutes if requested
  #INPUT ---
  #data: continuous df with 1st col: sample number, 2nd col: outcome 
  #      or list of trials; or list of slices
  #      note: a single outcome per trial is assumed
  #shuffle: if True, labels are shuffled (permutation tests) and overwritten
  #NOTES ---
  #case 1: no slices, 1 sample per trial (single df)
  #        - outcome is the 2nd column
  #        - will be simply overwritten with shuffled labels
  #case 2: no slices, multi sample per trial (df/trial list)
  #        - outcome is the 1st value from every trial
  #        - will loop over trials and replace with new labels, each repeated for nrow 
  #        - note: case 1 can be considered a special case with nrow = 1
  #case 3: slices, multi sample per trial (slice list)
  #        - outcome is extracted from 1st slice (identical for all others)
  #          with 1 value per trial (inferred from window size / sample number)
  #        - will loop over slices and replace with new labels, each repeated for samplenum
  #        - note: case 4 can be considered a special case with samplenum = 1
  #case 4: slices, single sample per trial (slice list)
  #        - outcome is the 2nd column of the 1st slice (identical for all others)
  #        - will loop over slices and replace with new labels
  #RETURNS ---
  #list with 3 elements: 
  # - data in case format 
  # - outcome vector (1 value per trial)
  if ( class(data) != "list" ) {
    # multi-sample data with dependency?
    independent = ( length( unique(data[,1]) ) == nrow(data) )
    sliced = FALSE
    if (independent) { #case 1: df with 1 sample per trial
      outcome = data[,2] 
    } else { #case 2a: df with multiple samples per trial
      nsamples = data.get_samplenum(data)
      outcome = data[ seq(1, nrow(data), nsamples), 2]
      data = data.trials.split(data, strip=F) #transform to trial list
    }
  } else { # a list
    nsamples = data.get_samplenum(data[[1]])
    independent = nsamples == 1 #TRUE if 1 sample per trial
    sliced = .is.sliced(data)
    if (sliced & !independent) { #case 3: slice list with multiple samples per trial
      outcome = data[[1]][seq( 1, nrow(data[[1]]), by=nsamples ), 2] #same order in every slice
    } else if (sliced & independent) { #case 4: slice list with 1 sample per trial
      outcome = data[[1]][,2] #same order in every slice
    } else { #not sliced
      #case 2b (or 1b if independent): as 2a (or 1) but a trial list
      outcome = sapply(data, "[[", 1, 2) #first sample of every trial df
    }
  }
  if (shuffle) { 
    outcome = sample(outcome) #shuffle outcome
    #change labels within the data
    if (!sliced) {
      if (class(data) != "list") { #case 1: df
        data[,2] = outcome #simply replace 2nd column
      } else { #case 2 (or case 1b): trial list
        data = lapply(1:length(data), function(i) {
          data[[i]][,2] = outcome[i]
          data[[i]]
        })
      }
    } else { #sliced
      if (!independent) { #case 3: slice list, multiple samples
        nsamples = sapply(data, data.get_samplenum) #sample nums per slice
        data = lapply(1:length(data), function(i) {
          data[[i]][,2] = rep(outcome, each=nsamples[i]) #repeat outcome for window size
          data[[i]]
        })
      } else { #case 4: slice list, single sample
        data = lapply(data, function(slice) {
          slice[,2] = outcome #1 value per trial
          slice
        })
      }
    }
  }
  return( list(data = data, outcome = as.factor(outcome)) )
}

data.train_test.split <- function(data, ratio=0.7, split=T) {
  ## splits data into training and test data set according to the ratio
  ## the split is balanced for outcome and stratified for trials
  #INPUT ---
  #data: df, list of trials/slices
  #ratio: percentage of data that goes into training
  #split: if True, data is returned as train and test set, 
  #       otherwise the the trial indices of the training set are returned
  #RETURNS ---
  #either a numeric vector with the trial indices (if split=F) or a list
  #with train/test set. If data was sliced, the output list will be a slice
  #list with train and test elements on each slice
  .eval_package("caret")
  temp = data.permute_labels(data, shuffle=F)
  trainidx = caret::createDataPartition(y=temp$outcome, p=ratio)[[1]]
  if (!split) return(trainidx)
  #split data into train and test set
  if ( .is.sliced(temp$data) ) { #data is a slice list
    #iterate the slices and split
    splitlist = setNames( lapply( seq_along(temp$data), function(slice) {
      #get slice data into trial format
      slice.data = data.trials.split( temp$data[[slice]] )
      list( train = data.check( slice.data[ trainidx ], aslist=F ), 
            test = data.check( slice.data[ -trainidx ], aslist=F ) )
    }), paste0("slice", seq_along(temp$data)) )
  } else if ( class(temp$data) == "data.frame" ) { #data is a df
    splitlist = list(train = temp$data[ trainidx, ], test = temp$data[ -trainidx, ])
  } else { #data is a trial list
    toDF = class(data) == "data.frame" #if input was a df, transform split back to df
    splitlist = list(train = data.check( temp$data[ trainidx ], aslist=!toDF ),
                     test = data.check( temp$data[ -trainidx ], aslist=!toDF ))
  }
  return( splitlist )
}

data.folds.split <- function(data, k=5, shuffle=F) {
  ## fold data for k-fold cross validation
  #INPUT ---
  #data: continuous df with 1st col: sample number, 2nd col: outcome 
  #      or list of trials; or list of slices
  #      note: a single outcome per trial is assumed
  #k: number of folds
  #independent + shuffle see data.permute_labels
  #if shuffle is True, labels in the data are changed!
  #RETURNS ---
  #list of k elements with sublists comprising:
  #  - training set (corresponding to k-1 folds)
  #  - test set (corresponding to 1 fold)
  #if data is sliced, the structure is {k, num slices, training&test set}
  #NOTES ---
  #case 1: single measurement df = every row is considered a trial
  #        the resulting output will be simply a split by outcome into training and test set
  #        -> behavioral experiment scenario
  #case 2: repeated measures, df or trial list = every list element is considered a trial
  #        the resulting output will be a split by outcome and merge the correponding samples
  #        of a trial (all of them) into a training and test set, respectively
  #        i.e. the training set and test set will comprise data of entirely different trials
  #        -> the probable choice for paradigms such as CSP or SpecCSP
  #case 3: repeated measures, slice list = every list element is an independent data set
  #        for each data set, the outcome of every trial is represented by a fixed number of samples
  #        the output (within each fold) will be a list corresponding to the number of slices
  #        with each slice containing a training and test set analogous to case 2
  #        i.e. if a trial is represented by 20 samples in a slice, all of these will be
  #        assigned either to the training or the test set!
  #        -> most common choice for decoding analyses
  #case 4: single measurement, slice list = every list element is an independent data set
  #        otherwise resembles case 1
  #        -> probable in decoding analyses with low sampling rates
  .eval_package("caret")
  #evaluate data and permute labels if requested:
  if ( class(data) != "list" ) {
    independent = ( length( unique(data[,1]) ) == nrow(data) )
    sliced = FALSE
  } else { # a list
    independent = data.get_samplenum(data[[1]]) == 1 #TRUE if 1 sample per trial
    sliced = .is.sliced(data)
    if (sliced & !independent) { #unslice if necessary to speed up folding
      temp = .data.unslice(data, params.out = T)
      window = temp$window; overlap = temp$overlap
      data = data.trials.split(temp$data, strip=F)
    }
  }
  #create fold indices based on outcome:
  temp = data.permute_labels(data, shuffle=shuffle)
  data = temp$data #if shuffle is TRUE, labels are different now
  indxFolds <- caret::createFolds(y=temp$outcome, k=k)
  #create data folds:
  foldedData = lapply(indxFolds, function(fold) {
    if (!sliced) {
      if (independent) { #case 1: df
        trainData = data[-fold,]
        testData = data[fold,]
        foldData = list(train=trainData, test=testData) #2 dfs
      } else { #case 2: trial list or df
        
        trainData = data[-fold]
        testData = data[fold]
        foldData = list(train=trainData, test=testData) #2 lists with dfs
      }
    } else { #sliced
      if (!independent) { #case 3: trial list (to be sliced)
        trainData = suppressWarnings( data.trials.slice(data[-fold], window, overlap) )
        testData = suppressWarnings( data.trials.slice(data[fold], window, overlap) )
        foldData = list(train=trainData, test=testData)
      } else { #case 4: slice list
        foldData = list(train = lapply( data, function(slice) slice[-fold,] ),
                        test = lapply( data, function(slice) slice[fold,] ))
        #2 lists with slice lists in either case
      }
    }
  })
  return( foldedData )
}

data.folds.fit <- function(foldData, ...) {
  ## calls data.fit_model to fit all the data folds produced by data.folds.split
  #INPUT ---
  #foldData: a list obtained from data.folds.split, ore more generally
  #          a list with its main structure corresponding to the folds,
  #          substructure [+ optionally a list of slices] a list of trainig & test set
  #RETURNS ---
  #a list with 2 elements: summary and foldFits
  #summary: df with the AUC value per fold (and slice if applicable)
  #foldFits: list with elements corresponding to the number of folds
  #          within the respective fold are again 2 elements: 
  #          result: already integrated in the summary
  #          fit: relevant information about the model, see output of data.fit_model
  #NOTES ---
  #generally, expected input is the output from data.folds.split
  #the input to data is assumed to be a list of independent folds;
  #expected structure: {folds, train & test set} or {folds, train & test set, slices}
  if ( class(foldData) != "list" | class(foldData[[1]]) != "list" ) {
    stop( "Expected data structure is a list of folds with sublists.\n", 
          "See output of data.folds.split or make use of data.fit_model directly." )
  }
  args.in = list(...)
  foldFits = lapply(foldData, function(fold) { #iterate over folds
    if ( .is.sliced( fold$train ) ) { #list contains slices
      fit = lapply( 1:length(fold$train), function(slice) {
        do.call(data.fit_model, modifyList( args.in, list(trainData=fold$train[[slice]], 
                                                          testData=fold$test[[slice]]) ))
      })
      result = data.frame( slice=1:length(fit), AUC=sapply(fit, "[[", 1) )
      
    } else { #list contains training and test set
      fit = do.call(data.fit_model, modifyList( args.in, list(trainData=fold$train, 
                                                              testData=fold$test) ))
      result = data.frame(AUC=fit$AUC)
    }
    list(result=result, fit=fit)
  })
  #prepare output
  summary = plyr::ldply( lapply(foldFits, "[[", "result"), .id="fold" )
  if ( "slice" %in% names(summary) ) {
    fits = lapply(foldFits, function(fold) {
      fit = setNames( lapply( 1:length(fold$fit), function(slice) fold$fit[[slice]] ),
                      paste0( "slice", 1:length(fold$fit) ))
    }) 
  } else {
    fits = lapply(foldFits, "[[", "fit")
  }
  return( list(summary = summary, fits = fits) )
}

data.fit_model <- function(trainData, testData, model="L2", scale=T, ...) {
  ## machine learning procedure: fits a model to a training set 
  ## and uses the model to predict the outcome of a test set
  #INPUT ---
  #trainData: data set containing training samples, outcome in 2nd col
  #testData: data set containing test samples, outcome in 2nd col
  #model: model to fit, defaults to L2 logistic regression, optional "LDA"
  #scale: if True, data is scaled during training, and the obtained parameters are used
  #       to also scale the test set
  #RETURNS ---
  #list with elements:
  #AUC: area under the curve, prediction accuracy of the model
  #CM: confusion matrix of predictions and true values in testData
  #fit: model fit with the respective list structure
  #currently implemented models:
  models = c("L2", "LDA") 
  if ( !toupper(model) %in% models ) {
    stop( "Undefined model selected. Current options are L2 or LDA." ) 
  }
  .eval_package( c("caret", "AUC") )
  args.in = list(...)
  trainData = data.check(trainData, aslist=F)
  testData = data.check(testData, aslist=F)
  infoidx = !is.datacol(trainData) #non-data columns
  if ( !is.factor(trainData[,2]) ) {
    trainData[,2] = as.factor(trainData[,2])
    testData[,2] = as.factor(testData[,2])
  }
  if (scale) {
    trainD = scale(trainData[, !infoidx ], center=T, scale=T)
    testD = scale(testData[, !infoidx ], attr(trainD,"scaled:center"), attr(trainD,"scaled:scale"))
    trainData = cbind( trainData[, infoidx ], trainD)
    testData = cbind( testData[, infoidx ], testD)
  }
  if ( toupper(model) == "L2" ) {
    .eval_package("LiblineaR")
    type = 7 #dual L2
    cost = 1 #cost parameter (inverse of regularization parameter)
    # check optional inputs
    if ( "type" %in% names(args.in) ) type = args.in$type
    if ( "cost" %in% names(args.in) ) cost = args.in$cost
    fit = LiblineaR::LiblineaR(data=trainData[, !infoidx ], 
                               target=trainData[,2], type=type, cost=cost) #dual L2
    fit$cost = cost
  } else if ( toupper(model) == "LDA" ) {
    .eval_package("MASS")
    fit = MASS::lda(x=trainData[, !infoidx ], grouping=trainData[,2])    
  }
  if (scale) {
    fit$scale = list( center=attr(trainD,"scaled:center"), 
                      scale=attr(trainD,"scaled:scale") )
  }
  #get predictions
  preds <- predict( fit, testData[, !infoidx ] )
  cm = caret::confusionMatrix( preds[[1]], testData[,2] )
  acc = AUC::auc( AUC::roc( preds[[1]], testData[,2] ) )
  return( list(AUC=acc, CM=cm, fit=fit) )
}

createSubjFolds <- function(subnum, k) {
  ## helper to create folds for a list of subject data sets
  ## caret's createFolds does not properly handle this scenario
  ## e.g. repeat createFolds(1:20, k=5) a few times
  #INPUT ---
  #subnum: number of subjects or vector with sub numbers
  #k: number of folds
  #RETURNS: the fold indices
  subs = unique(subnum)
  if (length(subnum) == 1) subs = 1:subnum
  subs = which( subs == subs ) #get indices in case numeric vector was supplied
  cuts = sample( cut(subs, breaks=k, labels=F, include.lowest=T) )
  folds = setNames( split(subs, cuts), paste0("Fold",1:k) )
  return(folds)
}
