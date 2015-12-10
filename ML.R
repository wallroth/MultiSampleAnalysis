## Machine Learning functions

#make helper functions accessible
source("myFuncs.R")

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
  #CSP: if True, decode with logvar of CSP components per slice
  #     note: it is recommended to band-pass filter the data for CSP
  if ( !slice & !.is.sliced(data) ) stop( "Data must be sliced for GAT to be meaningful." )
  .eval_package("caret")
  args.slice = .eval_ellipsis("data.trials.slice", ...)
  args.slice$split = F #save RAM
  k = .eval_ellipsis("data.folds.split", ...)$k #get k for the fold split
  args.fit = list(...)
  if (CSP) {
    .loadFuncs(CSP=T)
    args.csp = modifyList( .eval_ellipsis("CSP.apply", ...), list(baseline=NULL) )
    args.csp.feat = .eval_ellipsis("CSP.get_features", ...)
    #if features are to be log-transformed, turn off centering and scaling
    if ( eval(args.csp$logtransform) ) args.fit = modifyList( args.fit, list(scale=F) )
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
  outcome = sapply(data, "[[", 1, 2) #1st value of every trial
  dataidx = is.datacol(data[[1]])
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
          csp = do.call(CSP.apply, modifyList( args.csp, list(data=slice.train) ))
          slice.train = data.frame(sample=i, outcome=csp$outcome, csp$features)
          slice.test = data.frame( sample=i, outcome=outcome[ indxFolds[[foldnum]] ],
            #project with training set filters
            do.call(CSP.get_features, modifyList( args.csp.feat, list(data=slice.test, filters=csp$filters) )) )
        }
        #get the model fit
        fit = do.call(data.fit_model, modifyList( args.fit, list(trainData=slice.train, testData=slice.test) ))$fit
        if (CSP) {
          fit$csp = csp$filters #add training set CSP filters to fit list
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
            testData = do.call(CSP.get_features, modifyList( args.csp.feat, list(data=testData, filters=fit$csp) ))
            outcome.test = outcome[ indxFolds[[foldnum]] ] #1 value per trial
          } else if ( "scale" %in% names(fit) ) {
            #scale test set with scaling parameters of training set
            testData = scale(testData[, dataidx], fit$scale$center, fit$scale$scale)
          } else { # no scaling
            testData = testData[, dataidx]
          }
          preds = predict( fit, testData )
          if ( !is.factor(outcome.test) ) {
            outcome.test = as.factor(outcome.test)
          }
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
  return( .output.prepare(SSres, "subject") )
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
    args.csp = modifyList( .eval_ellipsis("CSP.apply", ...), list(baseline=NULL) )
    args.csp.feat = .eval_ellipsis("CSP.get_features", ...)
    #if features are to be log-transformed, turn off centering and scaling
    if ( eval(args.csp$logtransform) ) args.fit = modifyList( args.fit, list(scale=F) )
  }
  #evaluate data, get outcome
  Data.true = data.permute_labels(data, shuffle=F)
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
  if (CSP & !slice) {
    stop( "CSP without slicing can be done with more functionalities via pipe.CSP.", 
          "For this function, please supply sliced data or define a window." )
  }
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
          csp = do.call(CSP.apply, modifyList( args.csp, list(data=slice.train) ))
          slice.train = data.frame( sample=i, outcome=csp$outcome, csp$features )
          slice.test = data.frame( sample=i, outcome=sapply( slice.test, "[[", 1, 2 ),
                                   do.call(CSP.get_features, #project with training set filters
                                           modifyList( args.csp.feat, list(data=slice.test, filters=csp$filters) )) )
                                   
        }
        do.call(data.fit_model, modifyList( args.fit, list(trainData=slice.train, testData=slice.test) ))
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
      if (slice) {
        fits.true = slice.fitting(train.true, test.true)
        result.true = data.frame( label="true", slice=slicenums, AUC=sapply(fits.true, "[[", 1) )
        if (permute) {
          fits.rand = slice.fitting(train.rand, test.rand)
          result.rand = data.frame( label="random", slice=slicenums, AUC=sapply(fits.rand, "[[", 1) )
        }
      } else {
        #simple fitting procedure on trial data
        fits.true = do.call(data.fit_model, modifyList( args.fit, list(trainData=train.true, testData=test.true) ))
        result.true = data.frame(label="true", AUC=fits.true$AUC)
        if (permute) {
          fits.rand = do.call(data.fit_model, modifyList( args.fit, list(trainData=train.rand, testData=test.rand) ))
          result.rand = data.frame(label="rand", AUC=fits.rand$AUC)
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
  out = .output.prepare(result, "cv")
  #test for significance
  if (permute) {
    args.in = .eval_ellipsis("decoding.signtest", ...)
    if ( !slice ) {
      #simple independent t.test
      args.in = modifyList( args.in, list(group="") )
    } 
    out$significance = do.call(decoding.signtest, modifyList( args.in, list(result=out$summary) ))
  }
  return(out)
}

decoding.signtest <- function(result, dv="AUC", pair="label", group="slice", 
                          adjust="none", alpha=0.01) {
  ## does t-tests at all slices for significant difference in decodability 
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
  tests = sapply(slices, function(slice) {
    if (paired) { #average for subjects
      slice = aggregate(aggform, data=slice, mean)
    }
    test = t.test(form, data=slice, paired=paired)
    diff = ifelse(paired, unname( test$estimate[1] ), 
                  unname( test$estimate[1]-test$estimate[2] ))
    c(diff=diff, pval=test$p.value)
  })
  #correct for multiple NHT
  pvals = p.adjust(tests[2,], method=adjust)
  resultdf = as.data.frame( cbind( 1:length(pvals), meandiff=tests[1,], 
                                   pval=pvals, significant=(pvals<alpha) ) )
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
  return( list(data = data, outcome = outcome) )
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




