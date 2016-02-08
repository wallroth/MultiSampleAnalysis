## Machine Learning functions

decode.GAT <- function(data, method="within", decompose="none", repetitions=1, 
                       permutations=1, nCores=NULL, ...) {
  ## Generalization across time: how long does a process stay active, how many are there
  ## -> train model at one time point and test prediction accuracy on all other time points
  ## i.e. the data of one slice is the training set while the data of all other slices are 
  ## sequentially the test sets; to prevent overfitting, GAT is nested within CV loop
  #INPUT ---
  #data: continuous df, list of trials, slices, subjects
  #method: one of c("within","across","between")
  #        within: single subject decoding; either a data set of a single subject
  #                or a list of subject data sets that will be parallelized
  #        across: data of multiple subjects is concatenated where folds are
  #                stratified for subjects, i.e. train and test set contain
  #                trials of every subject equally
  #        between: data of multiple subjects is separated into train and test set
  #                 i.e. train and test set contain trials of different subjects
  #decompose: decomposition method to use, one of c("CSP","SpecCSP","SPoC")
  #repetitions: number of repetitions of the k-fold CV, a value >= 1
  #permutations: number of permutations on each fold, corresponds to the number
  #              of true fits, i.e. a value of 1 means as many permutations as there 
  #              are true fits; 2 means twice as many, etc.; 0 = no permutations
  #nCores: number of CPU cores to use for parallelization
  #        if NULL, automatic selection; if 1, sequential execution
  #        if an empty list, an externally registered cluster will be used
  #        note: if method is within and data is a list of subjects, these will be
  #              parallelized. Otherwise the slices will be parallelized.
  #              if slices should be parallelized for the former case, loop this function
  #              and let it parallelize the slices for each participant separately.
  #RETURNS ---
  #the GAT matrix, averaged over all CVs and each individual GAT matrix separately
  if ( repetitions < 1 ) stop( "Repetitions must be >= 1." )
  .eval_package("foreach", load=T); .eval_package("caret")
  args.in = as.list( match.call() )[-c(1,2)]
  #make sure variable names are evaluated before they are passed on
  args.in = lapply(args.in, function(x) if (class(x)=="name") eval(x) else x)
  method = tolower(method)
  decompose = toupper(decompose)
  if ( !method %in% c("within","across","between") ) {
    stop( "Undefined decoding method. Options are 'within', 'across' and 'between'." )
  }
  if ( !decompose %in% c("NONE","CSP","SPECCSP","SPOC") ) {
    stop( "Undefined decomposition method. Options are 'CSP', 'SpecCSP' and 'SPoC'." )
  }
  data = data.set_type(data)
  if ( method == "within" && "subjects" %in% attr(data, "type") ) { #parallelize subjects
    pcheck = .parallel_check(required=length(data), nCores=nCores)
    result = foreach(d=data, .combine=list, .multicombine=T) %dopar%
      do.call(decode.GAT, utils::modifyList( args.in, list(data=d, nCores=1) ))
    .parallel_check(output=pcheck)
    names(result) = paste0("subject", seq_along(data))
    if ( permutations > 0 ) {
      GAT.avg = list( GAT.true = Reduce("+", lapply(result, function(sub) sub[["GAT"]]$GAT.true))/length(result),
                      GAT.random = Reduce("+", lapply(result, function(sub) sub[["GAT"]]$GAT.random))/length(result) )
    } else {
      GAT.avg = Reduce("+", lapply(result, function(sub) sub[["GAT"]]))/length(result)
    }
    return( list( average = GAT.avg, individual = result ) )
  } else if ( method != "within" && !"subjects" %in% attr(data, "type") ) {
    stop( "Please supply the subject data as a list of data frames for this method." )
  } else if ( method == "across" ) { #get list of outcome info for fold split
    outcomes = lapply( data.permute_classes(data, shuffle=F), attr, "outcome" )
  } else if ( method == "within" ) { #get outcome info for fold split
    outcome = attr( data.permute_classes(data, shuffle=F), "outcome" )
  }
  ## retrieve input arguments and prepare
  if ( is.null(permutations) ) permutations = 0
  args.slice = .eval_ellipsis("data.split_slices", ...); args.slice$idx.out = T #save RAM
  if ( "slices" %in% attr(data, "type") && is.null(args.in$window) ) {
    #data is sliced and no new slicing parameters were supplied:
    data = .data.unslice(data)
    args.slice$window = attr(data, "window")
    args.slice$overlap = attr(data, "overlap")
  }
  k = .eval_ellipsis("data.split_folds", ...)$k #get k for the fold split
  tol = .eval_ellipsis("data.permute_classes", ...)$tol #get tol for the permutations
  decomp = decompose != "NONE"
  decomp.method = paste0("decompose.", grep(decompose, c("CSP","SpecCSP","SPoC"), ignore.case=T, value=T)[1])
  decomp.name = substr(decomp.method, 11, 20)
  args.proj = .eval_ellipsis("data.project", ...) #projection arguments for decomp 
  if ( decomp && eval(args.proj$approximate) && is.null(args.in$scale) ) args.in$scale = F
  if ( !is.null(args.in$baseline) ) {
    data = data.remove_samples(data, end=args.in$baseline)
    args.in$baseline = NULL
  }
  if ( nchar( as.character(args.slice$window) ) == 0 ) stop( "Please supply slice settings." )
  sidx = do.call(data.split_slices, modifyList( args.slice, list(data=data) )) #slice indices
  if ( method != "within" ) {
    if (any( diff( sapply(sidx, function(s) length(s$slice.idx) ) ) > 0 )) {
      stop( "Number of slices varies across subjects. Fix the sample numbers to be identical." )
    }
    sidx = sidx[[1]] #first subject generalizes to the rest
    measurements = is.datacol(data[[1]]) #1st subject generalizes to the rest
  } else {
    data = data.check(data, aslist=F)
    measurements = is.datacol(data) #measurement info on all trials
  }
  slicenums = unique(sidx$slice.idx)
  measurements[2] = F #ensure there is no mis-identification
  #initialize parallelization for slices (sequential if already parallelized within subjects)
  pcheck = .parallel_check(required=length(slicenums), nCores=nCores)
  data = data.split_trials(data, strip=F) #change into trial format

  ## internal function to fit train set ##
  slice_fit <- function(train) {
    nsamples = data.samplenum(train) #might vary across slices
    if (decomp) { #decompose slice data and overwrite set
      if ( eval(args.proj$approximate) ) { #1 value per trial
        train.outcome = train[ seq(1, nrow(train), nsamples), 2]
      } else { #copy the outcome column
        train.outcome = train[,2]
      }
      d = do.call(decomp.method, modifyList( args.in, list(data=train) ))
      train = data.frame(1, train.outcome, d$features)
    }
    fit = list( fit = do.call(data.fit, modifyList( args.in, list(train=train) )) )
    if (decomp) fit[[decomp.name]] = d[!grepl("features|outcome", names(d))] #add decomposition output to fit
    return(fit)
  }
  ## internal function to predict slice ## 
  slice_predict <- function(test, fit) {
    if (decomp) { #apply training set filters to test set
      test = do.call(data.project, modifyList( args.proj, list(data=test, weights=fit[[decomp.name]]$filters) ))
    }
    test.outcome = test[,2]
    test = unname( scale(test[, measurements], fit$fit$center, fit$fit$scale) )
    #get predictions
    preds <- predict(fit$fit, test)[[1]]
    preds = as.ordered(preds)
    if ( length(unique(test.outcome)) > 2 ) { #multiclass
      auc = pROC::multiclass.roc(test.outcome, preds)$auc[1]
    } else { #binary
      auc = pROC::roc(test.outcome, preds)$auc[1]
    }
    return(auc)
  }
  
  ######## GAT START ########
  cat("Starting GAT procedure . . .\n")
  GAT.result = list()
  for ( run in seq_len(repetitions) ) {
    reptime = proc.time()[3]
    cat("Run",run,"\n")
    ## To reduce RAM demands, grab data only as needed ##
    #get fold indices
    if ( method == "within" ) { #one subject
      indxFolds = caret::createFolds(y=outcome, k=k) #fold indices are trial numbers  
    } else if ( method == "across" ) { #multiple subjects
      indxFolds = lapply(outcomes, caret::createFolds, k=k) #fold indices are trial numbers per subject
    } else if ( method == "between" ) { #multiple subjects
      indxFolds = createSubjFolds(length(data), k) #fold indices are subject numbers    
    }
    #start iteration over folds
    fold.result = list()
    for ( foldnum in 1:k ) {
      cat("Fold",foldnum,"\n")
      #get training and test sets
      if ( method == "within" ) {
        #trials of one subject are assigned to training/test set
        train.true = data[ -indxFolds[[foldnum]] ] #training trials
        test.true = data[ indxFolds[[foldnum]] ] #test trials
      } else if ( method == "across" ) {
        #stratified fold split for subjects
        train = lapply(seq_along(data), function(sub) {
          data[[sub]][ -indxFolds[[sub]][[foldnum]] ]
        }) #training trials from all subjects
        test = lapply(seq_along(data), function(sub) {
          data[[sub]][ indxFolds[[sub]][[foldnum]] ]
        }) #test trials from all subjects
      } else if ( method == "between" ) {
        #training and test sets comprise all trials but from different subjects
        train = data[ -indxFolds[[foldnum]] ] #training subjects
        test = data[ indxFolds[[foldnum]] ] #test subjects
      }
      if ( method != "within" ) { #concatenate trials
        train.true = do.call(c, train)
        test.true = do.call(c, test)
      }
      if ( k <= 1 ) train.true = test.true #no CV
      
      ## part 1: obtain slice fits
      cat("Fitting . . .")
      fits.true = foreach(i=slicenums, .combine=list, .multicombine=T) %dopar%
        {
          sliceidx = sidx$trial.idx[ sidx$slice.idx == i ]
          slice.train = data.check(lapply( train.true, function(d) na.omit( d[ sliceidx, ] ) ), aslist=F)
          slice_fit(train=slice.train)
        }
      #do fold-wise permutations on train and test set
      fits.random = list()
      for ( perm in seq_len(permutations) ) {
        if ( method == "within" ) {
          train.random = data.permute_classes(train.true, shuffle=T, tol)
        } else { #don't switch labels between subjects
          train.random = do.call(c, lapply(train, data.permute_classes, shuffle=T, tol))
        }
        fits.random[[perm]] = foreach(i=slicenums, .combine=list, .multicombine=T) %dopar%
          {
            sliceidx = sidx$trial.idx[ sidx$slice.idx == i ]
            slice.train = data.check(lapply( train.random, function(d) na.omit( d[ sliceidx, ] ) ), aslist=F)
            slice_fit(train=slice.train)
          }
      }
      cat(" | Time elapsed (mins):", round((proc.time()[3]-reptime)/60,2), "\n")
      
      ## part 2: omnibus prediction of slices with model fits
      cat("Predicting . . .")
      result = list()
      for (s in slicenums) { #iterate the slices
        sliceidx = sidx$trial.idx[ sidx$slice.idx == s ]
        slice.test = data.check(lapply(test.true, function(d) d[ sliceidx, ] ), aslist=F)
        result[[s]] = foreach(i=slicenums, .combine=c, .multicombine=T) %dopar%
          slice_predict(test=slice.test, fit=fits.true[[i]]) #iterate the fits
      }
      GAT.true = do.call(rbind, result) #columns: fit, rows: generalization
      rownames(GAT.true) = paste0("test", slicenums)
      colnames(GAT.true) = paste0("train", slicenums)
      #permutations
      GAT.random = list()
      for ( perm in seq_len(permutations) ) {
        tmp = list()
        if ( method == "within" ) {
          test.random = data.permute_classes(test.true, shuffle=T, tol)
        } else { #don't switch labels between subjects
          test.random = do.call(c, lapply(test, data.permute_classes, shuffle=T, tol))
        }
        for (s in slicenums) { #iterate the slices
          sliceidx = sidx$trial.idx[ sidx$slice.idx == s ]
          slice.test = data.check(lapply(test.random, function(d) d[ sliceidx, ] ), aslist=F)
          tmp[[s]] = foreach(i=slicenums, .combine=c, .multicombine=T) %dopar%
            slice_predict(test=slice.test, fit=fits.random[[perm]][[i]]) #iterate the fits
        }
        GAT.random[[perm]] = do.call(rbind, tmp)
        rownames(GAT.random[[perm]]) = paste0("test", slicenums)
        colnames(GAT.random[[perm]]) = paste0("train", slicenums)
      }
      fold.result[[foldnum]] = list(GAT.true = GAT.true)
      if ( permutations > 0 ) {
        fold.result[[foldnum]]$GAT.random = setNames( GAT.random, paste0("permutation", seq_len(permutations)) )
      }
      cat(" | Time elapsed (mins):", round((proc.time()[3]-reptime)/60,2), "\n")
    } #end of fold
    GAT.result[[run]] = setNames(fold.result, paste0("fold",1:k))
  } #end of rep
  GAT.result = setNames(GAT.result, paste0("run", seq_len(repetitions)) )
  ######## GAT END ########
  .parallel_check(output=pcheck) #prints final time      
      
  #compile average for output
  GAT.avg = list( #aggregate repetitions, folds
    GAT.true = Reduce("+", lapply(GAT.result, function(run) { 
      Reduce("+", lapply(run, "[[", "GAT.true"))/k }) )/repetitions 
    )
  if ( permutations > 0 ) { #aggregate repetitions, folds, num permutations
    GAT.avg$GAT.random = Reduce("+", lapply(GAT.result, function(run) {
      Reduce("+", lapply(run, function(fold) {
        Reduce("+", fold[["GAT.random"]])/permutations }) )/k 
      }) )/repetitions
  }
  output = list( average = GAT.avg, individual = GAT.result )
  return( .unnest(output) ) #return and remove nesting where only 1 element
}

decode.eval <- function(train=NULL, test=NULL, fits=NULL, decompose="none", 
                        method="classification", nCores=NULL, package=NULL, ...) {
  ## a more flexible version of decode which leaves the CV to the user
  ## if train and not fits is supplied, fits will be generated; if additionally test is supplied,
  ## the obtained model fits are evaluated on the test set
  ## if test and fits is supplied without train, the fits will be evaluated on the test data
  ## if fits is a list corresponding to slices, test should either be sliced as well
  ## or the appropriate slicing parameters have to be provided
  ## if decomposition was used to create fits, decompose should be set accordingly and
  ## the arguments approximate/logtransform should be provided in case they differed from their defaults
  #INPUT ---
  #train, test: data frame, list of trials, slices, subjects
  #fits: either directly a model fit, a function output containing the fit as an element,
  #      or a list of slices, each containing a fit (e.g. as obtained via decode)
  #      or a list of subjects corresponding to subject data
  #nCores: number of CPU cores to use for parallelization
  #        if NULL, automatic selection; if 1, sequential execution
  #        if an empty list, an externally registered cluster will be used
  #        note: if data is a list of subjects, these will be parallelized.
  #              otherwise slices will be parallalized.
  #package: if parallelized & in cases where the predict function needs a certain library to work,
  #         supply the character string of the respective library
  #... : further arguments correspond to fit, slice and decomposition settings if applicable
  #RETURNS ---
  #a list with fits (if train was supplied) and summary, predictions (if test was supplied)
  doFit = is.null(fits) && !is.null(train) #a training set and no fits supplied
  doEval = !is.null(test) && xor( !is.null(train), !is.null(fits) ) #a test set to evaluate the fits and either train or fits
  if ( !doFit && !doEval ) stop( "Supply at least one of train or test and don't supply both train and fits." )
  method = tolower(method)
  if ( !method %in% c("classification", "regression") ) {
    stop( "Undefined evaluation method. Options are 'classification' and 'regression'." )
  }
  decompose = toupper(decompose)
  if ( !decompose %in% c("NONE","CSP","SPECCSP","SPOC") ) {
    stop( "Undefined decomposition method. Options are 'CSP', 'SpecCSP' and 'SPoC'." )
  }
  if (doFit) {
    type = attr(data.set_type(train), "type")
  } else { #evaluate only
    type = attr(data.set_type(test), "type")
  }
  if ( identical("subjects", type) && doEval && !doFit && (length(fits) != length(test)) ) {
    stop( "Number of fits does not seem to match number of subjects.",
          "Supply a list of fits that corresponds to a list of subject data." )
  }
  args.in = as.list( match.call() )[-1]
  args.in = args.in[!grepl("train|test|fits", names(args.in))] #remove main input
  #make sure variable names are evaluated before they are passed on
  args.in = lapply(args.in, function(x) if (class(x)=="name") eval(x) else x)
  if ( method == "regression" && is.null(args.in$model) ) args.in$model = "Reg"
  .eval_package("foreach", load=T)
  if ( identical("subjects", type) ) { #parallelize subjects
    req = ifelse(doFit, length(train), length(test)) #required CPU cores
    pcheck = .parallel_check(required=req, nCores=nCores)
    result = foreach(i=seq_len(req), .combine=list, .multicombine=T, .packages=package) %dopar%
      do.call(decode.eval, utils::modifyList( args.in, list( train=train[[i]], test=test[[i]], fits=fits[[i]], nCores=1 ) ))
    .parallel_check(output=pcheck)
    result = .output_merge(result)
    if ( doEval && !"result" %in% names(result) ) { #no slicing
      tmp = sapply(result[[1]], "[[", 1) #collect each subject's performance score
      result[[1]] = setNames( data.frame(seq_len(req), tmp, row.names=NULL), c("subject", names(result)[1]) )
      names(result)[1] = "result"
    }
    return(result)
  }
  #evaluate input arguments
  args.slice = .eval_ellipsis("data.split_slices", ...); args.slice$idx.out = T #save RAM
  if ( identical("slices", type) && is.null(args.in$window) ) {
    #data is sliced and no new slicing parameters were supplied:
    if (doFit) train = .data.unslice(train)
    if (doEval) test = .data.unslice(test)
    args.slice$window = ifelse(doFit, attr(train, "window"), attr(test, "window"))
    args.slice$overlap = ifelse(doFit, attr(train, "overlap"), attr(test, "overlap"))
  }
  slice = nchar( as.character(args.slice$window) ) > 0 #should data get sliced
  decomp = decompose != "NONE"
  decomp.method = paste0("decompose.", grep(decompose, c("CSP","SpecCSP","SPoC"), ignore.case=T, value=T)[1])
  decomp.name = substr(decomp.method, 11, 20)
  args.proj = .eval_ellipsis("data.project", ...) #projection arguments for decomp 
  if ( !is.null(args.in$baseline) ) {
    if (doFit) train = data.remove_samples(train, end=args.in$baseline)
    if (doEval) test = data.remove_samples(test, end=args.in$baseline)
    args.in$baseline = NULL
  }
  if (slice) {
    if (doFit) { #obtain slice indices on train
      sidx = do.call(data.split_slices, modifyList( args.slice, list(data=train) )) #slice indices
    } else { #else on test
      sidx = do.call(data.split_slices, modifyList( args.slice, list(data=test) )) #slice indices
    }
  } else { #imitate a single slice
    if (doFit) {
      train = data.check(train, aslist=F) #make sure data is a df
      nsamples = data.samplenum(train)
    } else { #doEval
      test = data.check(test, aslist=F) #make sure data is a df
      nsamples = data.samplenum(test)
    }
    sidx = list(trial.idx = 1:nsamples, slice.idx = rep(1, nsamples)) #every sample per trial
  }
  slicenums = unique(sidx$slice.idx)
  output = list() #prepare function output
  #initialize parallelization for slices (sequential if already parallelized within subjects)
  pcheck = .parallel_check(required=length(slicenums), nCores=nCores)
  
  ### fitting procedure
  if (doFit) {
    fittime = proc.time()[3]
    if ( decomp && eval(args.proj$approximate) && is.null(args.in$scale) ) args.in$scale = F
    train = data.split_trials(train, strip=F) #change into trial format
    ## internal function to fit train set ##
    slice_fit <- function(train) {
      nsamples = data.samplenum(train) #might vary across slices
      if (decomp) { #decompose slice data and overwrite set
        if ( eval(args.proj$approximate) ) { #1 value per trial
          train.outcome = train[ seq(1, nrow(train), nsamples), 2]
        } else { #copy the outcome column
          train.outcome = train[,2]
        }
        d = do.call(decomp.method, modifyList( args.in, list(data=train) ))
        train = data.frame(1, train.outcome, d$features)
      }
      fit = list( fit = do.call(data.fit, modifyList( args.in, list(train=train) )) )
      if (decomp) fit[[decomp.name]] = d[!grepl("features|outcome", names(d))] #add decomposition output to fit
      return(fit)
    }
    ## start fitting ##
    cat("Fitting . . .")
    fits = foreach(i=slicenums, .combine=list, .multicombine=T) %dopar%
      {
        sliceidx = sidx$trial.idx[ sidx$slice.idx == i ]
        slice.train = data.check(lapply( train, function(d) na.omit( d[ sliceidx, ] ) ), aslist=F)
        slice_fit(train=slice.train)
      }
    
    if (!slice) fits = list(fits)
    fits = setNames(fits, paste0("slice", slicenums))
    #add to function output
    if (slice) output$fits = fits else output$fit = fits[[1]]
    if (!slice && !decomp) output$fit = output$fit[[1]] #remove nesting
    cat(" | Time elapsed (mins):", round((proc.time()[3]-fittime)/60,2), "\n")
  } #end of doFit
  
  ### evaluation procedure
  if (doEval) {
    evaltime = proc.time()[3]
    .eval_package( c("caret", "pROC") )
    if ( "fit" %in% attr(fits, "type") || ( is.list(fits) && class(fits) != "list" ) ) {
      fits = list( list(fit=fits) ) #output of data.fit, mimick fit of 1 slice
    } else if ( "fit" %in% names(fits) ) {
      fits = list(fits) #wrap as fit of 1 slice
    }
    if ( length(slicenums) != length(fits) ) stop( "Number of fits does not seem to match number of slices." )
    test = data.check(test, aslist=F)
    if ( !exists("nsamples", inherits=F, mode="numeric") ) nsamples = data.samplenum(test)
    measurements = is.datacol(test); measurements[2] = F #make sure outcome is not mis-classified
    outcome = test[ seq(1, nrow(test), nsamples), 2 ] #1 value per trial
    if ( method == "classification" ) {
      outcome = as.ordered(outcome) 
    } else { #regression
      outcome = suppressWarnings( as.numeric(outcome) ) #suppress warnings in case of NAs
      if ( any( is.na(outcome) ) ) { #check if NAs were introduced due to non-numeric characters
        outcome = test[ seq(1, nrow(test), nsamples), 2 ]
        outcome = as.numeric( as.factor(outcome) )
      }
    }
    test = data.split_trials(test, strip=F) #change into trial format
    ## internal prediction function ## 
    slice_predict <- function(test, fit) {
      if (decomp) { #apply training set filters to test set
        test = do.call(data.project, modifyList( args.proj, 
                            list(data=test, weights=fit[[decomp.name]]$filters) ))
      }
      test = data.check(test, aslist=F)
      test = unname( scale(test, fit$fit$center, fit$fit$scale) )
      out = rep( outcome, each=nrow(test)/length(outcome) )
      #get predictions
      preds <- predict(fit$fit, test)[[1]]
      if ( method == "classification" ) {
        preds = as.ordered(preds)
        if ( length(unique(out)) > 2 ) { #multiclass
          auc = pROC::multiclass.roc(out, preds)$auc[1]
        } else { #binary
          auc = pROC::roc(out, preds)$auc[1]
        }
        res = list( AUC = auc, predictions = caret::confusionMatrix(preds, out,
                                                dnn=c("Predicted Label", "True Label"))$table )
      } else { #regression
        res = list( RMSE = sqrt( mean((out - preds)^2) ), #Root mean squared error
                    deviation = summary(out - preds ) ) #distribution of errors
      }
      return(res)
    }
    ## start evaluating ##
    cat("Predicting . . .")
    result = foreach(i=slicenums, .combine=list, .multicombine=T, .packages=package) %dopar%
      { 
        sliceidx = sidx$trial.idx[ sidx$slice.idx == i ]
        slice.test = lapply(test, function(d) d[ sliceidx, measurements ] )
        attr(slice.test, "type") = "trials"
        slice_predict(test=slice.test, fit=fits[[i]])
      }
    
    if (slice) { #prepare output
      output$result = data.frame(slice=slicenums)
      output$result[[ names(result[[1]])[1] ]] = sapply(result, "[[", 1) #add AUC/RMSE to result df
      #add predictions/deviation to output
      output[[ names(result[[1]])[2] ]] = setNames(lapply(result, "[[", 2), paste0("slice", slicenums))
    } else {
      output = c(output, result)
    }
    cat(" | Time elapsed (mins):", round((proc.time()[3]-evaltime)/60,2), "\n")
  } #end of doEval
  .parallel_check(output=pcheck) #prints final time
  #re-order to match typical output ordering
  if (doFit && doEval) output = output[c(2,3,1)] #AUC/result, predictions, fits
  return(output)
}


decode <- function(data, method="within", decompose="none", repetitions=1, 
                   permutations=1, verbose=F, nCores=NULL, ...) {
  ## the heart of the decoding package: decode with repeated k-fold CV
  ## in different variations (within/across/between subjects)
  #INPUT ---
  #data: continuous df, list of trials, slices, subjects
  #method: one of c("within","across","between")
  #        within: single subject decoding; either a data set of a single subject
  #                or a list of subject data sets that will be parallelized
  #        across: data of multiple subjects is concatenated where folds are
  #                stratified for subjects, i.e. train and test set contain
  #                trials of every subject equally
  #        between: data of multiple subjects is separated into train and test set
  #                 i.e. train and test set contain trials of different subjects
  #decompose: decomposition method to use, one of c("CSP","SpecCSP","SPoC")
  #repetitions: number of repetitions of the k-fold CV, a value >= 1
  #permutations: number of permutations on each fold, corresponds to the number
  #              of true fits, i.e. a value of 1 means as many permutations as there 
  #              are true fits; 2 means twice as many, etc.; 0 = no permutations
  #verbose: if True, fits are appended to output. Will take up a lot of memory
  #         if many slices/repetitions are performed.
  #nCores: number of CPU cores to use for parallelization
  #        if NULL, automatic selection; if 1, sequential execution
  #        if an empty list, an externally registered cluster will be used
  #        note: if method is within and data is a list of subjects, these will be
  #              parallelized. Otherwise the slices will be parallelized.
  #              if slices should be parallelized for the former case, loop this function
  #              and let it parallelize the slices for each participant separately.
  #RETURNS ---
  #summary: significance evaluation of the result
  #result: data.frame with the performance scores 
  #fits: list with model fits
  if ( repetitions < 1 ) stop( "Repetitions must be >= 1." )
  .eval_package("foreach", load=T); .eval_package("caret")
  args.in = as.list( match.call() )[-c(1,2)]
  #make sure variable names are evaluated before they are passed on
  args.in = lapply(args.in, function(x) if (class(x)=="name") eval(x) else x)
  args.test = .eval_ellipsis("decode.test", ...)
  method = tolower(method)
  decompose = toupper(decompose)
  if ( !method %in% c("within","across","between") ) {
    stop( "Undefined decoding method. Options are 'within', 'across' and 'between'." )
  }
  if ( !decompose %in% c("NONE","CSP","SPECCSP","SPOC") ) {
    stop( "Undefined decomposition method. Options are 'CSP', 'SpecCSP' and 'SPoC'." )
  }
  data = data.set_type(data)
  if ( method == "within" && "subjects" %in% attr(data, "type") ) { #parallelize subjects
    pcheck = .parallel_check(required=length(data), nCores=nCores)
    result = foreach(d=data, .combine=list, .multicombine=T) %dopar%
      do.call(decode, utils::modifyList( args.in, list(data=d, nCores=1) ))
    result = .output_merge(result)
    summary = do.call(decode.test, utils::modifyList( args.test, 
          list(result=result$result, between=ifelse("slice" %in% names(result$summary), "slice", ""), nCores=list() )))
    result = c(list(summary.overall = summary), result)
    .parallel_check(output=pcheck)
    return(result)
  } else if ( method != "within" && !"subjects" %in% attr(data, "type") ) {
    stop( "Please supply the subject data as a list of data frames for this method." )
  } else if ( method == "across" ) { #get list of outcome info for fold split
    outcomes = lapply( data.permute_classes(data, shuffle=F), attr, "outcome" )
  } else if ( method == "within" ) { #get outcome info for fold split
    outcome = attr( data.permute_classes(data, shuffle=F), "outcome" )
  }
  ## retrieve input arguments and prepare
  if ( is.null(permutations) ) permutations = 0
  args.slice = .eval_ellipsis("data.split_slices", ...); args.slice$idx.out = T #save RAM
  if ( "slices" %in% attr(data, "type") && is.null(args.in$window) ) {
    #data is sliced and no new slicing parameters were supplied:
    data = .data.unslice(data)
    args.slice$window = attr(data, "window")
    args.slice$overlap = attr(data, "overlap")
  }
  k = .eval_ellipsis("data.split_folds", ...)$k #get k for the fold split
  tol = .eval_ellipsis("data.permute_classes", ...)$tol #get tol for the permutations
  decomp = decompose != "NONE"
  decomp.method = paste0("decompose.", grep(decompose, c("CSP","SpecCSP","SPoC"), ignore.case=T, value=T)[1])
  decomp.name = substr(decomp.method, 11, 20)
  args.proj = .eval_ellipsis("data.project", ...) #projection arguments for decomp 
  if ( decomp && eval(args.proj$approximate) && is.null(args.in$scale) ) args.in$scale = F
  if ( !is.null(args.in$baseline) ) {
    data = data.remove_samples(data, end=args.in$baseline)
    args.in$baseline = NULL
  }
  slice = nchar( as.character(args.slice$window) ) > 0 #should data get sliced
  if (slice) {
    sidx = do.call(data.split_slices, modifyList( args.slice, list(data=data) )) #slice indices
    if ( method != "within" ) {
      if (any( diff( sapply(sidx, function(s) length(s$slice.idx) ) ) > 0 )) {
        stop( "Number of slices varies across subjects. Fix the sample numbers to be identical." )
      }
      sidx = sidx[[1]] #first subject generalizes to the rest
    }
  } else { #imitate a single slice
    if ( method == "within" ) {
      data = data.check(data, aslist=F) #make sure data is a df
      nsamples = data.samplenum(data)
    } else {
      nsamples = max( sapply(data, data.samplenum) ) #sample numbers may differ
    }
    sidx = list(trial.idx = 1:nsamples, slice.idx = rep(1, nsamples)) #every sample per trial
  }
  slicenums = unique(sidx$slice.idx)
  #initialize parallelization for slices (sequential if already parallelized within subjects)
  pcheck = .parallel_check(required=length(slicenums), nCores=nCores)
  data = data.split_trials(data, strip=F) #change into trial format
  
  ## internal functions ##
  slice_fit <- function(train, test) {
    ## internal function to fit train and test slice sets (lists of trials)
    if (decomp) {
      #decompose slice data and overwrite sets
      d = slice_decompose(train, test)
      train = d$train; test = d$test
    }
    fit = do.call(data.fit, modifyList( args.in, list(train=train, test=test) ))
    if (decomp) fit[[decomp.name]] = d$decomp #add decomposition output to fit
    return(fit)
  }
  slice_decompose <- function(train, test) {
    ## helper to prepare train and test set with decomposition procedure
    nsamples = data.samplenum(train) #might vary across slices
    if ( eval(args.proj$approximate) ) { #1 value per trial
      train.outcome = train[ seq(1, nrow(train), nsamples), 2]
      test.outcome = test[ seq(1, nrow(test), nsamples), 2]    
    } else { #copy the outcome column
      train.outcome = train[,2]
      test.outcome = test[,2]
    }
    d = do.call(decomp.method, modifyList( args.in, list(data=train) )) #decompose training set
    train = data.frame(1, train.outcome, d$features)
    test = data.frame(1, test.outcome, #apply training set filters to test set
                      do.call(data.project, modifyList( args.proj, list(data=test, weights=d$filters) )))
    return( list(train=train, test=test, decomp=d[!grepl("features|outcome", names(d))]) )
  }
  
  ######## FITTING START ########
  checktime = proc.time()[3]
  cat("Starting to decode . . .\n")
  result = list()
  for ( run in seq_len(repetitions) ) {
    reptime = proc.time()[3]
    cat("Run",run,"\n")
    ## To reduce RAM demands, grab data only as needed ##
    #get fold indices
    if ( method == "within" ) { #one subject
      indxFolds = caret::createFolds(y=outcome, k=k) #fold indices are trial numbers  
    } else if ( method == "across" ) { #multiple subjects
      indxFolds = lapply(outcomes, caret::createFolds, k=k) #fold indices are trial numbers per subject
    } else if ( method == "between" ) { #multiple subjects
      indxFolds = createSubjFolds(length(data), k) #fold indices are subject numbers    
    }
    #start iteration over folds
    fold.fits = list()
    for ( foldnum in 1:k ) {
      cat("fold",foldnum)
      #get training and test sets
      if ( method == "within" ) {
        #trials of one subject are assigned to training/test set
        train.true = data[ -indxFolds[[foldnum]] ] #training trials
        test.true = data[ indxFolds[[foldnum]] ] #test trials
      } else if ( method == "across" ) {
        #stratified fold split for subjects
        train = lapply(seq_along(data), function(sub) {
          data[[sub]][ -indxFolds[[sub]][[foldnum]] ]
        }) #training trials from all subjects
        test = lapply(seq_along(data), function(sub) {
          data[[sub]][ indxFolds[[sub]][[foldnum]] ]
        }) #test trials from all subjects
      } else if ( method == "between" ) {
        #training and test sets comprise all trials but from different subjects
        train = data[ -indxFolds[[foldnum]] ] #training subjects
        test = data[ indxFolds[[foldnum]] ] #test subjects
      }
      if ( method != "within" ) { #concatenate trials
        train.true = do.call(c, train)
        test.true = do.call(c, test)
      }
      if ( k <= 1 ) train.true = test.true #no CV
      fits.true = foreach(i=slicenums, .combine=list, .multicombine=T) %dopar%
        {
          sliceidx = sidx$trial.idx[ sidx$slice.idx == i ]
          slice.train = data.check(lapply( train.true, function(d) na.omit( d[ sliceidx, ] ) ), aslist=F)
          slice.test = data.check(lapply( test.true, function(d) na.omit( d[ sliceidx, ] ) ), aslist=F)
          slice_fit(slice.train, slice.test)
        }
      if ( !slice ) fits.true = list(fits.true)
      #do fold-wise permutations on train and test set
      fits.random = list()
      for ( perm in seq_len(permutations) ) {
        if ( method == "within" ) {
          train.random = data.permute_classes(train.true, shuffle=T, tol)
          test.random = data.permute_classes(test.true, shuffle=T, tol)
        } else { #don't switch labels between subjects
          train.random = do.call(c, lapply(train, data.permute_classes, shuffle=T, tol))
          test.random = do.call(c, lapply(test, data.permute_classes, shuffle=T, tol))
        }
        fits.random[[perm]] = foreach(i=slicenums, .combine=list, .multicombine=T) %dopar%
          {
            sliceidx = sidx$trial.idx[ sidx$slice.idx == i ]
            slice.train = data.check(lapply( train.random, function(d) na.omit( d[ sliceidx, ] ) ), aslist=F)
            slice.test = data.check(lapply( test.random, function(d) na.omit( d[ sliceidx, ] ) ), aslist=F)
            slice_fit(slice.train, slice.test)
          }
        if ( !slice ) fits.random[[perm]] = list( fits.random[[perm]] )
      }
      cat(" | Time elapsed (mins):", round((proc.time()[3]-reptime)/60,2), "\n")
      fits.true = setNames( fits.true, paste0("slice", slicenums) )
      fold.fits[[foldnum]] = list(fits.true, fits.random)
    } #end of folds
    fold.fits = setNames( fold.fits, paste0("fold", 1:k) )
    #create summary df
    acc.true = do.call(c, lapply(fold.fits, function(fit) sapply(fit[[1]], "[[", 1)) )
    summary = data.frame( run=run, fold=rep(1:k, each=length(slicenums)), 
                          slice=slicenums, label="true", acc=acc.true, row.names=NULL )
    if (permutations > 0) {
      acc.random = do.call(c, lapply(fold.fits, function(fit) sapply(fit[[2]], function(rep) sapply(rep, "[[", 1)) ) )
      summary = rbind(summary, data.frame( run=run, fold=rep(1:k, each=length(slicenums)*permutations), 
                                           slice=slicenums, label="random", acc=acc.random, row.names=NULL ))
    }
    #replace placeholder name with actual performance measurement name
    fnames = names(fold.fits$fold1[[1]]$slice1)[1:2] #1st fold, true fits, 1st slice
    names(summary)[5] = fnames[1]
    output = list(summary = summary, fits = lapply(fold.fits, "[[", 1)) #true label fits
    if (!verbose) { #replace with aggregated predictions
      output$fits = lapply(slicenums, function(s) { Reduce( "+", lapply(1:k, function(f) {
        output$fits[[f]][[s]][[ fnames[2] ]] }) )/k }) #average over folds per slice
    }
    result[[run]] = output
  } #end of runs
  result = setNames(result, paste0("run", seq_len(repetitions)) )
  ######## FITTING END ########
  ## compile output
  resultdf = plyr::rbind.fill( lapply(result, "[[", "summary") )
  resultdf = resultdf[ order(resultdf$label, resultdf$slice, resultdf$run, resultdf$fold), ]
  row.names(resultdf) = NULL
  if (!slice) resultdf = resultdf[,-3] #remove slice column
  if (repetitions < 2) resultdf = resultdf[,-1] #remove run column
  output = list(result = resultdf)
  if ( permutations > 0 ) { #significance testing
    cat("Evaluating result . . .\n")
    summary = do.call(decode.test, modifyList( args.test, 
                         list(result=resultdf, id="", between=ifelse(slice, "slice", ""), nCores=list()) ))
    output = c( list(summary=summary), output ) #add summary to first position
  }
  fits = lapply(result, "[[", "fits")
  if (verbose) { #add entire fits to output
    output$fits = fits
  } else { #add only aggregation of predictions to output
    output[[ fnames[2] ]] = setNames( lapply(slicenums, function(s) {
      Reduce("+", lapply(fits, "[[", s))/repetitions }), paste0("slice",slicenums) )
  }
  .parallel_check(output=pcheck) #prints final time
  return( .unnest(output) ) #return and remove nesting where necessary
}

.output_merge <- function(output, newStr = "subject") {
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
      # idxname = ifelse( grepl("fits\\.", elem), "", "fits." )
      out[[elem]] = temp #append to list
    }
  }
  return(out)
}

decode.test <- function(result, p.vals=T, bootstrap=F, dv="AUC", within.var="label", within.null="random",
                        between="slice", id="subject", adjust="none", alpha=0.05, iterations=10000, nCores=NULL) {
  ## non-parametric testing at all slices for significant difference in decodability 
  ## between true and random labels. all operations are stratified for id.
  #INPUT ---
  #result: df with columns specified by dv, pair and group
  #        or a list of such dfs corresponding to different subjects
  #p.vals: if True, p-values are computed by switching the labels of performance scores randomly
  #        and computing the mean within each permuted label category on each iteration. 
  #        p-value = the number of times the permuted mean was higher than the true mean
  #bootstrap: if True, resamples the random labels with replacement and extracts the
  #           confidence level on each iteration. Yields as many values as iterations
  #           which are averaged to get the final CI to which the true performance
  #           will be compared. Only sensible if too few permutations were performed.
  #           if False, the CIs are computed directly on the dv values of within.null
  #dv: dependent variable, column with the performance score; required
  #within.var: pair to compare, e.g. column with the labels (true vs. random); required
  #within.null: within.var's value to compute the CIs with (representing the null hypothesis)
  #between: variable that groups the dv into separate segments; can be left empty
  #id: variable that stratifies all computations into independent parts; can be left empty
  #adjust: correction method for the p-values, default is "none"
  #alpha: significance threshold to test p-values against, CIs are computed for this level
  #iterations: number of times to repeat the permutations and/or bootstrap 
  #RETURNS ---
  #a summary df with means, CIs, p-values, significance indication
  if ( any(alpha > 1) || any(alpha < 0) ) stop( "alpha must be in the range 0|1." )
  if ( class(result) == "list" ) result = plyr::ldply(result, .id=id)
  #check the columns
  if ( !within.var %in% names(result) ) stop( "Cannot find a column named '", within.var, "'." )
  if ( !dv %in% names(result) ) stop( "Cannot find a column named '", dv, "'." )
  between = ifelse( is.null(between), "", between ) #replace NULL argument with empty string
  if ( is.null(id) || !nzchar(id) ) {
    id = ".tmp"
    result[[id]] = 1 #temporary column to split if no id supplied
  }
  result = result[ order(result[[id]]), ] #sort by id for stratification later
  if ( any( c(bootstrap, p.vals) ) ) { #prepare parallelization if necessary
    .eval_package("foreach", load=T)
    if ( nzchar(between) ) { #parallelize time points (between)
      res.groups = split(result, result[,between])
    } else {
      res.groups = list(result) #list of 1 element for sequential backend
    }
    pcheck = .parallel_check(required=length(res.groups), nCores=nCores)
  }
  alpha = sort(alpha, decreasing=T) #in case of multiple values
  if (p.vals) { #get permutation-based p-values
    permutefun <- function(res) { #permutations under the null
      nulldiff = replicate(iterations, {
        #shuffle labels while retaining group structure (stratification for id)
        perm = as.vector( sapply(unique( res[[id]] ), function(sub) {
          sample( res[[within.var]][ res[[id]] == sub ] , replace=F ) }) )
        null1 = res[[dv]][ perm != within.null ] #'true' labels
        null2 = res[[dv]][ perm == within.null ] #'random' labels
        mean(null1) - mean(null2) #(null) difference after permutation
      })
      truediff = mean( res[[dv]][ res[[within.var]] != within.null ] ) - 
        mean( res[[dv]][ res[[within.var]] == within.null ] )
      #one-sided p-value: number of times nulldiff >= truediff / N
      return( mean( nulldiff >= truediff ) )
    }
    pvals = foreach(res=res.groups, .combine=c, .multicombine=T) %dopar%
      permutefun(res)
  }
  #get CIs
  if (bootstrap) { #obtain CI thresholds via bootstrap
    bootfun <- function(null) { #computes the mean of n bootstrapped quantiles
      CI = matrix( replicate(iterations, {
        #sampling with replacement stratified for id
        x = as.vector( sapply(unique( null[[id]] ), function(sub) {
          sample( null[[dv]][ null[[id]] == sub ] , replace=T ) }) )
        quantile(x, probs=1-alpha)
      }), nrow=length(alpha) ) #matrix with row=threshold, col=iteration
      return( rowMeans(CI) ) #mean threshold
    }
    CIs = foreach(res=res.groups, .combine=rbind, .multicombine=T) %dopar%
      bootfun( res[ res[[within.var]] == within.null, ] ) #only null performance
  } else { #obtain CI thresholds directly on null performance
    null = result[ result[[within.var]] == within.null, ] #only null performance
    if ( nzchar(between) ) { #calculate per time point
      CIs = as.matrix( aggregate( as.formula(paste(dv,"~",between)), data=null, FUN=quantile, probs=1-alpha )[[dv]] )
    } else { #no grouping at all
      CIs = quantile( null[[dv]], probs=1-alpha )
    }
  }
  if ( !nzchar(between) ) CIs = matrix(CIs, ncol=length(alpha)) #make sure indexing works on columns for output
  if ( any( c(bootstrap, p.vals) ) ) .parallel_check(output=pcheck)
  #collect result summary df
  summary = aggregate( as.formula( paste(dv,"~",sub("^\\+", "", paste0(between,"+",within.var))) ), data=result, mean )
  #get mean difference
  meandiff = summary[[dv]][ summary[[within.var]] != within.null ] - summary[[dv]][ summary[[within.var]] == within.null ]
  #reduce summary to the final parameters
  summary = summary[ summary[[within.var]] != within.null, ]
  summary[["effect"]] = meandiff
  for ( i in seq_along(alpha) ) {
    summary[[ paste0("CI",100-alpha[i]*100) ]] = CIs[,i] #add CIs to summary
  }
  if (p.vals) {
    summary$p = pvals
    if ( tolower(adjust) != "none" ) summary$p.adj = p.adjust( pvals, method=adjust )
    sign.sum = sapply(summary[[ncol(summary)]], function(p) sum(p<alpha)) #count pvals < alpha (might be multiple alphas)
    summary$sign = sapply(sign.sum, function(s) paste0( rep("*",s), collapse="" )) #significance indicator
  }
  summary = summary[, !names(summary) %in% within.var] #remove obsolete within column
  return( summary )
}

data.permute_classes <- function(data, shuffle=T, tol=1) {
  ## extracts class labels from data and permutes if requested
  #INPUT ---
  #data: df or list of trials, slices, subjects
  #      with 1st col: sample number, 2nd col: outcome
  #      note: a single outcome per trial is assumed
  #shuffle: if True, labels are shuffled and overwritten
  #tol: if shuffle is True, the permutation success is verified
  #     by checking for each class if half the labels +/- tol * that
  #     have switched their position; e.g. with 100 trials for a class,
  #     50 of those have to be switched for perfect permutation
  #     if tol = 0, exactly 50 have to be switched
  #     if tol = 1, it can be any number (0:100)
  #     with tol at .2, the tolerance range will be from 30:70
  #     i.e. (50-100*.2):(50+100*.2)
  #RETURNS ---
  #data as input, with shuffled labels if shuffle = T
  #outcome (1 value per trial) as attribute
  data = data.set_type(data)
  if ( "subjects" %in% attr(data, "type") ) {
    data = lapply(data, data.permute_classes, shuffle, tol)
    attr(data, "type") = "subjects"
    return(data)
  } 
  if ( class(data) != "list" ) { #df/matrix
    nsamples = data.samplenum(data)
    outcome = data[ seq(1, nrow(data), by=nsamples), 2 ]
  } else { #list
    nsamples = data.samplenum(data[[1]])
    if ( "trials" %in% attr(data, "type") ) {
      outcome = unname( sapply(data, "[[", 1, 2) ) #first sample of every trial df
    } else if ( "slices" %in% attr(data, "type") ) {
      outcome = data[[1]][ seq(1, nrow(data[[1]]), by=nsamples), 2 ] #same order in every slice
    }
  }
  if (shuffle) { #shuffle outcome 
    #half the labels (+/- tol) of each class have to switch positions
    #set tol to 1 if no such verification is desired
    classes = table(outcome)
    success = F
    while ( !success ) {
      tmp = sample(outcome) #shuffle
      success = all( sapply( names(classes), function(class) { #iterate class labels
        #range in which permutaiton was successful:
        permute.range = floor(classes[class]/2 - classes[class]/2*tol):ceiling(classes[class]/2 + classes[class]/2*tol)
        sum( which(outcome == class) %in% which(tmp == class) ) %in% permute.range
      }) )
      #make sure even with tol=1 the shuffle is not identical to the true labels:
      if ( identical(tmp, outcome) ) success=F
    }
    outcome = tmp
    #change labels within the data
    if (class(data) != "list") { #df/matrix
      data[,2] = rep(outcome, each=nsamples)
    } else { #list
      if ( "trials" %in% attr(data, "type") ) {
        data = lapply( seq_along(data), function(i) {
          data[[i]][,2] = outcome[i] #replace
          data[[i]] #return
        })
        attr(data, "type") = "trials"
      } else if ( "slices" %in% attr(data, "type") ) {
        nsamples = sapply(data, data.samplenum) #sample nums per slice
        data = lapply(1:length(data), function(i) {
          data[[i]][,2] = rep(outcome, each=nsamples[i]) #repeat outcome for window size
          data[[i]]
        })
        attr(data, "type") = "slices"
      }
    }
  } #end of shuffle
  attr(data, "outcome") = outcome
  return(data)
}

data.split_partitions <- function(data, ratio=0.7, idx.out=F) {
  ## splits data into training and test data set according to the ratio
  ## the split is balanced for outcome and stratified for trials
  #INPUT ---
  #data: df, list of trials/slices
  #ratio: percentage of data that goes into training
  #idx.out: if False, data is returned as train and test set, 
  #         otherwise the the trial indices of the training set are returned
  #RETURNS ---
  #either a numeric vector with the trial indices (if idx.out=T) or a list
  #with train/test set. If data was sliced, the output list will be a slice
  #list with train and test elements on each slice
  .eval_package("caret")
  data = data.permute_classes(data, shuffle=F)
  if ( "subjects" %in% attr(data, "type") ) {
    data = lapply(data, data.split_partitions, ratio, split)
    attr(data, "type") = "subjects"
    return(data)
  }
  outcome = attr(data, "outcome")
  trainidx = caret::createDataPartition(y=outcome, p=ratio)[[1]]
  if (idx.out) return(trainidx)
  #split data into train and test set
  if ( "slices" %in% attr(data, "type") ) { #data is a slice list
    #iterate the slices and split
    splitlist = setNames( lapply( seq_along(data), function(slice) {
      #get slice data into trial format
      slice.data = data.split_trials( data[[slice]] )
      splitlist = list( train = data.check( slice.data[ trainidx ], aslist=F ), 
                        test = data.check( slice.data[ -trainidx ], aslist=F ) )
      attr(splitlist, "type") = "fold"
      splitlist #return
    }), paste0("slice", seq_along(data)) )
    attr(splitlist, "type") = "slices"
  } else { #data is a df/trial list
    data = data.split_trials(data)
    splitlist = list(train = data.check( data[ trainidx ], aslist=F ),
                     test = data.check( data[ -trainidx ], aslist=F ))
    attr(splitlist, "type") = "fold"
  }
  attr(splitlist, "test.trials") = which(!seq_along(outcome) %in% trainidx)
  return(splitlist)
}

data.split_folds <- function(data, k=5, idx.out=F) {
  ## fold data for k-fold cross validation
  #INPUT ---
  #data: df, list of trials, slices, subjects
  #k: number of folds
  #idx.out: if True, fold test trial indices are returned
  #         otherwise, data is split into folds
  #NOTE:
  #assumes a single outcome per trial; data is output in df format;
  #if k == 1 only indices are returned (which corresponds to all trials)
  #RETURNS ---
  #if not idx.out, a list with k elements, each containing a train and test df
  #list has type attribute "folds", within each fold the attribute
  #test.trials indicates the trial indices that went into the test set
  #sliced data is simply unsliced and returned back as dfs
  .eval_package("caret")
  data = data.permute_classes(data, shuffle=F) #get type and outcome
  if ( "subjects" %in% attr(data, "type") ) {
    data = lapply(data, data.split_folds, k, split)
    attr(data, "type") = "subjects"
    return(data)
  } 
  #create fold indices based on outcome:
  indxFolds = caret::createFolds(y=attr(data, "outcome"), k=k)
  if ( idx.out || k <= 1 ) return(indxFolds)
  #create data folds:
  data = data.split_trials(data, strip=F)
  folded.data = lapply(indxFolds, function(fold) {
    fold.data = list( train = data.check( data[-fold], aslist=F ),
                      test = data.check( data[fold], aslist=F ) )
    attr(fold.data, "type") = "fold"
    attr(fold.data, "test.trials") = fold
    fold.data #return to lapply
  })
  attr(folded.data, "type") = "folds"
  return(folded.data)
}

data.fit <- function(train, test=NULL, model="LogReg", 
                     center=T, scale=T, ...) {
  ## machine learning procedure: fits a model to a training set 
  ## and uses the model to predict the outcome of a test set
  #INPUT ---
  #train: data set containing training samples, outcome in 2nd col
  #test: data set containing test samples, outcome in 2nd col
  #model: model to fit, defaults to dual L2 logistic regression
  #       available options are: "Reg", "LogReg", "SVM", "LDA", "rLDA"
  #center: if True, data is centered during training
  #scale: if True, data is scaled during training
  #note: the obtained parameters are used to also center/scale the test set
  #RETURNS ---
  #list with elements:
  #AUC: area under the curve, prediction accuracy of the model (if test)
  #predictions: confusion table of predictions and true values (if test)
  #fit: model fit with the respective list structure
  models = c("Reg", "LogReg", "SVM", "LDA", "rLDA") #currently implemented models
  if ( !toupper(model) %in% toupper(models) ) {
    stop( "Undefined model selected. Current options are: ", 
          paste(models, collapse=", "), "." ) 
  }
  train = data.check(train, aslist=F)
  train.outcome = train[,2] #save outcome
  measurements = is.datacol(train)
  if ( toupper(model) != "REG" ) {
    train.outcome = as.factor(train.outcome)
  } else { #regression
    train.outcome = as.numeric(train.outcome)
    measurements[2] = F #set continuous outcome as non-data column
  }
  train = scale(train[, measurements], center, scale) #center/scale
  #save scaling parameters from training set:
  if (center) center = attr(train, "scaled:center")
  if (scale) scale = attr(train, "scaled:scale")
  #fit the model:
  if ( toupper(model) %in% toupper(c("Reg", "LogReg", "SVM")) ) {
    .eval_package("LiblineaR")
    args.in = .eval_ellipsis(LiblineaR::LiblineaR, ...)
    if ( is.null( list(...)$type ) ) { #set type if none specified
      if ( toupper(model) == "LOGREG" ) {
        args.in$type = 7 #dual L2 log reg
      } else if ( toupper(model) == "SVM" ) {
        args.in$type = 1 #dual L2 SVM
      } else {
        args.in$type = 12 #dual L2 reg
      }
    } 
    fit = do.call( LiblineaR::LiblineaR, modifyList(args.in, 
                        list(data=train, target=train.outcome)) )
  } else if ( toupper(model) == "LDA" ) {
    .eval_package("MASS")
    fit = MASS::lda(x=train, grouping=train.outcome)    
  } else if ( toupper(model) == "RLDA" ) { #robust regularized LDA
    .eval_package("rrlda", load=T) #sadly does not work without attaching
    args.in = .eval_ellipsis("rrlda", ...)
    fit = do.call( rrlda, modifyList(args.in, list(x=train, grouping=train.outcome)) )
  }
  #save scaling parameters
  fit$center = center
  fit$scale = scale
  attr(fit, "type") = "fit"
  if ( is.null(test) ) return(fit)
  #if test set is supplied, do prediction
  .eval_package( c("caret", "pROC") )
  test = data.check(test, aslist=F)
  test.outcome =  test[,2]
  test = scale(test[, measurements], center, scale) #center/scale
  #get predictions
  preds <- predict(fit, test)[[1]]
  if ( toupper(model) != "REG" ) {
    preds = as.ordered(preds)
    if ( length(unique(test.outcome)) > 2 ) {
      auc = pROC::multiclass.roc(test.outcome, preds)$auc[1]
    } else {
      auc = pROC::roc(test.outcome, preds)$auc[1]
    }
    out = list( AUC = auc, predictions = caret::confusionMatrix(preds, test.outcome,
                                    dnn=c("Predicted Label", "True Label"))$table )
  } else { #regression
    test.outcome = as.numeric(test.outcome) #making sure outcome is numeric
    out = list( RMSE = sqrt( mean((test.outcome - preds)^2) ), #Root mean squared error
                deviation = summary( test.outcome - preds ) )#distribution of errors
  }
  out$fit = fit
  return(out)
}



createSubjFolds <- function(subnum, k) {
  ## helper to create folds for a list of subject data sets
  ## caret's createFolds does not properly handle this scenario
  ## e.g. repeat createFolds(1:20, k=5) a few times
  #INPUT ---
  #subnum: number of subjects or vector with sub numbers
  #k: number of folds
  #RETURNS: the fold indices
  if (length(subnum) == 1) subs = seq_len(subnum) else subs = seq_along(subnum)
  if (k <= 1) return( list(Fold1 = subs) )
  if (k > length(subs)) stop( "Cannot have more folds than subjects." )
  cuts = sample( cut(subs, breaks=k, labels=F, include.lowest=T) )
  folds = setNames( split(subs, cuts), paste0("Fold",1:k) )
  return(folds)
}
