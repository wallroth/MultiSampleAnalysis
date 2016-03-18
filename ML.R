#ML DT

decode <- function(data, method="within", decompose="none", repetitions=0, permutations=1, nCores=NULL, ...) {
  ## the heart of the decoding package: decode with (repeated) k-fold CV (within/across/between subjects)
  #INPUT ---
  #method: one of c("within","across","between")
  #        within: single subject decoding; each subject data is treated independently from the next
  #        across: data of multiple subjects is analyzed together where folds are stratified for subjects, 
  #                i.e. train and test set contain trials of every subject in equal numbers
  #        between: data of multiple subjects is separated into train and test set
  #                 i.e. train and test set contain trials of different subjects
  #decompose: decomposition method to use, one of c("CSP","SpecCSP","SPoC")
  #repetitions: number of repetitions of the k-fold CV, a value >= 0
  #permutations: number of permutations on each fold, corresponds to the number
  #              of true fits, i.e. a value of 1 means as many permutations as there 
  #              are true fits; 2 means twice as many, etc.; 0 = no permutations
  #nCores: number of CPU cores to use for parallelization
  #        if NULL, automatic selection; if 1, sequential execution
  #        if an empty list, an externally registered cluster will be used
  #RETURNS ---
  #a DT showing the performance scores per run/fold/slice/label (whichever is applicable)
  .eval_package("foreach", load=T)
  method = tolower(method); decompose = toupper(decompose)
  if ( !method %in% c("within","across","between") ) stop( "Undefined decoding method. Options are 'within', 'across' and 'between'." )
  if ( !decompose %in% c("NONE","CSP","SPECCSP","SPOC") ) stop( "Undefined decomposition method. Options are 'CSP', 'SpecCSP' and 'SPoC'." )
  data = data.check(data)
  #retrieve arguments for various functions later down the line (exclude function name and data from evaluation)
  args.in = lapply( as.list(match.call())[-c(1,2)], eval ) #eval input before foreach in case variable names were passed
  args = do.call( .eval_ellipsis, modifyList(args.in, list(funStr=c("fold.idx","slice.idx","data.permute_classes","data.project","model.fit"))) )
  args$pass = args.in[ !names(args.in) %in% c( names(formals(decode)), do.call(c, lapply(args, names)) ) ] #additional arguments
  ## evaluate model specification once to reduce overhead from multiple model.fit calls ##
  models = c("Reg", "LogReg", "SVM", "LDA", "SDA") #currently implemented models
  args$model.fit$model = toupper(args$model.fit$model)
  if ( !args$model.fit$model %in% toupper(models) ) {
    stop( "Undefined model selected. Current options are: ", paste(models, collapse=", "), "." ) 
  }
  args$model.fit[c("center","scale")] = lapply(args$model.fit[c("center","scale")], eval) #evaluate logicals
  if (args$model.fit$model != "REG") { #classification
    .eval_package("pROC")
    args$model.fit$multiClass = data[, length(unique(outcome))] > 2
  }
  args$model.fit$nm = setdiff( names(data), key(data) )
  if ( args$model.fit$model %in% c("REG", "LOGREG", "SVM") ) {
    .eval_package("LiblineaR")
    if ( is.null(args.in$type) ) { #set type if none specified
      if ( args$model.fit$model == "LOGREG" ) args$pass$type = 0 #primal L2 log reg
      else if ( args$model.fit$model == "SVM" ) args$pass$type = 2 #primal L2 SVM
      else args$pass$type = 11 #primal L2 reg
    }
    funargs = as.list( formals(LiblineaR::LiblineaR) )
  } else if ( args$model.fit$model == "LDA" ) {
    .eval_package("MASS")
    funargs = list()
  } else if ( args$model.fit$model == "SDA" ) { #shrinkage LDA
    .eval_package("sda")
    funargs = modifyList( as.list(formals(sda::sda)), list(verbose=F) )
  }
  args$model.fit$fitfun = modifyList( funargs, args$pass[ names(args$pass) %in% names(funargs) ], keep.null=T )
  ## end model evaluation ##
  if (method == "between") args$fold.idx$fold.type = "subjects"
  args$decomp = decompose != "NONE"
  args$within = method == "within"
  if (args$decomp) { #a decomposition method is set
    args$data.project[c("approx","logscale")] = lapply(args$data.project[c("approx","logscale")], eval) #evaluate logicals
    args$pass$silent = T #suppress messages from decompose methods
    args$decomp.method = paste0("decompose.", grep(decompose, c("CSP","SpecCSP","SPoC"), ignore.case=T, value=T)[1])
    if ( args$data.project$approx && args$slice.idx$window == 1 ) {
      cat( "With a window size of 1 approximation of trial variance is not possible; setting to FALSE.\n" )
      args$data.project$approx = F
    }
    if ( args$data.project$approx && args$data.project$logscale && is.null(args.in$scale) ) {
      cat( "Deactivating scaling via SD because log-transformation of trial variance is set.\n" )
      args$model.fit$scale = F
    }
  }
  if ( nchar(args$slice.idx$window) > 0 ) { #get slice indices
    allequal = data[, .N, by=.(subject,trial)][, all(diff(N) == 0)] #do all trials have the same number of samples
    if (!allequal) stop( "Assumption of non-varying sample numbering across subjects is violated." )
    minmax = data[c(1,.N), sample] #get first/last sample number (sorted = min/max)
    slices = do.call( slice.idx, modifyList(args$slice.idx, list(n=minmax)) ) #one slice setting is applied to all subjects
  } else { #no slicing, every sample is used
    slices = matrix( data[, unique(sample)], ncol=1 )
  }
  repetitions = repetitions + 1 #a single run equals 1 'repetition'
  cl = parBackend(nCores, required=ncol(slices)) #initialize parallel backend (if not already)
  ## fit function for a single slice ##
  fit.slice <- function(train.slice, test.slice, args) {
    if (!args$within) { #collapse data across/between subjects
      return( .model.fit(train=train.slice, test=test.slice, args=args$model.fit) )
    } else { #keep subject data separate
      perfs = lapply(train.slice[, unique(subject)], function(subj) {
        .model.fit(train=train.slice[.(subj)], test=test.slice[.(subj)], args=args$model.fit)
      })
      return( set( rbindlist(perfs), j="subject", value=train.slice[, unique(subject)] ) )
    }
  }
  ### START FITTING ###
  run.performance = list() #collect run results
  for ( run in seq_len(repetitions) ) {
    runtime = proc.time()[3]
    if (repetitions > 1) cat("Run",run,"\n")
    #obtain the folds for the CV scheme:
    folds = do.call( fold.idx, modifyList(args$fold.idx, list(data=data)) )
    fold.performance = list() #collect fold results
    for ( foldnum in folds[, unique(fold)] ) { #iterate the folds
      cat("Fold",foldnum)
      #subset data into test/training sets:
      if (method == "between") test = data[ folds[fold == foldnum, .(subject)] ]
      else test = data[ folds[fold == foldnum, .(subject, trial)] ] #within/across
      train = data[!test] #mutually exclusive data sets irrespective of method
      if (nrow(train) == 0) train = test #no CV
      #fit a model with unchanged (true) labels:
      performance = foreach(i=1:ncol(slices), .combine=function(...) rbindlist(list(...), idcol="slice"), 
                            .multicombine=T, .maxcombine=ncol(slices)+1) %dopar%
        fit.slice(train[ sample %in% slices[,i] ], test[ sample %in% slices[,i] ], args)
      if ( truelength(performance) == 0 ) alloc.col(performance, 10) #fix to 0 allocated columns error when parallelized without slicing
      set(performance, j="label", value="true") #add label info
      #fit model(s) with shuffled labels:
      for ( perm in seq_len(permutations) ) { #repeat for number of permutations
        trainP = data.permute_classes(train, tol=args$data.permute_classes$tol) #permute training set labels
        testP = data.permute_classes(test, tol=args$data.permute_classes$tol) #permute test set labels
        perform.rand = foreach(i=1:ncol(slices), .combine=function(...) rbindlist(list(...), idcol="slice"), 
                               .multicombine=T, .maxcombine=ncol(slices)+1) %dopar%
          fit.slice(trainP[ sample %in% slices[,i] ], testP[ sample %in% slices[,i] ], args)
        if ( truelength(perform.rand) == 0 ) alloc.col(perform.rand, 10) #fix to 0 allocated columns error when parallelized without slicing
        performance = rbind(performance, set(perform.rand, j="label", value="random")) #bind to existing performance DT
      }
      fold.performance[[foldnum]] = performance #add to fold list
      cat(" | Time elapsed (mins):", round((proc.time()[3]-runtime)/60,2), "\n")
    } #end of fold
    run.performance[[run]] = rbindlist(fold.performance, idcol=if (length(fold.performance)>1) "fold") #bind to one DT
  } #end of run
  ### END FITTING ###
  performance = rbindlist(run.performance, idcol=if (repetitions>1) "run") #bind to one DT
  keys = intersect( c("subject","slice","run","label","fold"), names(performance) ) #key existing info columns
  setcolorder( performance, union(keys, names(performance)) ) #re-arrange column order
  setkeyv(performance[, label := factor(label, levels=c("true","random"))], keys) #sort by key order
  parBackend(cl) #close backend (if initialized by function) and print final time
  return( setattr(performance, "call", match.call())[] )
}  

model.fit <- function(train, test=NULL, model="LogReg", center=T, scale=T, ...) {
  ## machine learning procedure: fits a model to a training set 
  ## and uses the model to predict the outcome of a test set
  #INPUT ---
  #train: data set containing training samples
  #test: data set containing test samples
  #model: model to fit, defaults to primal L2 logistic regression
  #       available options are: "Reg", "LogReg", "SVM", "LDA", "SDA"
  #center: if True, data is centered during training
  #scale: if True, data is scaled during training
  #note: the obtained parameters are used to also center/scale the test set
  #RETURNS ---
  #list with elements:
  #fit: model fit with the respective list structure
  #if test set was supplied, attribute performance
  models = c("Reg", "LogReg", "SVM", "LDA", "SDA") #currently implemented models
  model = toupper(model)
  if ( !model %in% toupper(models) ) {
    stop( "Undefined model selected. Current options are: ", paste(models, collapse=", "), "." ) 
  }
  nm = setdiff( names(train), key(train) )
  train.outcome = train[, outcome]
  train = scale(train[, .SD, .SDcols=nm], center, scale)
  #save scaling parameters from training set:
  if (center) center = attr(train, "scaled:center")
  if (scale) scale = attr(train, "scaled:scale")
  #fit the model:
  if ( model %in% toupper(c("Reg", "LogReg", "SVM")) ) {
    .eval_package("LiblineaR")
    args.in = modifyList( as.list(formals(LiblineaR::LiblineaR)), list(...), keep.null=T )
    if ( is.null( list(...)$type ) ) { #set type if none specified
      if ( model == "LOGREG" ) {
        args.in$type = 0 #primal L2 log reg
      } else if ( model == "SVM" ) {
        args.in$type = 2 #primal L2 SVM
      } else {
        args.in$type = 11 #primal L2 reg
      }
    }
    fit = do.call( LiblineaR::LiblineaR, modifyList(args.in, list(data=train, target=train.outcome)) )
  } else if ( model == "LDA" ) {
    .eval_package("MASS")
    fit = MASS::lda(x=train, grouping=train.outcome)    
  } else if ( model == "SDA" ) { #shrinkage LDA
    .eval_package("sda")
    args.in = modifyList( as.list(formals(sda::sda)), list(...), keep.null=T )
    args.in = args.in[ names(args.in) %in% names(formals(sda::sda)) ]
    fit = do.call( sda::sda, modifyList(args.in, list(Xtrain=train, L=train.outcome, verbose=F)) )
  }
  #save scaling parameters
  setattr(fit, "center", center)
  setattr(fit, "scale", scale)
  if ( is.null(test) ) return(fit)
  #if test set is supplied, do prediction
  .eval_package("pROC")
  test.outcome =  test[, outcome]
  test = scale(test[, .SD, .SDcols=nm], center, scale)
  #get predictions
  preds <- predict(fit, test, verbose=F)[[1]]
  if ( model != "REG" ) {
    preds = as.ordered(preds)
    if ( length(unique(test.outcome)) > 2 ) {
      auc = pROC::multiclass.roc(test.outcome, preds, direction="<", algorithm=3)$auc[1]
    } else {
      auc = pROC::auc(test.outcome, preds, direction="<", algorithm=3)[1]
    }
    CM = table(preds, test.outcome) #confusion matrix
    setattr( fit, "performance", data.table(AUC = auc, correct = sum(diag(CM)), incorrect = sum(CM[upper.tri(CM)], CM[lower.tri(CM)])) )
  } else { #regression
    residuals = test.outcome - preds
    setattr( fit, "performance", data.table(RMSE = sqrt(mean(residuals^2)), R2 = 1-var(residuals)/var(test.outcome)) )
  }
  return(fit)
}

fold.idx <- function(data, k=NULL, n.classes=NULL, fold.type="trials") {
  #k: number of folds, if between 0 and 1 regarded as percentage of data for test set; if NULL: LOO
  #n.classes: if outcome is a continous variable, n.classes should bet set to set the number of categories
  #           if outcome is a float, n.classes will default to 5 if not specified
  #fold.type: trials|subjects - specifies if folds are built with trials or subject ids
  fold.type = tolower(fold.type)
  if ( !fold.type %in% c("trials","subjects") ) stop( "Undefined fold type. Options are 'trials' or 'subjects'." )
  data = data.check(data)
  
  if (fold.type == "trials") { #create for each subject folds with trials balanced for the classes in outcome
    dat = data[, .(outcome = outcome[1]), by=.(subject,trial)] #reduce to necessary info
    #if no k specified: Leave-One-Out based on minimum outcome occurence across all subjects
    if ( is.null(k) || k <= 0 ) k = dat[, .N, by=.(subject,outcome)][, min(N)] 
    folds = dat[, {
      if ( (is.numeric(outcome) && !is.null(n.classes)) || is.float(outcome) ) { #create grouping for numeric vector
        if (is.null(n.classes)) n.classes = 5 #if float and no n specified
        #create classes according to n quantiles
        outcome = cut(outcome, breaks=quantile(outcome, probs=seq(0,1,length.out=n.classes)), labels=F, include.lowest=T)
      }
      class.tab = table(outcome) #unique classes and counts
      if (k <= 1) { #percentage of trials per outcome which go into test set
        test = c() #trial vector
        for (i in 1:length(class.tab)) {
          test = c(test, sample( trial[ outcome==names(class.tab)[i] ], size=round(class.tab[i]*k) )) #sample trial index  
        }
        .(fold = 1, trial = test, outcome = outcome[trial %in% test]) #return
      } else { #fold split
        if ( any(class.tab%/%k < 1) ) stop( "Cannot have more folds than trials per class." )
        folds = rep(0, length(outcome)) #fold vector
        for (i in 1:length(class.tab)) {
          num.reps = class.tab[i]%/%k #num times of class per fold
          folds[ outcome==names(class.tab)[i] ] = c( sample(rep(1:k, num.reps)), sample(1:k, class.tab[i]%%k) ) #append if any are left
        }
        .(fold = folds, trial = trial, outcome = outcome) #return
      }
    }, by=subject]
    return( setkeyv(folds, c("subject","fold","trial"))[] )
    
  } else { #create folds with all trials from a subset of the subjects
    folds = data[, .(subject = unique(subject))][, {
    if ( is.null(k) || k <= 0 ) k = length(subject) #Leave-One-Subject-Out
    if (k <= 1) { #percentage of subjects which go into test set
      .(fold = 1, subject = sample( subject, size=round(length(subject)*k) )) #return
    } else { #fold split
      if (k > length(subject)) stop( "Cannot have more folds than subjects." )
      .(fold = sample(cut(subject, breaks=k, labels=F, include.lowest=T)), subject = subject)
    } }]
    return( setkeyv(folds, c("fold", "subject"))[] )
  }
}

data.permute_classes <- function(data, tol=1) {
  ## extracts class labels from data and permutes if requested
  #tol: the permutation success is verified by checking for each class if half the labels +/- tol * that
  #     have switched their position; e.g. with 100 trials for a class, 50 of those have to be switched for perfect permutation
  #     if tol = 0, exactly 50 have to be switched; if tol = 1, it can be any number (1:100)
  #     with tol at .2, the tolerance range will be from 30:70, i.e. (50-100*.2):(50+100*.2)
  #RETURNS: data with shuffled outcome labels
  data = data.check( copy(data) )
  keys = key(data)
  permuted = data[, .(outcome = outcome[1]), by=.(subject,trial)][, { 
    classes = table(outcome)
    success = F
    while ( !success ) {
      tmp = sample(outcome) #shuffle
      #check how many class labels switched positions
      success = all( sapply(1:length(classes), function(i) { #iterate class labels
        #range in which permutation was successful:
        range = floor(classes[i]/2 - classes[i]/2*tol):ceiling(classes[i]/2 + classes[i]/2*tol)
        sum( which(outcome == names(classes)[i]) %in% which(tmp == names(classes)[i]) ) %in% range
      }) )
      #make sure even with tol=1 the shuffle is not identical to the true labels:
      if ( identical(tmp, outcome) ) success=F
    }
    .(outcome = tmp) }, by=subject]
  return( setkeyv(data[, outcome := permuted[.GRP,outcome], by=.(subject,trial)], keys)[] )
}


#### private

.model.fit <- function(train, test, args) {
  ## model fit with less overhead for repeated execution
  train.outcome = train[, outcome]
  train = scale(train[, .SD, .SDcols=args$nm], args$center, args$scale)
  #save scaling parameters from training set:
  if (args$center) center = attr(train, "scaled:center")
  else center = F
  if (args$scale) scale = attr(train, "scaled:scale")
  else scale = F
  #fit the model:
  if ( args$model %in% c("REG", "LOGREG", "SVM") ) {
    fit = do.call( LiblineaR::LiblineaR, modifyList(args$fitfun, list(data=train, target=train.outcome)) )
  } else if ( args$model == "LDA" ) {
    fit = MASS::lda(x=train, grouping=train.outcome)    
  } else if ( args$model == "SDA" ) { #shrinkage LDA
    fit = do.call( sda::sda, modifyList(args$fitfun, list(Xtrain=train, L=train.outcome)) )
  }
  test.outcome =  test[, outcome]
  test = scale(test[, .SD, .SDcols=args$nm], center, scale)
  #get predictions
  preds <- predict(fit, test, verbose=F)[[1]]
  if ( args$model != "REG" ) {
    preds = as.ordered(preds)
    if (args$multiClass) {
      auc = pROC::multiclass.roc(test.outcome, preds, direction="<", algorithm=3)$auc[1]
    } else {
      auc = pROC::auc(test.outcome, preds, direction="<", algorithm=3)[1]
    }
    CM = table(preds, test.outcome) #confusion matrix
    return( data.table(AUC = auc, correct = sum(diag(CM)), incorrect = sum(CM[upper.tri(CM)], CM[lower.tri(CM)])) )
  } else { #regression
    residuals = test.outcome - preds
    return( data.table(RMSE = sqrt(mean(residuals^2)), R2 = 1-var(residuals)/var(test.outcome)) )
  }
}
