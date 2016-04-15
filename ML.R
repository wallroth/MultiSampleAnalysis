#ML DT

decode <- function(data, method="within", decompose="none", GAT=F, repetitions=0, permutations=1, CV=NULL, nCores=NULL, ...) {
  ## the heart of the decoding package: decode with (repeated) k-fold CV (within/across/between subjects)
  #INPUT ---
  #method: one of c("within","across","between")
  #        within: single subject decoding; each subject data is treated independently from the next
  #        across: data of multiple subjects is analyzed together where folds are stratified for subjects, 
  #                i.e. train and test set contain trials of every subject in equal numbers
  #        between: data of multiple subjects is separated into train and test set
  #                 i.e. train and test set contain trials of different subjects
  #decompose: decomposition method to use, one of c("CSP","SpecCSP","SPoC")
  #GAT: if True, extend to "Generalization across time" procedure: 
  #     additionally to the standard decoding approach, the training time model fits are 
  #     applied to all other slices (test times) to generalize the information at one time point to others;
  #     the GAT diagonal (train time = test time) corresponds to the default decoding scheme (i.e. GAT=F)
  #repetitions: number of repetitions of the k-fold CV, a value >= 0
  #permutations: number of permutations on each fold, corresponds to the number
  #              of true fits, i.e. a value of 1 means as many permutations as there 
  #              are true fits; 2 means twice as many, etc.; 0 = no permutations
  #nCores: number of CPU cores to use for parallelization
  #        if NULL, automatic selection; if 1, sequential execution
  #        if an empty list, an externally registered cluster will be used
  #RETURNS ---
  #a DT showing the performance scores per subject/slice/test/label/run/fold (whichever is applicable)
  .eval_package(c("foreach","fastmatch"), load=T)
  "%in%" <- function(x, table) fmatch(x, table, nomatch = 0) > 0 #overwrite match with fmatch
  method = tolower(method); decompose = toupper(decompose)
  if ( !method %in% c("within","across","between") ) stop( "Undefined decoding method. Options are 'within', 'across' and 'between'." )
  if ( !decompose %in% c("NONE","CSP","SPECCSP","SPOC") ) stop( "Undefined decomposition method. Options are 'CSP', 'SpecCSP' and 'SPoC'." )
  if ( !is.null(CV) && repetitions > 0 ) {
    cat( "Ignoring the requested number of repetitions because a predefined CV scheme is supplied.\n" )
    repetitions = 0
  }
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
  nm = setdiff( names(data), key(data) )
  if (args$model.fit$model != "REG") { #classification
    .eval_package("pROC")
    args$model.fit$multiClass = data[, uniqueN(outcome) > 2]
    if (!is.ordered(data$outcome)) setkeyv( data[, outcome:=ordered(outcome)], setdiff(names(data), nm) )
  }
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
  if ( method == "between" && is.null(CV) ) args$fold.idx$fold.type = "subjects"
  else if ( method == "across" && is.null(args$fold.idx$n.folds) && is.null(CV) ) { #LOO: minimum occurence across subjects
    n = args$fold.idx$n.classes
    if ( is.null(n) && !is.float(data$outcome) ) { #outcome is grouped by its value
      args$fold.idx$n.folds = data[sample == sample[1], outcome, by=.(subject,trial)][, .N, by=.(subject,outcome)][, min(N)]  
    } else { #outcome is grouped by quantiles
      if (is.null(n)) n = 4 #default 5 quantiles for float numbers (= 4 groups)
      args$fold.idx$n.folds = data[, uniqueN(trial)/n, by=subject][, min(floor(V1))] #number of values distributed evenly across n groups
    }
  }
  if (permutations > 0) {
    pkeys = c("subject","trial","sample") #keys to resort permuted labels with
    ptol = args$data.permute_classes$tol
    train.allow = eval(args$data.permute_classes$allow.identical)
    #if not specified check if LOO and if so, allow identicals (otherwise 100% switching is enforced)
    test.allow = ifelse( is.null(args.in$allow.identical) && is.null(args$fold.idx$n.folds) &&
         (is.null(CV) || CV[, .(N=uniqueN(trial)), by=.(subject,fold)][, min(N) <= 2]), T, train.allow )
  }
  args$DECOMPOSE = decompose != "NONE"
  args$within = method == "within"
  if (args$DECOMPOSE) { #a decomposition method is set
    args$data.project[c("approx","logscale")] = lapply(args$data.project[c("approx","logscale")], eval) #evaluate logicals
    args$data.project$npattern = ifelse(is.null(args.in$npattern), 3, args.in$npattern)
    args$data.project$nm = nm
    args$data.project$nmV = paste0("V",1:length(nm))
    args$DECOMPOSE.method = paste0("decompose.", grep(decompose, c("CSP","SpecCSP","SPoC"), ignore.case=T, value=T)[1])
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
    if (GAT) stop( "Data must be sliced for GAT. Please define a window size, cf. slice.idx" )
    slices = matrix( data[, unique(sample)], ncol=1 )
  }
  repetitions = repetitions + 1 #a single run equals 1 'repetition'
  dataMat = unname( as.matrix(data[, .SD, .SDcols=nm]) ) #transform data to matrix (measurements only)
  dataInfo = data[, .SD, .SDcols=!nm] #info columns DT with subset information
  slicelist = lapply(1:ncol(slices), function(i) { #pre-create slice idx list
    idx = data[ sample %in% slices[,i], which=T ]
    if (args$within) split(idx, dataInfo$subject[idx]) #nest slices for subjects
    else idx
  })
  subIDs = unique(dataInfo$subject)
  cl = parBackend(nCores, required=ncol(slices)) #initialize parallel backend (if not already)
  
  ## fit function for a single slice ##
  fit.slice <- function(slice) {
#     if (args$DECOMPOSE) {
#       features = .DECOMPOSE(train.slice, test.slice, args=args$data.project)
#       test.slice = features[ setkeyv(test.slice[, .(trial=unique(trial)), by=subject], c("subject","trial")) ]
#       train.slice = features[!test.slice]
#       if (nrow(train.slice) == 0) train.slice = test.slice #no CV
#     }
    if (!args$within) { #collapse data across/between subjects
      return( setDT(.PREDICT(.FIT(slice$train, cvoutcome$train.outcome, args$model.fit),
                             slice$test, cvoutcome$test.outcome, args$model.fit)) )
    } else { #slice is nested for subjects to keep subject data separate
      perfs = lapply(seq_along(testIDs), function(sub) {
        .PREDICT(.FIT(slice[[sub]]$train, cvoutcome[[sub]]$train.outcome, args$model.fit),
                 slice[[sub]]$test, cvoutcome[[sub]]$test.outcome, args$model.fit)
      })
      return( set(rbindlist(perfs), j="subject", value=testIDs) )
    }
  }
  fit.GAT <- function(slice) {
    if (!args$within) { #collapse data across/between subjects
      fit = .FIT(slice, cvoutcome$train.outcome, args$model.fit) #train at a slice
      return( set( rbindlist( lapply(testdata, function(test) { #predict at all slices
        .PREDICT(fit, test, cvoutcome$test.outcome, args$model.fit)
      })), j="test", value=seq_along(testdata)) )
    } else { #slice is nested for subjects to keep subject data separate
      return( rbindlist(lapply(seq_along(testIDs), function(sub) { #iterate subjects
        fit = .FIT(slice[[sub]], cvoutcome[[sub]]$train.outcome, args$model.fit) #train at a slice
        set( set( rbindlist( lapply(testdata[[sub]], function(test) { #iterate slices for prediction
          .PREDICT(fit, test, cvoutcome[[sub]]$test.outcome, args$model.fit) 
        }) ), j="test", value=seq_along(testdata[[sub]]) ), j="subject", value=testIDs[sub] )
      })) )
    }
  }
  
  ### START FITTING ###
  run.performance = list() #collect run results
  for ( run in seq_len(repetitions) ) {
    runtime = proc.time()[3]
    if (repetitions > 1) cat("Run",run,"\n")
    if (is.null(CV)) { #obtain the folds for the CV scheme:
      folds = do.call( fold.idx, modifyList(args$fold.idx, list(data=data)) )
    } else {
      folds = CV
    }
    fold.performance = list() #collect fold results
    for ( foldnum in folds[, unique(fold)] ) { #iterate the folds
      cat("Fold",foldnum)
      #subset data into test/training sets:
      testIdx = dataInfo[ folds[fold == foldnum, .(subject, trial)], nomatch=0, which=T ]
      trainIdx = dataInfo[!testIdx, which=T] #mutually exclusive data sets
      if (length(trainIdx) == 0) trainIdx = testIdx #no CV
      if (args$within) { #subjects with fewer trials may not be part of the fold
        testIDs = unique(dataInfo$subject[testIdx]) #ID values
        posIDs = which(subIDs %in% testIDs) #ID positions in slicelist
      }
      if (GAT) { #create empty list structure for test data (every worker will require the full test set)
        if (args$within) testdata = list()[subIDs] #main structure: subjects, nested for slices
        else testdata = list()[seq_along(slicelist)] #slice list
      }
      #create the sliced train/test sets outside foreach to reduce memory demand
      cvdata = lapply(seq_along(slicelist), function(slice) { #iterate the slices
        if (args$within) { #nested for subjects
          lapply(posIDs, function(sub) { #index into respective subject element
            idx = slicelist[[slice]][[sub]] #grab indices for subject at slice time 
            if (GAT) {
              testdata[[sub]][[slice]] <<- dataMat[ idx[idx %in% testIdx], ] #test set is used at every slice
              dataMat[ idx[idx %in% trainIdx], ] #only train set is specific to slice
            } else { #train and test are specific to every slice
              list(train = dataMat[ idx[idx %in% trainIdx], ], test = dataMat[ idx[idx %in% testIdx], ])
            }
          })
        } else { #across/between
          idx = slicelist[[slice]]
          if (GAT) {
            testdata[[slice]] <<- dataMat[ idx[idx %in% testIdx], ]
            dataMat[ idx[idx %in% trainIdx], ]
          } else {
            list(train = dataMat[ idx[idx %in% trainIdx], ], test = dataMat[ idx[idx %in% testIdx], ])
          }
        }
      })
      #remove NULL list elements from testdata so positions are identical across lists
      if (GAT && args$within && length(posIDs) < length(subIDs)) testdata = testdata[ !sapply(testdata, is.null) ]
      #create outcome list separately as it remains identical across slices
      if (args$within) {
        cvoutcome = lapply(slicelist[[1]][posIDs], function(idx) {
          list(train.outcome = dataInfo$outcome[ idx[idx %in% trainIdx] ], 
               test.outcome = dataInfo$outcome[ idx[idx %in% testIdx] ])
        })
      } else {
        idx = slicelist[[1]]
        cvoutcome = list(train.outcome = dataInfo$outcome[ idx[idx %in% trainIdx] ], 
                         test.outcome = dataInfo$outcome[ idx[idx %in% testIdx] ])
      }
      #fit a model with unchanged (true) labels:
      performance = foreach(slice=cvdata, .combine=function(...) rbindlist(list(...), idcol="slice"),
                            .multicombine=T, .maxcombine=length(cvdata)+1) %dopar%
                            {
                              if (GAT) fit.GAT(slice) else fit.slice(slice)
                            }
      if ( truelength(performance) == 0 ) alloc.col(performance, 10) #fix to 0 allocated columns error when parallelized without slicing
      set(performance, j="label", value="true") #add label info
      
      #fit model(s) with shuffled labels:
      for ( perm in seq_len(permutations) ) { #repeat for number of permutations
        #permute training and test set labels separately and re-order by setting the key (so indices remain valid)
        dataInfoP = setkeyv(rbind( data.permute_classes(dataInfo[trainIdx], ptol, train.allow), 
                                   data.permute_classes(dataInfo[testIdx], ptol, test.allow) ), pkeys)
        #update outcome
        if (args$within) {
          cvoutcome = lapply(slicelist[[1]][posIDs], function(idx) {
            list(train.outcome = dataInfoP$outcome[ idx[idx %in% trainIdx] ], 
                 test.outcome = dataInfoP$outcome[ idx[idx %in% testIdx] ])
          })
        } else {
          idx = slicelist[[1]]
          cvoutcome = list(train.outcome = dataInfoP$outcome[ idx[idx %in% trainIdx] ], 
                           test.outcome = dataInfoP$outcome[ idx[idx %in% testIdx] ])
        }
        perform.rand = foreach(slice=cvdata, .combine=function(...) rbindlist(list(...), idcol="slice"), 
                               .multicombine=T, .maxcombine=length(cvdata)+1) %dopar%
                               {
                                 if (GAT) fit.GAT(slice) else fit.slice(slice)             
                               }
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
  keys = intersect( c("subject","slice","test","run","label","fold"), names(performance) ) #key existing info columns
  setcolorder( performance, union(keys, names(performance)) ) #re-arrange column order
  setkeyv(performance[, label := factor(label, levels=c("true","random"))], keys) #sort by key order
  parBackend(cl) #close backend (if initialized by function) and print final time
  return( setattr(performance, "call", match.call())[] )
}  

decode.test <- function(result, test="wilcox", one.sided=T, dv="AUC", within.var="label", within.null="random",
                        between="slice", id="subject", adjust="fdr", alpha=0.05, bootstrap=F, iterations=1000) {
  ## significance testing and permutation confidence intervals of decoding result
  #INPUT ---
  #result: df with columns specified in dv, id, between, within.var
  #test: wilcox, t./-test, none
  #      NOTE: care that fold splits or permutations within folds were not repeated, otherwise
  #            the variance can be made arbitrarily small which will then lead to grossly underestimated p-values
  #one.sided: if True, alternative hypothesis is set to "greater", otherwise "two.sided"
  #dv: dependent variable, column with the performance score; required
  #within.var: pair to compare, e.g. column with the labels (true vs. random); required
  #within.null: within.var's value that represents the null distribution (used to compute the CIs with)
  #between: variable that groups the dv into separate segments (e.g. slices); can be left empty
  #id: variable that represents the independent parts within the comparisons
  #    most likely one of 'subject' or 'fold' (for p-values in single subject case)
  #adjust: correction method for the p-values, defaults to FDR with dependencies; can be "none"
  #alpha: significance threshold to test p-values against, CIs are computed for this level; can be multiple values
  #bootstrap: if True, resamples the dv of the within.null group with replacement per id and computes
  #           the CI on each iteration, of which the mean is returned and finally aggregaed over id
  #           if False, the CIs are computed directly on the dv values of within.null  
  #iterations: number of times to repeat the bootstrap 
  #RETURNS ---
  #a summary df with means, CIs, p-values, significance indication
  test = tolower(test)
  if (!test %in% c("wilcox", "t-test", "t.test", "none")) test = "none"
  if (!is.data.table(result)) result = setDT(copy(result), key=c(id, between, within.var))
  #rename (or add if missing) columns
  result = setnames(copy(result), old=c(dv,within.var), new=c("dv","wv"))
  if (!is.null(id) && nzchar(id)) setnames(result, old=id, new="id") else result[, id:=0]
  if (!is.null(between) && nzchar(between)) setnames(result, old=between, new="bv") else result[, bv:=0]
  alpha = sort(alpha, decreasing=T) #in case of multiple values
  aggres = result[, .(dv=mean(dv)), by=.(id,bv,wv)] #aggregate dv over the groups
  if (test != "none") { #compute p-values
    hypothesis = ifelse(one.sided, "greater", "two.sided")
    ptab = aggres[, {
      if (test == "wilcox") .(p = wilcox.test(dv ~ wv, paired=T, alternative=hypothesis)$p.value)
      else .(p = t.test(dv ~ wv, paired=T, alternative=hypothesis)$p.value)
    }, by=bv]
    if (tolower(adjust) != "none") ptab[, p.adj:=p.adjust(p, method=adjust)]
  }
  CI = result[wv == within.null, { 
    if (bootstrap) { #bootstrapped CI of the null distribution
      x = replicate(iterations, quantile(sample(dv, replace=T), 1-alpha)) #repeated sampling with replacement
      if (length(alpha)>1) split(rowMeans(x), 1:length(alpha)) else .("1"=mean(x))
    } else { #CI without bootstrapping
      split(quantile(dv, 1-alpha), 1:length(alpha))
    }
  }, by=.(id,bv)][, split(colMeans(.SD),1:length(alpha)), .SDcols=as.character(1:length(alpha)), by=bv]
  setnames(CI, old=as.character(1:length(alpha)), new=paste0("CI",100*(1-alpha)) ) #rename columns
  #generate output, merge with CI
  output = aggres[, .(mean=mean(dv[wv!=within.null]), effect = mean(dv[wv!=within.null]-dv[wv==within.null])), by=bv][CI, on="bv"]
  if (test != "none") { #if p-values were computed, merge and indicate significance (using sum in case of multiple alphas)
    output = output[ptab, on="bv"][, sign := paste0( rep("*", ifelse(tolower(adjust)!="none", sum(p.adj < alpha, na.rm=T), 
                                                                     sum(p < alpha, na.rm=T))), collapse=""), by=bv]
  }
  if (is.null(between) || !nzchar(between)) set(output, j="bv", value=NULL) #remove obsolete bv column
  else setnames(output, old="bv", new=between) #change back to original name
  return(output[])
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
    resids = test.outcome - preds
    R2 = sum(resids^2)/sum((test.outcome-mean(test.outcome))^2) #SSE/SST
    #1 - [(n-1)/(n-p)] * (SSE/SST) #adjustment doesn't work when n < p?
    setattr( fit, "performance", data.table(RMSE = sqrt(mean(resids^2)), R2=1-R2, N=length(resids)) )
  }
  return(fit)
}

fold.idx <- function(data, n.folds=NULL, n.classes=NULL, fold.type="trials") {
  #n.folds: number of folds, if between 0 and 1 regarded as percentage of data for test set; if NULL: LOO
  #n.classes: if outcome is a continous variable, n.classes should bet set to set the number of categories
  #           if outcome is a float, n.classes will default to 5 (=4 groups) if not specified
  #fold.type: trials|subjects - specifies if folds are built with trials or subject ids
  fold.type = tolower(fold.type)
  if ( !fold.type %in% c("trials","subjects") ) stop( "Undefined fold type. Options are 'trials' or 'subjects'." )
  data = data.check(data)

  if (fold.type == "trials") { #create for each subject folds with trials balanced for the classes in outcome
    folds = data[sample == sample[1], {
      if ( (is.numeric(outcome) && !is.null(n.classes)) || is.float(outcome) ) { #create grouping for numeric vector
        if (is.null(n.classes)) n.classes = 5 #if float and no n specified
        #create classes according to n quantiles, results in n-1 groups
        outcome = cut(outcome, breaks=quantile(outcome, probs=seq(0,1,length.out=n.classes)), labels=F, include.lowest=T)
      }
      class.tab = table(outcome) #unique classes and counts
      class.tab = class.tab[class.tab > 0] #in case of undropped factor levels
      #if n.folds not specified: Leave-One-Out based on minimum outcome occurence
      if ( is.null(n.folds) || n.folds <= 0 ) k = min(class.tab)
      else k = n.folds
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
  } else { #create folds with all trials from a subset of the subjects
    folds = setkeyv(data[, .(subject = unique(subject))][, {
    if ( is.null(n.folds) || n.folds <= 0 ) n.folds = length(subject) #Leave-One-Subject-Out
    if (n.folds <= 1) { #percentage of subjects which go into test set
      .(fold = 1, subject = sample( subject, size=round(length(subject)*n.folds) )) #return
    } else { #fold split
      if (n.folds > length(subject)) stop( "Cannot have more folds than subjects." )
      .(fold = sample(cut(subject, breaks=n.folds, labels=F, include.lowest=T)), subject = subject)
    } }], c("fold", "subject"))[data[, .(trial=unique(trial)), by=subject]]
  }
  return( setkeyv(folds, c("subject","fold","trial"))[] )
}

data.permute_classes <- function(data, tol=1, allow.identical=F) {
  ## extracts class labels from data and permutes if requested
  #tol: the permutation success is verified by checking for each class if half the labels +/- tol * that
  #     have switched their position; e.g. with 100 trials for a class, 50 of those have to be switched for perfect permutation
  #     if tol = 0, exactly 50 have to be switched; if tol = 1, it can be any number (1:100)
  #     with tol at .2, the tolerance range will be from 30:70, i.e. (50-100*.2):(50+100*.2)
  #allow.identical: if False, the input will be excluded from possible permutations
  #RETURNS: data with shuffled outcome labels
  data = data.check( copy(data) )
  keys = key(data)
  permuted = data[sample == sample[1], {
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
      if ( !allow.identical && identical(tmp, outcome) ) success=F
    }
    .(trial = trial, outcome = tmp) 
  }, by=subject]
  duplIDs = permuted[, .(N=.N/uniqueN(trial)), by=subject] #check for duplicate trial IDs
  if (any(duplIDs$N > 1)) { #duplicates
    if (duplIDs[, uniqueN(N)] > 1) { #different number of duplicates across subjects
      if ( is.float(duplIDs$N) ) { #no consistency, go through every trial
        setkeyv(permuted, c("subject","trial"))  
        data[, outcome := rep_len(permuted[.BY, outcome], .N), by=.(subject,trial)] #slow
      } else { #every trial ID repeated the same number of times
        ns = data[, uniqueN(sample)]
        data[, outcome := permuted[, {
          x = duplIDs[.BY,N]
          start = seq(1, .N, by=x)
          end = start+x-1
          .( y = outcome[ apply( mapply(seq, start, end), 2, rep, times=ns) ] )
        } , by=subject]$y ]
      }
    } else { #number of duplicates is constant across subjects
      n = duplIDs$N[1]
      start = seq(1, nrow(permuted), by=n)
      end = start+n-1
      data[, outcome := permuted$outcome[ apply( mapply(seq, start, end), 2, rep, times=uniqueN(sample)) ] ]
    }
  } else { #all unique trial IDs
    data[, outcome := rep(permuted$outcome, each=uniqueN(sample))] #much faster
  }
  return( setkeyv(data, keys[ !keys %in% "outcome" ])[] ) #don't re-order for outcome
}


#### private

.FIT <- function(train, train.outcome, args) {
  ## model fit with less overhead for repeated execution
  train = scale(train, args$center, args$scale)
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
  setattr(fit, "center", center)
  return( setattr(fit, "scale", scale) )
}

.PREDICT <- function(fit, test, test.outcome, args) {
  ## model prediction where test is a matrix
  test = scale(test, attr(fit, "center"), attr(fit, "scale"))
  preds = predict(fit, test, verbose=F)[[1]]
  if ( args$model != "REG" ) {
    preds = ordered(preds, levels(test.outcome))
    if (args$multiClass) {
      auc = pROC::multiclass.roc(test.outcome, preds, direction="<", algorithm=3)$auc[1]
    } else {
      auc = pROC::auc(test.outcome, preds, direction="<", algorithm=3)[1]
    }
    CM = table(preds, test.outcome) #confusion matrix
    return( list(AUC = auc, correct = sum(diag(CM)), incorrect = sum(CM[upper.tri(CM)], CM[lower.tri(CM)])) )
  } else { #regression
    resids = test.outcome - preds
    R2 = sum(resids^2)/sum((test.outcome-mean(test.outcome))^2) #SSE/SST
    #1 - [(n-1)/(n-p)] * (SSE/SST) #adjustment doesn't work when n < p?
    return( list(RMSE = sqrt(mean(resids^2)), R2=1-R2, N=length(resids)) )
  }
}

.DECOMPOSE <- function(train, test, args) {
  ## data decomposition / feature transformation
  C = setkeyv( train[, {
    mat = as.matrix(.SD)
    C = cov( mat - matrix( colMeans(mat), nrow=.N, ncol=length(args$nm), byrow=T ) ) #with centering
    # C = cov(mat) #without centering
    lapply(1:length(args$nm), function(ch) C[,ch])
  }, .SDcols = args$nm, by=.(subject,trial,outcome)][, .idx := 1:length(args$nm)][, {
    lapply(.SD, mean) }, .SDcols=args$nmV, by=.(subject,outcome,.idx)], c("subject","outcome") )
#   C = setkeyv( train[, { 
#     # mat = as.matrix(.SD)
#     # C = cov( mat - matrix( colMeans(mat), nrow=.N, ncol=length(nm), byrow=T ) )
#     C = cov(.SD) #without centering
#     lapply(1:length(args$nm), function(ch) C[,ch])
#   }, .SDcols=args$nm, by=.(subject,outcome)], c("subject","outcome") )
  
  dat = C[, { #create spatial filters with training set covariance matrices
    C1 = as.matrix(.SD[ 1:(.N/2) ])
    C2 = as.matrix(.SD[ (.N/2+1):.N ])
    VD = eigen(C1+C2, symmetric=T)
    V = VD$vectors
    d = VD$values
    r = sum(d > 10^-6*d[1])
    if (r < 2*args$npattern) stop( "Too few data columns to calculate", 2*args$npattern, "filters." )
    P = diag(d[1:r]^-0.5) %*% t(V[,1:r])
    S1 = P %*% C1 %*% t(P)
    R = eigen(S1, symmetric=T)$vectors
    W = t(P) %*% R
    filters = W[,c(1:args$npattern, (ncol(W)-args$npattern+1):ncol(W))]
    #create features for ML by applying training set filters to both train and test set
    rbind(train,test)[subject == .BY$subject, { 
      features = as.matrix(.SD) %*% filters
      lapply(seq_len(2*args$npattern), function(f) {
        if (args$approx) {
          feat = var(features[,f])
          ifelse(args$logscale, ifelse(feat > 0, log(feat), -70), feat) #replace -Inf with -70 (~ log(4^-31))
        } else {
          features[,f]
        }
      })
    }, .SDcols=args$nm, by=.(subject,trial,outcome)][, subject:=NULL]
  }, .SDcols=args$nmV, by=subject]
  setkeyv(dat, c("subject","trial","outcome"))
  return(dat)
}


# .DECOMPOSE <- function(train, test, args) {
#   C = setkeyv( train[, {
#     mat = as.matrix(.SD)
#     C = cov( mat - matrix( colMeans(mat), nrow=.N, ncol=length(args$nm), byrow=T ) ) #with centering
#     # C = cov(mat) #without centering
#     lapply(1:length(args$nm), function(ch) C[,ch])
#   }, .SDcols = args$nm, by=.(subject,trial,outcome)][, .idx := 1:length(args$nm)][, {
#     lapply(.SD, mean) }, .SDcols=args$nmV, by=.(subject,outcome,.idx)], c("subject","outcome") )
#   
#   dat = C[, { #create spatial filters with training set covariance matrices
#     C1 = as.matrix(.SD[ 1:(.N/2) ])
#     C2 = as.matrix(.SD[ (.N/2+1):.N ])
#     VD = eigen(C1+C2, symmetric=T)
#     V = VD$vectors
#     d = VD$values
#     r = sum(d > 10^-6*d[1])
#     if (r < 2*args$npattern) stop( "Too few data columns to calculate", 2*args$npattern, "filters." )
#     P = diag(d[1:r]^-0.5) %*% t(V[,1:r])
#     S1 = P %*% C1 %*% t(P)
#     R = eigen(S1, symmetric=T)$vectors
#     W = t(P) %*% R
#     filters = W[,c(1:args$npattern, (ncol(W)-args$npattern+1):ncol(W))]
#     # filters = apply(filters, 2, function(w) w/sqrt(t(w)%*%((C1+C2)/2)%*%w)) #scale to unit variance
#     # filters = apply(filters, 2, function(w) w/sqrt(sum(w^2))) #normalize by 2-norm (euclidian distance for vectors)
#     #create features for ML by applying training set filters to both train and test set
#     rbind(train,test)[subject == .BY$subject, { 
#       features = as.matrix(.SD) %*% filters
#       lapply(seq_len(2*args$npattern), function(f) features[,f])
#     }, .SDcols=args$nm, by=.(subject,trial,outcome)][, subject:=NULL]
#   }, .SDcols=args$nmV, by=subject][, sample := c( train[, sample], test[, sample] )]
#   setkeyv(dat, c("subject","trial","sample","outcome"))
#   return(dat)
# }
