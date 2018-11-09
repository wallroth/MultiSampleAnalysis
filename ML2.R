## NEW FUNCTIONS
#decision to 'naive' probability: 1/(1+exp(d)) #uncalibrated -> see Platt Scaling

#### Convenience functions ----
pipeline <- function(data, slices, CV, GAT=F, grid=F, prob=F, scaling="CV", par.folds=F, ...) {
  ## pipeline of fitting and prediction, returning a data table output
  # GAT: generalization across time (cf. function GAT)
  # grid: grid search (cf. function grid_search) - ToDo
  # prob: return the probability/certainty of the decision (cf. LiblineaR/kernlab)
  # scaling: CV (on train set, applied to test set), full (slicewise on the full set), off
  # par.folds: parallelize across folds (T) or time (F)
  # ... : further arguments to tfit (kernel, type, C, balance, params)
  tstart = Sys.time()
  #args = do.call( .eval_ellipsis, modifyList(list(...), list(funStr=c("tfit","grid_search"))) )
  args = lapply(.eval_ellipsis("tfit", ...)$tfit[c("kernel","type","C","balance","params")], eval)
  scaling = tolower(scaling); z = scaling=="cv"
  PREP = .prepare_(data=data, slices=slices, CV=CV, GAT=GAT, z=scaling=="full", balance=args$balance, cut=T)
  SLICED = PREP$SLICED; y = PREP$y; itrain = PREP$itrain; itest = PREP$itest; w = PREP$w
  if (GAT) CVSLICED = PREP$CVSLICED
  lvls = sort(unique(y)); multi = length(lvls) > 2
  if (multi) lvls = apply(combn(lvls, 2), MARGIN=2, paste0, collapse="v") #for kernlab decision output
  else lvls = "" #only one value output from kernlab
  n = nrow(slices); m = ncol(slices)
  #internal output handler
  assemble <- function(y_pred, y_test, y_prob, y_dec) {
    out = list(ACC=as.integer(y_test==y_pred), target=y_test, predicted=y_pred)
    if (!is.null(y_dec)) {
      nm = colnames(y_dec)
      if (multi && is.null(nm)) colnames(y_dec) = lvls #kernlab output
      out$decision = y_dec
    }
    if (!is.null(y_prob)) {
      if (multi && !is.null(kernel)) colnames(y_prob) = lvls
      out$probability = y_prob
    }
    return(as.data.table(out))
  }
  
  if (grid) { #ToDo
    # params = args$grid_search$params; settings = lapply(args$grid_search$settings[-1], function(s) eval(s))
    # defaults = list(k=5, shuffle=F, stratified=T, within=T, optimize='AUC', crit=.95)
    # settings = c(settings, defaults[!names(defaults) %in% names(settings)])
    # ROI = settings$ROI
    # if (is.null(ROI)) ROI = 1:ncol(slices) #full time series
    return()
    #...
    
  } else { #no grid search
    kernel = args$kernel; type = args$type; C = args$C; params = args$params
    if (par.folds) { 
      #parallelize folds
      tab = foreach(k=1:length(itrain), .multicombine=T, .maxcombine=length(itrain)+1, .combine=function(...) 
        rbindlist(list(...), idcol="fold", use.names=T, fill=T)) %dopar%
        {
          scoreDT = rbindlist( lapply(SLICED, function(X) {
            f = fit_(X=X[ -itrain[[k]],], y=y[ -itrain[[k]] ], 
                     kernel=kernel, C=C, z=z, w=w, type=type, prob=prob, params=params)
            y_test = y[ itest[[k]] ]
            if (GAT) {
              y_test = rep(y_test, m)
              p = predict_(f, CVSLICED[[k]], prob=prob, z=z)
            } else {
              p = predict_(f, X[ itest[[k]],,drop=F], prob=prob, z=z)  
            }
            assemble(y_pred=p$predictions, y_test=y_test, y_prob=p$probabilities, y_dec=p$decisionValues)
          }), use.names=T, fill=T, idcol=if (GAT) "train" else "slice")
        }
      if (!"fold" %in% names(tab)) tab[, fold:=unique(CV$fold)]
    } else { #parallelize slices
      tab = foreach(X=SLICED, .multicombine=T, .maxcombine=length(SLICED)+1, .combine=function(...) 
        rbindlist(list(...), use.names=T, fill=T,idcol=if (GAT) "train" else "slice")) %dopar%
        {
          scoreDT = rbindlist( lapply(1:length(itest), function(k) {
            f = fit_(X=X[ -itrain[[k]],], y=y[ -itrain[[k]] ], 
                     kernel=kernel, C=C, z=z, w=w, type=type, prob=prob, params=params)
            y_test = y[ itest[[k]] ]
            if (GAT) {
              y_test = rep(y_test, m)
              p = predict_(f, CVSLICED[[k]], prob=prob, z=z)
            } else {
              p = predict_(f, X[ itest[[k]],,drop=F], prob=prob, z=z)  
            }
            assemble(y_pred=p$predictions, y_test=y_test, y_prob=p$probabilities, y_dec=p$decisionValues)
          }), idcol="fold", use.names=T, fill=T)
        }
    }
  } #end of non-grid decoding
  #add info
  trls = lapply(itest, function(i) PREP$trials[i])
  if (GAT) {
    tab[, test:=rep(1:m, each=length(trls[[fold]])), by=.(train,fold)]
    tab[, trial:=trls[[fold]], by=.(train,test,fold)]
    neword = c("train","test","fold","trial")
  } else {
    tab[, trial:=trls[[fold]], by=.(slice,fold)]
    if (n > 1) tab[, sample:=slices[,slice], by=.(slice,fold,trial)]
    neword = c("slice",if(n>1)"sample","fold","trial")
  }
  print( Sys.time()-tstart )
  return( setcolorder(tab, union(neword, names(tab)))[] )
}

pipe_permutations <- function(data, slices, permutations, GAT=F, scaling="full", summarize=T, par.perm=F, ...) {
  ## permutation pipeline, returning a data table output
  # permutations: output of permute_ or a single number specifying num permutations
  # GAT: generalization across time (cf. function GAT)
  # scaling: CV (on train set, applied to test set), full (slicewise on the full set), off
  # summarize: aggregate results across folds for a lean output
  # ... : further arguments to tfit (type, C, balance), create_CV, permute_, and score_
  tstart = Sys.time()
  args = do.call( .eval_ellipsis, modifyList(list(...), list(funStr=c("tfit","create_CV","permute_","score_"))) )
  CVset = lapply(args$create_CV[-1], function(s) eval(s))
  type = args$tfit$type; C = args$tfit$C; balance = eval(args$tfit$balance)
  kernel = args$tfit$kernel; params = eval(args$tfit$params)
  scaling = tolower(scaling); z = scaling=="cv"
  #get a CV scheme for preparation
  CV = create_CV(data, k=CVset$k, within=CVset$within, shuffle=CVset$shuffle, stratified=CVset$stratified)
  PREP = .prepare_(data=data, slices=slices, CV=CV, GAT=GAT, z=scaling=="full", balance=balance, cut=T)
  SLICED = PREP$SLICED; y = PREP$y; itrain = PREP$itrain; itest = PREP$itest; w = PREP$w
  if (GAT) CVSLICED = PREP$CVSLICED
  n = nrow(slices); m = ncol(slices)
  if (!is.matrix(permutations)) { #create permutation scheme
    permutations = permute_(y, n=permutations, min.tol=args$permute_$min.tol)
  }
  if (summarize) { #summarize output to save space
    lvls = sort(unique(y))
    multi = if (length(lvls) > 2) T else NULL #if multi-class compute average score
    binary = eval(args$score_$binary) #binary vs OvR - only relevant for multi-class
    if (par.perm) {
      if (GAT) byvar = c("train","test") else byvar = "slice" 
    } else {
      if (GAT) byvar = "test" else byvar = ""
    }
    summarizer <- function(result, lvls, multi, binary, byvar) {
      return( result[, {
        out = score_(target, predicted, lvls, hardvote=T, multi.avg=multi, binary=binary)
        if (is.logical(multi)) { correct = sum(ACC); cm = list(correct=correct, incorrect=.N-correct) } 
        else { cm = cm_(factor(predicted, levels=lvls), factor(target, levels=lvls)) }
        c(out, cm)
      }, by=byvar] )
    }
  }
  #internal CV handler
  doCV = CV[, uniqueN(fold)>1]
  CVdata = data[sample==1,key(data),with=F]
  perm_CV <- function(CVdata, y, doCV, n, CVset) {
    set(CVdata, j="outcome", value=y)
    CV = create_CV(CVdata, k=CVset$k, within=CVset$within, shuffle=CVset$shuffle, stratified=CVset$stratified)
    if (doCV) {
      itest = split(1:length(y), rep(CV$fold, each=n))
      return( list(fit=itest, predict=itest) )
    } else { #no CV
      return( list(fit=list(length(y)+1), predict=list(1:length(y))) )
    }
  }
  #internal output handler
  assemble <- function(y_pred, y_test) { 
    return( list(ACC=as.integer(y_test==y_pred), target=y_test, predicted=y_pred) )
  }
  if (par.perm) { #parallelize permutations
    .eval_package("iterators", load=T)
    tab = foreach(y=iter(permutations, by="column"), .multicombine=T, .maxcombine=ncol(permutations)+1, 
                  .combine=function(...) rbindlist(list(...), idcol="run")) %dopar%
                  {
                    icv = perm_CV(CVdata, y, doCV, n, CVset)
                    if (GAT) {
                      CVSLICED = lapply(icv$predict, function(i) { #ordered according to slices
                        do.call(rbind, lapply(SLICED, function(XX) XX[i,,drop=F]))
                      })
                    }
                    result = rbindlist( lapply(SLICED, function(X) {
                      rbindlist( lapply(1:length(icv$fit), function(k) {
                        f = fit_(X=X[ -icv$fit[[k]],], y=y[ -icv$fit[[k]] ], prob=F,
                                 kernel=kernel, C=C, z=z, w=w, type=type, params=params)
                        y_test = y[ icv$predict[[k]] ]
                        if (GAT) {
                          y_test = rep(y_test, m)
                          p = predict_(f, CVSLICED[[k]], prob=F, z=z)
                        } else {
                          p = predict_(f, X[ icv$predict[[k]],,drop=F], prob=F, z=z)  
                        }
                        assemble(y_pred=p$predictions, y_test=y_test) 
                      }), idcol="fold" )
                    }), idcol=if (GAT) "train" else "slice" )
                    if (GAT) result[, test:=rep(1:m, each=length(icv$predict[[fold]])), by=.(train,fold)]
                    if (summarize) summarizer(result, lvls, multi, binary, byvar)
                    else result
                  }
  } else { #parallelize slices
    ptab = list()
    for (r in 1:ncol(permutations)) {
      cat(ifelse(r%%50,".",r),sep="")
      y = permutations[,r] #replace y with current permutation
      ptab[[r]] = foreach(X=SLICED, .multicombine=T, .maxcombine=length(SLICED)+1, .combine=function(...) 
        rbindlist(list(...), idcol=if (GAT) "train" else "slice")) %dopar%
        {
          icv = perm_CV(CVdata, y, doCV, n, CVset)
          if (GAT) {
            CVSLICED = lapply(icv$predict, function(i) { #ordered according to slices
              do.call(rbind, lapply(SLICED, function(XX) XX[i,,drop=F]))
            })
          }
          result = rbindlist( lapply(1:length(icv$fit), function(k) {
            f = fit_(X=X[ -icv$fit[[k]],], y=y[ -icv$fit[[k]] ], prob=F,
                     kernel=kernel, C=C, z=z, w=w, type=type, params=params)
            y_test = y[ icv$predict[[k]] ]
            if (GAT) {
              y_test = rep(y_test, m)
              p = predict_(f, CVSLICED[[k]], prob=F, z=z)
            } else {
              p = predict_(f, X[ icv$predict[[k]],,drop=F], prob=F, z=z)  
            }
            assemble(y_pred=p$predictions, y_test=y_test)
          }), idcol="fold" )
          if (GAT) result[, test:=rep(1:m, each=length(icv$predict[[fold]])), by=fold]
          if (summarize) summarizer(result, lvls, multi, binary, byvar)
          else result
        }
    }
    tab = rbindlist(ptab, idcol="run")
    cat("\n")
  }
  if (GAT) setcolorder(tab, union(c("run","train","test"), names(tab)))
  print( Sys.time()-tstart )
  return(tab)
}

summarize <- function(tab, groupby=NULL, ...) {
  ## sumarize pipe output
  # tab: output from pipeline
  # groupby: character vector specifying grouping to aggregate by (first item is parallelized)
  # ...: further arguments to score_ (multi.avg, binary)
  tstart = Sys.time()
  args = lapply(.eval_ellipsis("score_", ...)$score_[c("multi.avg","binary")], eval)
  lvls = tab[, sort(unique(target))]
  multi = length(lvls) > 2
  if (multi) {
    if (is.null(args$multi.avg)) args$multi.avg = T
    else if (!args$multi.avg) {
      if (args$binary) lvlnm = apply(combn(lvls, 2), MARGIN=2, paste0, collapse="v")
      else lvlnm = sapply(lvls, function(l) paste0(l,"vR")) #OvR
    }
  }
  #probabilities present in tab?
  probnm = grep("prob", names(tab), value=T)
  prob = length(probnm) > 0
  #decision values present in tab?
  dvnm = grep("decision", names(tab), value=T)
  dv = length(dvnm) > 0
  #make iterator for parallelization
  if (is.null(groupby)) { groupby = ""; itab = list(tab) }
  else itab = split(tab, by=groupby[1]) 
  stab = foreach(tab=itab, .multicombine=T, .maxcombine=length(itab)+1, 
                 .combine=function(...) rbindlist(list(...))) %dopar%
                 {
                   tab[, {
                     out = score_(target=target, predicted=predicted, lvls=lvls, hardvote=T, 
                                  multi.avg=args$multi.avg, binary=args$binary)
                     if (multi) { 
                       correct = sum(ACC)
                       cm = list(correct=correct, incorrect=.N-correct)
                       if (!args$multi.avg) { #out is a matrix
                         out = c(as.list(out["AUCh",]), as.list(out["PRh",]))
                         names(out) = c(paste0("AUCh.",lvlnm), paste0("PRh.",lvlnm))
                       }
                     } else { cm = cm_(factor(predicted, levels=lvls), factor(target, levels=lvls)) }
                     out = c(out, cm)
                     if (prob) { out = c(out, lapply(.SD[, probnm, with=F], mean, na.rm=T)) }
                     if (dv) {
                       if (multi) {
                         dvals = as.matrix(.SD[,dvnm,with=F])
                       } else { #binary
                         if (length(dvnm) == 1) { #one column only (likely kernlab)
                           dvals = .SD[[dvnm]]
                         } else { #2 columns, probably LiblineaR - often inconsistent ordering
                           dvcol = apply(.SD[,dvnm,with=F], MARGIN=2, function(x) all(x==0))
                           if ( all(!dvcol) ) { #fix this by merging and switching signs to one class
                             tmp = .SD[,dvnm,with=F]
                             tmp[tmp[[1]]==0, (dvnm[1]) := -1*.SD, .SDcols=dvnm[2]]
                             dvals = tmp[[1]]
                           } else {
                             dvals = .SD[[ dvnm[!dvcol] ]]
                           }
                         }
                         pos.lab = predicted[dvals>0][1]
                         if (is.na(pos.lab)) pos.lab = lvls[-which(lvls == tab$predicted[1])]
                         if (!pos.lab == lvls[2]) dvals = -dvals #make sure the positive class is 2nd
                       }
                       soft = score_(target=target, predicted=dvals, lvls=lvls, hardvote=F, 
                                     multi.avg=args$multi.avg, binary=args$binary)
                       if (multi && !args$multi.avg) {
                         soft = c(as.list(soft["AUCs",]), as.list(soft["PRs",]))
                         names(soft) = c(paste0("AUCs.",lvlnm), paste0("PRs.",lvlnm))
                       }
                       out = c(out, soft)
                     }
                     out
                   }, by=groupby]
                 }
  print( Sys.time()-tstart )
  return(stab)
}


tfit <- function(data, slices, CV, kernel=NULL, type=NULL, C=1, balance=F, scaling="CV", prob=F, params=list()) {
  ## fit a model at each time step
  # slices: time steps as defined by slice.idx
  # CV: cross-validation to use as defined by create_CV
  # kernel: if NULL -> LiblineaR (fast), else kernlab's ksvm (cf. kernel parameter)
  # type: numeric, default 2 (cf. LiblineaR) or character, default 'C-svc' (cf. kernlab)
  # C: regularization constant, smaller = stronger regularization (= smaller weights)
  # balance: weight classes according to their distribution, if T weighted by ratio, can be named vector
  # scaling: CV (on train set, applied to test set), full (slicewise on the full set), off
  # prob: compute probability model if using kernlab (ignored if LiblineaR)
  # params: set bias for LiblineaR, otherwise cf. kpar @ kernlab
  # RETURNS a list with the structure: time points x folds | containing fits of class LiblineaR or ksvm
  tstart = Sys.time()
  scaling = tolower(scaling); z = scaling=="cv"
  PREP = .prepare_(data=data, slices=slices, CV=CV, GAT=F, z=scaling=="full", balance=balance, cut=T)
  SLICED = PREP$SLICED; y = PREP$y; itest = PREP$itrain; w = PREP$w
  if (ncol(slices)>1) {  ## DECODE - parallelized across time
    fits = foreach(X=SLICED, .combine=list, .multicombine=T, .maxcombine=ncol(slices)+1) %dopar%
      lapply(itest, function(test) fit_(X=X[-test,], y=y[-test], kernel=kernel, 
                                        C=C, z=z, w=w, type=type, prob=prob, params=params))
  } else { #save some computation time by parallelizing across folds
    X = SLICED[[1]]
    fits = foreach(test=itest, .combine=list, .multicombine=T, .maxcombine=length(itest)+1) %dopar%
      fit_(X=X[-test,], y=y[-test], kernel=kernel, C=C, z=z, w=w, type=type, prob=prob, params=params)
  }
  print( Sys.time()-tstart )
  return(fits)
}

grid_search <- function(data, slices, CV, scaling="CV", balance=F, prob=F, params=list(),
                        settings=list(k=5, shuffle=F, stratified=T, within=T, 
                                      optimize='AUC', crit=.95, ROI.max=NULL, ROI.min=NULL)) {
  ## fit a model at each time step with optimization of model and C parameter choice
  # slices: time steps as defined by slice.idx
  # CV: cross-validation to use as defined by create_CV
  # scaling: CV (on train set, applied to test set), full (slicewise on the full set), off
  # balance: weight classes during learning according to their distribution (i.e. give higher weight to rare classes)
  # prob: compute probability model if using kernlab (ignored if LiblineaR)
  # params: list of parameters to optimize over
  # - contains lists where each element is named for the kernel, "liblinear" is an option and assumed if no kernels are named
  # - e.g. list(liblinear=list(bias=c(0,1), type=1), vanilladot=list(C=c(.01,1)), rbfdot=list(sigma=c(1,"auto"), C=c(.01,1)))
  # settings: list of values which define the inner cross-validation and optimization criteria
  # - k, shuffle, stratified: inner folds params (cf. create_CV)
  # - optimize: score which should be optimized by the type/C param
  # - crit: performance criterion to calculate across time points -> .5 = median, 1 = maximum, 0 = minimum
  # - ROI.max: region of interest for optimization when you don't want to optimize over the full ts (e.g. exclude a baseline)
  # - ROI.min: region to minimize jointly with ROI.max, e.g. baseline period
  # RETURNS the optimized fits and a result data table showing the performance for each parameter combination
  # "fits" are a list with the structure: time points x folds
  tstart = Sys.time()
  if ( !is.list(params[[1]]) ) { #no nesting for kernels, assuming LiblineaR
    params = list(params)
    names(params) = "liblinear"
  }
  defaults = list(k=5, shuffle=F, stratified=T, within=T, optimize='AUC', crit=.95, ROI.max=NULL, ROI.min=NULL)
  settings = c(settings, defaults[!names(defaults) %in% names(settings)])
  midx = ifelse(settings$optimize == "PR", 2, 1) #1: AUC, 2: PR
  hard = !grepl("s", settings$optimize)
  ROI.max = settings$ROI.max; ROI.min = settings$ROI.min
  if (is.null(ROI.max)) ROI.max = 1:ncol(slices) #full time series
  # else if (is.null(ROI.min)) ROI.min = (1:ncol(slices))[-ROI.max]
  scaling = tolower(scaling); z = scaling=="cv"
  PREP = .prepare_(data=data, slices=slices, CV=CV, GAT=F, z=scaling=="full", balance=balance, cut=T)
  SLICED = PREP$SLICED; y = PREP$y; itest = PREP$itest; itrain = PREP$itrain; w = PREP$w
  lvls = sort(unique(y)) #for AUC scoring
  multavg = if (length(lvls) > 2) T else NULL
  #expand CV to num samples per slice and add sample info for inner CV creation
  CV = setkeyv(CV[rep(1:.N, each=nrow(slices))][, sample:=1:nrow(slices)], c("subject","trial","sample"))
  #if no outer CV duplicate for technical reasons to enable inner CV creation
  if ( CV[, uniqueN(fold)==1] ) CV = rbind(CV,copy(CV)[,fold:=fold+1L])
  #create inner loop of the CV
  itest = lapply(1:length(itest), function(k) {
    inCV = create_CV(CV[ -itest[[k]] ], k=settings$k, shuffle=settings$shuffle,
                     stratified=settings$stratified, within=settings$within)
    inCV = inCV[CV[ -itest[[k]], .(subject,trial,outcome)], on=.(subject,trial,outcome), nomatch=0] #resort
    list( train=itrain[[k]], inner=split(seq_along(y)[ -itrain[[k]] ], inCV$fold) )
  })
  #inner fit function which handles param input
  gridfit <- function(X,y,param,prob=F) {
    p = as.list(param)[!is.na(param)]
    kernel = p$kernel; C=p$C; type=p$type
    if (kernel == 'liblinear') kernel=NULL
    if (is.null(C)) C = 1
    if (!is.null(p$sigma) && grepl("auto",p$sigma)) p = list() #defaults to automatic in fit_
    return( fit_(X,y,kernel=kernel,C=C,z=z,w=w,type=type,prob=prob,
                 params=p[!names(p) %in% c("kernel","C","type")]) )
  }
  ## DECODE | parallelize folds so that we can optimize across time
  out = foreach(k=itest, .combine=list, .multicombine=T, .maxcombine=length(itest)+1) %dopar%
  {
    y_test = y[unlist(k$inner)] #order of predictions needs to be matched
    result = rbindlist(lapply(params, expand.grid), idcol="kernel", fill=T)[, kernel:=tolower(kernel)]
    for (p in 1:nrow(result)) { #iterate the parameter space for the region of interest
      #maximize
      ymax = lapply(SLICED[ROI.max], function(X) { #iterate time points
        lapply(k$inner, function(kk) { #iterate inner folds
          train = c(k$train,kk)
          fit = gridfit(X=X[-train,], y=y[-train],param=result[p]) 
          if (hard) predict_(fit, X[kk,,drop=F], prob=F, z=z)$predictions
          else predict_(fit, X[kk,,drop=F], prob=F, z=z)$decisionValues
        })
      })
      if (hard) ymax = lapply(ymax, unlist)
      else ymax = lapply(ymax, do.call, what=rbind)
      scores_max = sapply(ymax, function(p) {
        score_(target=y_test, predicted=p, lvls=lvls, hardvote=hard, multi.avg=multavg)[[midx]] })
      #minimize if ROI.min is non-NULL
      if (!is.null(ROI.min)) { 
        ymin = lapply(SLICED[ROI.min], function(X) { #iterate time points
          lapply(k$inner, function(kk) { #iterate inner folds
            train = c(k$train,kk)
            fit = gridfit(X=X[-train,], y=y[-train],param=result[p]) 
            if (hard) predict_(fit, X[kk,,drop=F], prob=F, z=z)$predictions
            else predict_(fit, X[kk,,drop=F], prob=F, z=z)$decisionValues
          })
        })
        if (hard) ymin = lapply(ymin, unlist)
        else ymin = lapply(ymin, do.call, what=rbind)
        scores_min = sapply(ymin, function(p) {
          score_(target=y_test, predicted=p, lvls=lvls, hardvote=hard, multi.avg=multavg)[[midx]] })
      }
      #criterion calculated across time and added to result
      result[p, max:=quantile(scores_max, probs=settings$crit)]
      if (!is.null(ROI.min)) result[p, ':=' (min=quantile(scores_min, probs=settings$crit),minvar=sd(scores_min))]
    }
    #score the result and create fits with best params
    if (!is.null(ROI.min)) result[, score:=max-abs(.5-min)] else result[, score:=max]
    result[, best:=""][which.max(score), best:="*"]
    p = rbindlist(lapply(params, expand.grid), idcol="kernel", fill=T)[, kernel:=tolower(kernel)]
    fits = lapply(SLICED, function(X) 
      gridfit(X=X[-k$train,], y=y[-k$train], param=p[result[which.max(score), which=T]], prob=prob))
    list(fits=fits, result=result) #return
  } #end of foreach
  print( Sys.time()-tstart )
  #if there was no CV foreach removes the outer nest of the list
  if ( all(c("fits", "result") %in% names(out)) ) return(out)
  #else... format output to time x fold
  result = rbindlist(lapply(out, "[[", "result"), idcol="fold")
  fits = lapply(1:length(SLICED), function(i) lapply(out, function(k) k$fits[[i]]))
  return(list(fits=fits, result=result))
}


tpredict <- function(data, slices, CV, fits, scaling="CV", prob=F, coupler="minpair", par.folds=T) {
  ## at each time step predict the outcome given a model
  # slices: time steps as defined by slice.idx (used during fitting)
  # CV: cross-validation to use as defined by create_CV (used during fitting)
  # fits: model fits to use as defined by tfit/grid_search
  # scaling: CV (on train set, applied to test set), full (slicewise on the full set), off
  # prob: return the probability/certainty of the decision (for LogReg if LiblineaR or if prob = T @ tfit with kernlab)
  # coupler: cf. couple @ kernlab - used for class probabilitiy estimates for > 2 classes
  # par.folds: parallelize across folds (T) or time (F) - unless very few folds it is usually faster with T
  # RETURNS predictions, a list of time points with the predictions per fold
  tstart = Sys.time()
  scaling = tolower(scaling); z = scaling=="cv"
  PREP = .prepare_(data=data, slices=slices, CV=CV, GAT=F, z=scaling=="full")
  SLICED = PREP$SLICED; itest = PREP$itest
  if ( class(fits[[1]]) != 'list' & ncol(slices)>1 ) fits = lapply(fits, list) #no CV nesting of fits
  ## PREDICT
  if (par.folds) {  #parallelize folds 
    if (ncol(slices)>1) {
        CVfits = lapply(1:length(itest), function(k) { #resort the iterator for folds
          lapply(1:length(SLICED), function(i) list(fit=fits[[i]][[k]], test=itest[[k]]))
        })
      out = foreach(fold=CVfits, .combine=list, .multicombine=T, .maxcombine=length(CVfits)+1) %dopar%
        lapply(1:length(fold), function(i) 
          predict_(fit=fold[[i]]$fit, X_test=SLICED[[i]][fold[[i]]$test,,drop=F], prob=prob, z=z, coupler=coupler))
      #resort out to match fit structure of time x fold instead of fold x time
      if ( length(itest) > 1 ) out = lapply(1:length(SLICED), function(i) lapply(out, "[[", i))
    } else { #no slicing
      CVfits = lapply(1:length(itest), function(k) { list(fit=fits[[k]], test=itest[[k]]) })
      out = foreach(fold=CVfits, .combine=list, .multicombine=T, .maxcombine=length(CVfits)+1) %dopar%
        predict_(fit=fold$fit, X_test=SLICED[[1]][fold$test,,drop=F], prob=prob, z=z, coupler=coupler)
    }
  } else { #parallelize time
    CVfits = lapply(1:length(SLICED), function(i) { #add time info to the iterator
      lapply(1:length(itest), function(k) list(fit=fits[[i]][[k]], test=itest[[k]], slice=i) )
    })
    out = foreach(cvf=CVfits, .combine=list, .multicombine=T, .maxcombine=length(CVfits)+1) %dopar%
      lapply(cvf, function(k) predict_(fit=k$fit, X_test=SLICED[[k$slice]][k$test,,drop=F], 
                                       prob=prob, z=z, coupler=coupler))
  }
  print( Sys.time()-tstart )
  return(out)
}

GAT <- function(data, slices, CV, fits, scaling="CV", prob=F, coupler="minpair") {
  ## Generalization across time: given a model, predict the outcome at all time steps
  # slices: time steps as defined by slice.idx (used during fitting)
  # CV: cross-validation to use as defined by create_CV (used during fitting)
  # fits: model fits to use as defined by tfit/grid_search
  # scaling: CV (on train set, applied to test set), full (slicewise on the full set), off
  # prob: return the probability/certainty of the decision (for LogReg if LiblineaR or if prob = T @ tfit with kernlab)
  # coupler: cf. couple @ kernlab - used for class probabilitiy estimates for > 2 classes
  # RETURNS predictions, a list with the structure: training time x folds x generalizations
  # if a fold contains multiple trials the rows in the generalization matrix are trial x timepoint
  # e.g. with 2 trials and 2 timepoitns: trial 1, time 1, trial 1, time 2, trial 2, time 1, trial 2, time 2
  tstart = Sys.time()
  scaling = tolower(scaling); z = scaling=="cv"
  PREP = .prepare_(data=data, slices=slices, CV=CV, GAT=T, z=scaling=="full")
  CVSLICED = PREP$CVSLICED
  if ( class(fits[[1]]) != 'list' ) fits = lapply(fits, list) #no CV nesting of fits
  ## PREDICT - outer nest is training time
  out = foreach(f=fits, .combine=list, .multicombine=T, .maxcombine=length(fits)+1) %dopar%
  {
    lapply(1:length(f), function(k) { #iterate folds
      predict_(fit=f[[k]], X_test=CVSLICED[[k]], prob=prob, z=z, coupler=coupler) #generalizations
    })
  }
  print( Sys.time()-tstart )
  return(out) # out has the structure: training time x fold x generalizations
}


formatDT <- function(data, slices, CV, predictions, GAT=F) {
  ## make a formatted data table with the predictions from tpredict/GAT
  # CV: cross-validation as defined by create_CV (used for predictions)
  # predictions: list of predictions as defined by tpredict/GAT
  # GAT: whether predictions were generated by tpredict or GAT
  # RETURNS a data table with 1 performance instance per row
  tstart = Sys.time()
  PREP = .prepare_(data=data, slices=slices, CV=CV)
  y = PREP$y; itest = PREP$itest #y sorted as in CV
  lvls = sort(unique(y)); multi = length(lvls) > 2
  if (multi) lvls = apply(combn(lvls, 2), MARGIN=2, paste0, collapse="|") #for kernlab decision output
  else lvls = "" #only one value output from kernlab
  n = nrow(slices); m = ncol(slices) #for test info
  #internal output handlers
  assemble <- function(y_pred, y_test, y_prob, y_dec) {
    out = list(ACC=as.integer(y_test==y_pred), target=y_test, predicted=y_pred)
    kernel = NULL
    if (!is.null(y_dec)) {
      nm = colnames(y_dec)
      if (is.null(nm)) {
        colnames(y_dec) = lvls #kernlab output
        kernel = T
      }
      out$decision = y_dec
    }
    if (!is.null(y_prob)) {
      if (multi && !is.null(kernel)) colnames(y_prob) = lvls
      out$probability = y_prob
    }
    return(as.data.table(out))
  }
  make_output <- function(folds) {
    if ( "predictions" %in% names(folds) ) folds = list(folds) #insert nesting (dropped by foreach)
    scoreDT = rbindlist( lapply(1:length(folds), function(k) {
      p = folds[[k]]
      y_test = y[ itest[[k]] ]
      if (GAT) y_test = rep(y_test, m)
      assemble(y_pred=p$predictions, y_test=y_test, y_prob=p$probabilities, y_dec=p$decisionValues)
    }), idcol="fold", use.names=T, fill=T )
    return(scoreDT)
  }
  #iterate training times
  tab = foreach(folds=predictions, .combine=function(...) 
    rbindlist(list(...), use.names=T, fill=T, idcol=if (GAT) "train" else "slice"), 
    .multicombine=T, .maxcombine=length(predictions)+1) %dopar% make_output(folds)
  
  if (ncol(slices)==1) tab[, slice:=1L]
  #add info
  trls = lapply(itest, function(i) PREP$trials[i])
  if (GAT) {
    tab[, test:=rep(1:m, each=length(trls[[fold]])), by=.(train,fold)]
    tab[, trial:=trls[[fold]], by=.(train,test,fold)]
    neword = c("train","test","fold","trial")
  } else {
    tab[, trial:=trls[[fold]], by=.(slice,fold)]
    if (n > 1) tab[, sample:=slices[,slice], by=.(slice,fold,trial)]
    neword = c("slice",if(n>1)"sample","fold","trial")
  }
  print( Sys.time()-tstart )
  return( setcolorder(tab, union(neword, names(tab)))[] )
}


tadjust <- function(data, tmin=2, tdiff=1, tvar="slice", group="") {
  ## temporal adjustment for multiple NHT
  # tmin: minimum num of consecutively significant timepoints (p<.05)
  # tdiff: time steps between samples
  # tvar: variable in data which codes the time info
  # group: variable(s) by which the adjustment should be stratified
  setnames(data, tvar, "slice")
  out = data[p<.05, {
    if (length(slice) >= tmin) {
      consecs = zoo::rollapply(slice, width=tmin, FUN=function(x) if (all(diff(x)==tdiff)) x else NULL)
      .SD[slice%in%as.vector(consecs)]
    } else {
      .SD[slice%in%NULL]
    }
  }, by=group]
  setnames(data, "slice", tvar) #change back
  return( setnames(out, "slice", tvar)[] )
}


create_CV <- function(data, k=NULL, within=T, shuffle=F, stratified=T) {
  ## create a CV scheme (all samples of a trial are assigned to the same fold)
  # k: number of folds, if between 0 and 1 percentage of data in the test set; if NULL: LO(S)O
  # within: make a within/across subjects scheme or a between subjects scheme
  # shuffle: randomly assign folds or choose in order of appearance (chronological)
  # stratified: balance outcome across folds, else random assignment (ignored if within=F)
  # - for LOO with stratification you get as many folds as the least common class
  # - without stratification LOO is simply one trial per fold
  # RETURNS a CV scheme as a data table with columns subject|fold|trial|outcome
  data = .check_(data)[sample==1, key(data), with=F] #only keep the relevant stuff for this
  if (within) { #within/across
    folds = data[, {
      class_tab = table(outcome) #unique classes and counts
      class_tab = class_tab[class_tab > 0] #in case of undropped factor levels
      if (is.null(k)) { #LOO if k not specified
        if (stratified) k = min(class_tab) else k = .N
      } 
      if (k < 1) { #percentage of trials per outcome which go into test set
        if (stratified) {
          ntrials = round(class_tab*k)
          trial_by_outcome = split(trial, outcome)
          test = as.vector( sapply(names(trial_by_outcome), function(y) {
            if (shuffle) sample(trial_by_outcome[[y]], size=ntrials[y])
            else trial_by_outcome[[y]][1:ntrials[y]] #chronological order
          }) )
        } else { #not stratified
          if (shuffle) test = sample(trial, size=round(.N*k))
          else test = trial[1:round(.N*k)]
        }
        .(fold = 1, trial = test, outcome = outcome[test]) #return
      } else { #fold split
        if (stratified) {
          if ( any(class_tab%/%k < 1) ) stop( "Cannot have more folds than trials per class." )
          folds = rep(0, .N) #fold vector
          for (i in 1:length(class_tab)) {
            nreps = class_tab[i]%/%k #num times of class per fold
            rest = class_tab[i]%%k
            if (shuffle) {
              folds[ outcome==names(class_tab)[i] ] = c( sample(rep(1:k, nreps)), sample(1:k, rest) ) #append if any are left
            } else { #chronological order
              nreps = rep(nreps, k)
              if (rest > 0) nreps[1:rest] = nreps[1:rest] + 1 #left-overs need to be allocated chronologically
              folds[ outcome==names(class_tab)[i] ] = rep(1:k, nreps)
            }
          }
        } else { #non-stratified
          folds = rep(1:k, each=.N/k)
          rest = .N-length(folds)
          if (rest > 0) folds = sort(c(folds, 1:rest))
          if (shuffle) folds = sample(folds)
        }
        .(fold = folds, trial = trial, outcome = outcome) #return
      }
    }, by=subject]
  } else { #between
    #create folds with all trials from n subjects
    folds = setkeyv( data[, .(subject = unique(subject))][, {
      if (is.null(k)) k = .N #LOSO
      if (k < 1) { #percentage of subjects which go into test set
        if (shuffle) subs = sample( subject, size=round(.N*k) )
        else subs = subject[1:round(.N*k)]
        .(fold = 1, subject = subs) #return
      } else { #fold split
        if (k > .N) stop( "Cannot have more folds than subjects." )
        folds = rep(1:k, each=.N/k)
        rest = .N-length(folds)
        if (rest > 0) folds = sort(c(folds, 1:rest))
        if (shuffle) folds = sample(folds)
        .(fold = folds, subject = subject)
      }
    }], c("fold", "subject"))[data[, .(subject,trial,outcome)], on=.(subject)]
  }
  return( setkeyv(folds, c("subject","fold","trial"))[] )
}

sliceidx <- function(n, window=1) {
  ## create indices for slicing of time series data
  #n: number of samples per trial or vector specifying min/max sample number or data
  #window: number of samples per slice that are contributed by each trial
  #RETURNS: a matrix where column = slice number, rows = slice samples
  if (is.data.frame(n)) n = .check_(n)[, range(sample)]
  min.n = ifelse(length(n) > 1, min(n), 1)
  max.n = max(n)
  if ( window > max.n-min.n+1 ) stop( "Window must be smaller or equal the number of slices." )
  start = seq(min.n, 1+max.n-window, by=window)
  end = seq(min.n+window-1, max.n, by=window)
  idx = mapply(seq, start, end) 
  #ensure idx is a matrix and not a vector when window = 1
  if ( !is.matrix(idx) ) idx = matrix(idx, nrow=window, byrow=T)
  return(idx) 
}


normalize <- function(data, baseline, average=F, scale=F) {
  ## baseline normalization
  # baseline: vector with sample numbers that represent baseline or single number setting endpoint of baseline
  # average: subtract the average across trials (T) or subtract for each trial individually (F)
  # scale: perform scaling in addition to centering
  data =.check_(data)
  if ( data[, uniqueN(sample)] <= 1 ) stop( "Normalization requires at least 2 samples per trial. ")
  if (length(baseline) == 1) baseline = 1:baseline
  nm = setdiff( names(data), key(data) )
  if (average) byvar = "subject" else byvar = c("subject","trial")
  m = data[sample %fin% baseline, lapply(.SD, mean), .SDcols = nm, by=byvar]
  if (scale) s = data[sample %fin% baseline, lapply(.SD, sd), .SDcols = nm, by=byvar]
  #expand m/s for mathematical operation between equally sized data.tables
  N = data[, .N, by=byvar][, N]
  m = m[rep(1:.N,times=N)]
  X = data[, nm, with=F] - m[, nm, with=F]
  if (scale) {
    s = s[rep(1:.N,times=N)]
    X = X / s[, nm, with=F]
  }
  return( .check_(cbind( data[, !nm, with=F], X)) )
}



#### Core functions ----
fit_ <- function(X, y, kernel=NULL, C=1, z=T, w=NULL, type=NULL, prob=T, params=list()) {
  # type: numeric (cf. LiblineaR) or character (cf. kernlab)
  # params: bias (LiblineaR) or cf. kpar of kernlab
  if (is.null(type)) type = ifelse(is.null(kernel), 2, "C-svc")
  if (is.logical(w) && w) w = .balance_(y)
  if (is.null(kernel)) { #LiblineaR
    if (is.null(params$bias)) params$bias = 1
    if (z) X = scale(X)
    fit = LiblineaR::LiblineaR(data=X, target=y, type=type, cost=C, wi=w, bias=params$bias)
    if (z) fit$center = attr(X, "scaled:center")
    if (z) fit$scale = attr(X, "scaled:scale")
  } else { #kernlab
    if (length(params) == 0 && kernel %in% c("rbfdot","laplacedot")) params = "automatic"
    fit = kernlab::ksvm(X, y, kernel=kernel, type=type, prob.model=prob, C=C,
                        fit=F, scaled=z, kpar=params, class.weights=w)
  }
  return(fit)
}


predict_ <- function(fit, X_test, prob=F, z=T, coupler="minpair") {
  if (class(fit) == "LiblineaR") {
    if (z) X_test = scale(X_test, fit$center, fit$scale)
    p = predict(fit, X_test, proba=prob, decisionValues=T)
    i = 0
    while ( any(is.na(p$decisionValues)) && i <= 100 ) { #weird inconsistent NaN behavior
      p = predict(fit, X_test, proba=prob, decisionValues=T)
      i = i+1
    }
    if (i > 100) p$predictions[ is.na(p$decisionValues[,1]) ] = NaN #couldn't fix
  } else { #kernlab presumed
    p = list(predictions = kernlab::predict(fit, X_test, type="response"), 
             decisionValues = kernlab::predict(fit, X_test, type="decision"))
    if (prob) p$probabilities = kernlab::predict(fit, X_test, type="probabilities", coupler=coupler)
  }
  return(p)
}


score_ <- function(target, predicted, lvls, hardvote=T, multi.avg=NULL, binary=T) {
  #hardvote: predicted is label outputs (doesn't matter if LiblineaR or kernlab)
  #multi.avg: if non-NULL computes multi-class AUC (return mean across comparisons if T)
  # -- if binary: compute pair-wise (if hardvote = F assumes kernlab pairwise decision scores)
  # -- if not binary: compute OvR (if hardvote = F assumes LiblineaR OvR decision scores)
  if (is.null(multi.avg)) { #binary 
    if (hardvote) {
      scores = list(AUCh = auc_(target, predicted, levels=lvls),
                    PRh = aucsoft_(target, predicted, pos.label=lvls[2], ROC=F))
    } else { #soft vote (predicted = decision scores)
      scores = list(AUCs = aucsoft_(target, predicted, pos.label=lvls[2], ROC=T),
                    PRs = aucsoft_(target, predicted, pos.label=lvls[2], ROC=F))
    }
  } else { #multiclass
    if (binary) { #iterate binary pairs (native kernlab solution)
      pairs = combn(lvls, m=2)
      if (hardvote) { #doesn't matter if labels come from LiblineaR or kernlab
        scores = apply(pairs, MARGIN=2, function(x) {
          idx = target %in% x
          r = as.integer(target[idx] == x[2]) #change to 0/1
          p = as.integer(predicted[idx] == x[2]) #change to 0/1
          c(AUCh = auc_(r, p, levels=c(0,1)), PRh = aucsoft_(r, p, pos.label=1, ROC=F))
        })
      } else { #softvote - assuming kernlab output with pairwise decision scores
        scores = sapply(1:ncol(pairs), function(i) {
          idx = target %in% pairs[,i]
          r = as.integer(target[idx] == pairs[2,i]) #change to 0/1
          p = predicted[idx,i]
          c(AUCs = aucsoft_(r, p, pos.label=1, ROC=T), PRs = aucsoft_(r, p, pos.label=1, ROC=F))
        })
      }
    } else { #OvR (native LiblineaR solution)
      if (hardvote) {
        scores = sapply(lvls, function(x) {
          r = as.integer(target == x)
          p = as.integer(predicted == x)
          c(AUCh = auc_(r, p, levels=c(0,1)), PRh = aucsoft_(r, p, pos.label=1, ROC=F))
        })
      } else { #softvote ovR - assuming LiblineaR output with OvR decision score columns
        scores = sapply(seq_along(lvls), function(i) {
          r = as.integer(target == lvls[i])
          p = predicted[,i] #assumes matching order between lvls and columns
          c(AUCs = aucsoft_(r, p, pos.label=1, ROC=T), PRs = aucsoft_(r, p, pos.label=1, ROC=F))
        })
      }
    }
    if (multi.avg) scores = as.list(rowMeans(scores)) #average over binary comparisons
  }
  return(scores)
}

auc_ <- function(target, predicted, levels) {
  #credits to https://github.com/benhamner/Metrics/blob/master/R/R/metrics.r
  r <- frank(predicted)
  yidx = target == levels[2]
  n_pos <- sum(yidx)
  n_neg <- length(target) - n_pos
  auc <- (sum(r[yidx]) - n_pos*(n_pos+1)/2) / (n_pos*n_neg)
  return(auc)
}

aucsoft_ <- function(target, decisions, pos.label, ROC=T) {
  posidx = target == pos.label
  pos_scores = decisions[posidx]
  neg_scores = decisions[!posidx]
  if (ROC) auc = PRROC::roc.curve(pos_scores, neg_scores)$auc #ROC curve
  else auc = PRROC::pr.curve(pos_scores, neg_scores)$auc.integral #precision-recall curve
  return(auc)
}

cm_ <- function(y_pred, y_test) { #confusion matrix
  cm = table(y_pred, y_test, exclude=NULL)
  tp = cm[2,2] #true positives / hits
  tn = cm[1,1] #true negatives / correct rejection
  fp = cm[2,1] #false positives / false alarms
  fn = cm[1,2] #false negatives / misses
  return(list(tp=tp, tn=tn, fp=fp, fn=fn)) #ACC = (tn+tp)/(tp+tn+fp+fn)
}

permute_ <- function(y, n=1, min.tol=0, max.tol=1-min.tol) {
  ## generate n unique permutations of the target vector y
  # min/max.tol: min/max percentage of y instances that have to be changed
  tab = table(y)
  perms = replicate(n, sample(y))
  lower = sum(tab*min.tol) #min num that needs to be switched
  upper = sum(tab*max.tol) #max num that needs to be switched
  crit = c(lower, upper)
  switched = perms != y
  check = colSums(switched) %between% crit & (!duplicated(cbind(y,perms), MARGIN=2))[-1]
  while (!all(check)) { #fix the insufficient/duplicate permutations
    pidx = which(!check)
    perms[,pidx] = replicate(length(pidx), sample(y))
    switched[,pidx] = perms[,pidx] != y
    check = colSums(switched) %between% crit & (!duplicated(cbind(y,perms), MARGIN=2))[-1]
  }
  return(perms)
}


#### Private functions ----
.check_ <- function(data) {
  cn = c("subject", "trial", "sample", "outcome")
  if ( is.data.table(data) && all(cn %in% key(data)) ) return(data) #everything OK
  #else...
  if ( !all(cn %in% names(data)) ) stop( "Info columns missing. Please refer to data.setinfo." )
  if ( !is.data.table(data) ) data = as.data.table(data) #if not a data table
  if ( !all(cn %in% key(data)) ) setkeyv(data, cn) #if keys are not set
  if ( !identical(names(data)[1:length(cn)], cn) ) setcolorder( data, union(cn, names(data)) ) #order columns
  return(data) #data changed in place if already a DT
}

.convert_ <- function(y) {
  if (is.factor(y) | is.character(y)) {
    y = as.character(y)
    if ( !all(grepl("[0-9]", y)) ) { #non-numeric values
      y = frank(y, ties.method="dense")-1L #assign numeric values 0:c
    }
  }
  return( as.numeric(y) )
}

.slice_ <- function(data, slices) {
  nm = setdiff( names(data), key(data) )
  X = unname( as.matrix(data[, .SD, .SDcols=nm]) ) #transform to matrix (measurements only)
  SLICED = lapply(1:ncol(slices), function(i) { #split into slices
    X[ data[ sample %fin% slices[,i], which=T ], ]
  })
  return(SLICED)
}

.prepare_ <- function(data, slices, CV, GAT=F, z=F, balance=F, cut=T) {
  data = .check_(data)
  if (cut) SLICED = .slice_(data, slices) else SLICED=NULL #split into slices
  if (z && cut) SLICED = lapply(SLICED, scale) #scale once
  S1 = data[ sample %fin% slices[,1], key(data), with=F ] #first slice
  CV = CV[S1, on=.(subject,trial,outcome)] #resort according to data and expand CV to n samples
  y = .convert_(S1$outcome) #S1 contains outcome according to n samples
  if ( CV[, uniqueN(fold)>1] ) {
    itest = split(1:length(y), CV$fold)
    itrain = itest
  } else { #no CV
    itest = list(1:length(y)) #include all samples
    itrain = list(length(y)+1) #exclusion of out-of-range index does nothing
  }
  if (is.logical(balance) && balance) { 
    w = T
  } else if (is.numeric(balance)) { #numeric
    names(balance) = sort(unique(y))
    w = balance
  } else { #FALSE
    w = NULL
  }
  if (!is.null(w)) y = factor(y) #otherwise kernlab fails
  PREP = list(SLICED=SLICED, y=y, itest=itest, itrain=itrain, w=w, trials=CV$trial)
  if (GAT) {
    PREP$CVSLICED = lapply(itest, function(i) { #ordered according to slices
      do.call(rbind, lapply(SLICED, function(X) X[i,,drop=F])) 
    })
  }
  return(PREP)
}

.balance_ <- function(y) {
  ratio = 1-table(y)/length(y)
  return(ratio/max(ratio)) #w
}


.eval_package(c("data.table","foreach","fastmatch"), load=T)