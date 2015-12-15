# utilities: often needed helper functions

is.datacol <- function(data, tol = .Machine$double.eps^0.5) {
  ## check if columns of data are measurements (floating numbers)
  ## used to exclude info columns in a df (e.g. sample number, outcome)
  sapply( 1:ncol(data), function(i) {
    if ( is.numeric( data[,i] ) ) {
      #are there non-integer values in the column?
      return( sum( abs( data[,i] - round(data[,i]) ) < tol ) < nrow(data) ) 
    }
    return( FALSE ) #non-numeric
  })
}

.is.sliced <- function(data) {
  ## check if data is a slice list
  if ( class(data) == "list" & length(data) > 1 ) {
    nsamples = max(data[[1]][,1])-min(data[[1]][,1])+1
    return( nrow(data[[1]]) != nsamples ) #TRUE if sliced
  }
  return(FALSE) #not a list
}

data.append_info <- function(data, trials=NULL, outcome=NULL) {
  ## function to append required info columns (sample numbers and/or outcome) to raw data
  ## which only contains measurements (or at least not the required info cols)
  #INPUT ---
  #data: df or list of trials
  #trials: int, if specified sample numbers will be added to 1st col accordingly
  #outcome: vector with 1 value per trial, if specified will be added to 2nd col
  #NOTES ---
  #if the sample number is in the 1st col and only outcome should be added, trials should
  #not be specified, otherwise samples will be added again.
  
  data = data.check(data, aslist=F)
  if ( !is.null(outcome) ) {
    if ( !is.null(trials) ) { #trials specified
      if ( (nrow(data)/trials) %% 1 > 0 ) {
        stop( "Either your specified trial number is incorrect or your trials differ in length." )
      }
      if ( length(outcome) != trials ) {
        stop( "Length of the outcome must match number of trials." )
      }
      outcome = rep( outcome, each=nrow(data)/trials )
      if ( length(outcome) < nrow(data) ) {
        stop( "Number of trials does not match the length of the data." )
      }
      data = cbind(outcome, data)
    } else if ( is.datacol(data)[1] ) { #sample number is not present in df and no trials specified
      stop( "No sample numbers found in 1st column and no trials specified.\n", 
            "Without proper trial info, outcome cannot be appended." )
    } else { #sample number is present in df
      nsamples = data.get_samplenum(data)
      outcome = rep( outcome, each=nsamples )
      if ( length(outcome) < nrow(data) ) {
        stop( "Length of the outcome must match number of trials." )
      }
      data = cbind(sample=data[,1], outcome, data[,-1]) #append to 2nd col
    }
  }
  if ( !is.null(trials) ) { #trials specified
    if ( (nrow(data)/trials) %% 1 > 0 ) {
      stop( "Either your specified trial number is incorrect or your trials differ in length." )
    }
    samples = rep( 1:(nrow(data)/trials), trials )
    if ( length(samples) < nrow(data) ) {
      stop( "Number of trials does not match the length of the data." )
    }
    data = cbind(sample=samples, data)
  }
  return( data.check(data, aslist=T, strip=F, transform=.transformed) )
}

data.get_samplenum <- function(data) {
  ## helper function to retrieve sample numbers per trial
  ## data is expected to have trials of equal length and continuous sample numbers
  if ( is.datacol(data)[1] ) {
    stop( "It seems your data does not have sample numbers in its first column.\n", 
          "Please provide that information. Use data.append_info at your convenience." )
  }
  nsamples = max(data[,1])-min(data[,1])+1
  if ( !isTRUE(all.equal( rep( min(data[,1]):max(data[,1]), nrow(data)/nsamples ), data[,1] )) ) {
    stop( "Perhaps your first column does not contain sample numbers.\n",
          "Other possible causes may be that your sample numbers are not evenly increasing,\n", 
          "or your trials are not of equal length. Please meet these assumptions." )
  }
  return( nsamples )
}

data.trials.split <- function(data, strip=T) {
  ## transform data into a list of trials according to sample numbering (1st col)
  ## returns a list of dfs with each element corresponding to a single trial
  #INPUT ---
  #strip: non-numeric columns are stripped (e.g. sample numbers, outcome)
  
  ## make sure data is a df 
  # can't call data.check because it might interfere with .transformed of previous call!
  if ( !is.data.frame(data) ) { #if not a df
    if ( class(data) == "list" ) { #if a list...
      #find out if trial list or slice list:
      if ( .is.sliced(data) ) {
        data = .data.unslice(data) 
      } else {
        return(data) #already trial list
      }
    } else {
      data = as.data.frame(data) #matrix/vector to df
    }
  }
  ##
  nsamples = data.get_samplenum(data) #continuous data requires 1st col with sample info
  if ( nsamples == nrow(data) ) { #no multi samples
    trials = data[,1]
  } else {
    trials = findInterval(1:nrow(data), seq(1, nrow(data), by=nsamples))  
  }
  if (strip) {
    return( split(data[, is.datacol(data)], trials) )
  }
  return( split(data, trials) )
}

data.check <- function(data, aslist=T, strip=F, transform=T) {
  ## check data format and return designated type
  #INPUT ---
  #data: df, matrix or list
  #aslist: True = return as list, False = df
  #strip: in case of aslist=T - see data.trials.split
  #transform: used to return to original data format using .transformed
  #RETURNS ---
  #data 
  #.transformed: bool in parent environment indicating if data was .transformed (True) or not (False)
  if (!transform) { #used if function is called with .transformed but no previous transformation happend
    suppressWarnings( rm(.transformed, pos=1) ) #remove from parent environment
    return(data)
  } 
  .transformed <<- FALSE
  if (aslist) {
    if ( class(data) != "list" ) { #if not a list
      .transformed <<- TRUE
      return( data.trials.split(data, strip=strip) ) #transform to list 
    }
    if ( strip ) { #even if a list strip non-data cols
      #FIX: possible misidentification in cases with only 1 sample per trial
      #obtain datacol info with higher sample number by rebinding temporarily
      datacols = is.datacol( as.data.frame( data.table::rbindlist(data) ) )
      return( lapply( data, function(d) d[, datacols] ) )
    }
    return(data) #already a list
  } 
  if ( !is.data.frame(data) ) { #if not a df
    if ( class(data) == "list" ) { #if a list...
      #find out if trial list or slice list:
      if ( .is.sliced(data) ) return( .data.unslice(data) )
      #retransform notice only if it was a trial list:
      .transformed <<- TRUE
      return( as.data.frame( data.table::rbindlist(data) ) ) #transform list to df
    }
  }
  return( as.data.frame(data) ) #df or matrix/vector to df
}

data.trials.demean <- function(data, end, start=1) {
  ## remove the mean of a range of data points from the other samples of a trial
  ## generally used to subtract the average of a baseline period, i.e. 1:baseline_end
  #INPUT ---
  #data: continuous df or list of trials, returns data in the same format
  #end: last sample of the baseline (range to subtract)
  #start: first sample of the baseline, defaults to first sample of the trial
  
  data = data.check(data, aslist=T, strip=F) #transform to list
  normdata = lapply(data, function(trial) {
    dcol = is.datacol(trial)
    baseline = scale(trial[start:end, dcol], scale=F)
    demeaned = scale(trial[ c( 0:(start-1), (end+1):nrow(trial) ), dcol], 
                     attr(baseline,"scaled:center"), scale=F)
    temp = cbind(trial[, !dcol], rbind( demeaned[0:(start-1),,drop=F], baseline, 
                                        demeaned[start:nrow(demeaned),,drop=F] ) ) 
    if ( sum(dcol) == 1 ) {
      #to prevent command as col name in case of single column data:
      names(temp)[is.datacol(temp)] = names(trial)[dcol]
    }
    temp #return to lapply
  })
  return( data.check(normdata, aslist=F, transform=.transformed) ) #transform back to df if needed
}

data.trials.remove_samples <- function(data, end, start=1, relabel=T) {
  ## remove datapoints from trials, specified by sample number (1st col)
  ## generally used to remove a baseline before analysis
  #INPUT ---
  #data: continuous df or list of trials, returns data in the same format
  #end: last sample of the baseline (range to remove)
  #start: first sample of the baseline, defaults to first sample of the trial
  #relabel: if True, changes sample numbering accordingly
  
  data = data.check(data, aslist=F) #transform to df
  nsamples = data.get_samplenum(data)
  if ( nsamples < end) {
    warning( "The number of the last sample is smaller than the specified end." )
    end = nsamples
  }
  trialnum = nrow(data)/nsamples
  data = data[data[,1]<start | data[,1]>end,] #remove range
  #fix sample numbers to be continuous again:
  if (relabel) {
    data[,1] = rep(1:( nsamples-end+start-1 ), trialnum) 
  }
  return( data.check(data, aslist=T, strip=F, transform=.transformed) ) #transform to list if needed
}

data.resample <- function(data, old.srate, new.srate) {
  ## resample data by trial-wise linear interpolation
  #INPUT ---
  #data: continuous df or trial list
  #old.srate: current sampling rate
  #new.srate: desired sampling rate
  
  data = data.check(data, aslist=T, strip=F) #transform to list
  #don't interpolate across boundaries:
  resampled = lapply(data, function(trial) {
    nsamples = data.get_samplenum(trial)
    n = ceiling( nsamples / (old.srate/new.srate) ) #new number of samples
    data.frame(samples=1:n, outcome=trial[1:n,2], 
               sapply( apply(trial[,is.datacol(trial)], 2, approx, n=n), "[[", 2) )
    })
  return( data.check(resampled, aslist=F, transform=.transformed) ) #transform to df if needed
}

data.trials.slice <- function(data, window, overlap=0, split=T) {
  ## group consecutive data points according to the specified window size into separate slices
  ## i.e. data points from different trials but same time intervals are grouped together 
  ## does several control checks to prevent unwanted results
  #INPUT ---
  #data: continuous df or list of trials, expects ordered (within trials) sample number in 1st col
  #window: number of samples per slice
  #overlap: >= 0 number of overlapping samples per slice, i.e. to have a sliding time window
  #split: if True, data is split. Else, only indices to split with are returned
  #       "slides" and "splitidx" which can be used like this:
  #       data[ slides[ splitidx == 1 ], ] #to get the first slice
  #       or split(data[slides,], splitidx) #get all slices at once
  #       if overlap is < 1, split(data, splitidx) is sufficient
  #       if data is a list of trials, slides has the indices, splitidx the intervals, use:
  #       i.e. data[[1]][ slides[splitidx==1] ,] | lapply(data, function(d) d[ slides[splitidx==1] ,])
  #RETURNS ---
  #a list of which the elements correspond to the number of slices
  #within each list element is a df with the grouped samples across different trials
  
  data = data.check(data, aslist=F) #transform to df
  nsamples = data.get_samplenum(data)
  if ( window >= nsamples ) {
    stop( "Window must be smaller than the number of samples per trial." )
  }
  #check if overlap makes sense:
  if ( overlap >= window ) {
    stop( "Overlap must be smaller than window or you can't slide across the data." )
  }
  if (overlap < 1) { #create simple time bins without overlap
    #check if slices will be of equal length:
    frac = (nsamples/window) %% 1
    if ( frac > 0) {
      warning( paste("Slices are not of equal length. Last slice contains only", 
                     paste0(round(frac*100), "% of samples relative to the other slices.")) )
    }    
    #slice for a single trial:
    slice = findInterval(1:nsamples, c(seq(1, nsamples, by=window), Inf))
    if (split) return( split(data, slice) ) #slice is recycled over the df
    return( list(slides=1:nrow(data), splitidx=slice) )
  }
  #sliding window (with overlap):
  startidx = seq(1, nsamples-1, by=window-overlap) #get sequence of start indices for the slices
  endidx = seq(window, nsamples, by=window-overlap) #get sequence of end indices for the slices
  imbalance = length(startidx)-length(endidx)
  endidx = c(endidx, rep( nsamples, imbalance )) #append nsamples as many times as needed
  #sliding time window intervals:
  slide = unlist( mapply( seq, startidx, endidx ) ) #single trial sample indices
  trialnum = nrow(data)/nsamples #number of trials
  #slide indices for whole df:
  slides = rep(slide, trialnum) + #repeat slide for all trials and add...
    rep( nsamples, trialnum*length(slide) ) * #the max sample number (repeated for all samples and trials) ...
    rep( 0:(trialnum-1), each=length(slide) ) #multiplied by the trial position (0:N-1)
  splitidx = rep( 1:length(startidx), c( rep( window, length(startidx)-imbalance ), #intervals to split with
                                         nsamples+1-rev(rev(startidx)[0:imbalance]) ) ) #correct for imbalance in last slide(s), if any
  if (imbalance) {
    warning( paste("Sliding Window: Your last", imbalance, "slice(s) contain(s)", 
                   "less samples than the previous", length(endidx)-imbalance, "slices.") )    
  }
  if (split) {
    temp = data[slides,] #create a new df with updated length due to overlap
    return( split(temp, splitidx) ) 
  }
  return( list(slides=slides, splitidx=splitidx) )
}

.data.unslice <- function(data, params.out=F) {
  ## unslice data that has been sliced (output as df)
  window = sapply( data, function(d) length( unique(d[,1]) ) )
  overlap = sum( unique( data[[1]][,1] ) %in% unique( data[[2]][,1] ) )
  slicelens = sapply( data, nrow )
  trialnum = slicelens[1]/window[1]
  fullidx = sum( window >= window[1] )
  
  unsliceidx = c( rep( rep(T, window[1]), trialnum ), #everything from 1st slice
                  #take the non-overlapping samples of all following "full" slices
                  rep( rep( c( rep(F, overlap), rep(T, window[1]-overlap) ), trialnum ), fullidx-1 ) )
  #in case of imbalance:
  if ( fullidx < length(data) ) {
    unsliceidx = c( unsliceidx, 
                    #add first imbalanced slice with different window size
                    rep( c( rep(F, overlap), rep(T, window[fullidx+1]-overlap) ), trialnum ) )
  }
  data = as.data.frame( data.table::rbindlist(data) )
  if ( nrow(data) > length(unsliceidx) ) {
    #concatenate FALSE for the rest of the slices if there were more
    unsliceidx = c( unsliceidx, rep( F, nrow(data)-length(unsliceidx) ) )
  }
  #now subset to get back original data:
  data = data[unsliceidx,]
  #regroup samples in the correct order corresponding to the trials
  trialsplit = c( rep(1:trialnum, each=window[1]), #first slice full window
                  rep( rep(1:trialnum, each=window[1]-overlap), fullidx-1 )) #other slices non-overlapping samples
  if ( fullidx < length(window) ) {
    #add last (imbalanced) slice samples:
    trialsplit = c( trialsplit, rep(1:trialnum, each=window[fullidx+1]-overlap) )
  }
  #now get trial data format:
  data = split(data, trialsplit)
  data = as.data.frame( data.table::rbindlist(data), row.names=NULL)
  if ( !params.out )  return( data )
  return( list(data=data, window=unname(window[1]), overlap=overlap) )
}

data.reduce_levels <- function(data, labels, col=2) {
  ## function to reduce the outcome factor to fewer levels
  ## subset multiclass data into binary data according to specified labels for 1 vs 1 classification
  #INPUT ---
  #data: continuous df or list of trials
  #labels: multi element vector corresponding to the class labels you wish to retain
  #col: column with the outcome (class labels), defaults to 2nd column
  #RETURNS ---
  #a subset of the original data with only the specified class labels
  
  data = data.check(data, aslist=F) #transform to df
  #subset data according to specified labels:
  data = data[ data[,col] %in% labels, ]
  #in case the outcome was a factor:
  try( data[,col] <- droplevels(data[,col]), silent=T ) 
  return( data.check(data, aslist=T, strip=F, transform=.transformed) ) #transform to list if needed
}

data.collapse_levels <- function(data, labels, col=2) {
  ## function to collapse specified factor labels
  ## preferable choice over data.reduce_levels if you want to keep all samples 
  ## and some factor levels are more similar than others
  #INPUT ---
  #data: continuous df or list of trials
  #labels: multi element vector or list with multi-element vectors
  #        corresponding to the class labels you wish to merge,
  #        e.g. c(1,2) changes the labels to 1
  #        e.g. list( c(1,2), c(3,4) ) changes labels to 1 and 3 (4 labels->2 labels)
  #col: column with the outcome (class labels), defaults to 2nd column
  #RETURNS ---
  #data of the same size with reduced number of class labels
  #note: label 0 will be avoided as a label for the merged classes

  data = data.check(data, aslist=F) #transform to df
  outcome = data[,col]
  if ( !is.list(labels) ) { labels = list(labels) }
  allLabs = do.call(c, labels ) #all specified labels
  if ( length( unique(allLabs) ) < length(allLabs) ) {
    stop( "Don't specify labels multiple times. This will likely lead to unintended results." )
  }  
  temp = rep(0, length(outcome)) #vector to save positions and labels
  for (lab in labels) {
    if (length(lab) < 2) {
      stop( paste0("Problem with label ", as.character(lab), 
                   ". Specify at least 2 labels to collapse.") )
    }
    collapsed = as.numeric( outcome %in% lab ) #merge labels
    newlab = ifelse(lab[1] > 0, lab[1], lab[2]) #new label for the combination (not 0)
    collapsed[collapsed>0] = newlab
    print( paste("Collapsed factor levels", 
                 paste(as.character(lab),collapse=" & "), "to", newlab) )
    temp = temp + collapsed #save position and label in temp
  }
  #change labels in data
  data[ outcome %in% allLabs, col] = temp[ outcome %in% allLabs ]
  try( data[,col] <- droplevels(data[,col]), silent=T ) #in case of factor
  return( data.check(data, aslist=T, strip=F, transform=.transformed) ) #transform to list if needed
}  


.eval_ellipsis <- function(funStr, ...) {
  ## helper to evaluate input to ellipsis (...) before handing it to a function
  require(utils, quietly=T)
  args.in = list(...)
  fun.args = formals( funStr ) #function arguments
  use.args = modifyList(fun.args, args.in) #overwrite matching arguments 
  use.args = use.args[ names(use.args) %in% names(fun.args) ] #use only legal arguments
  return( use.args )
}

.eval_package <- function(pkgStrings, load=F) {
  ## helper to check if package(s) are present
  for (pkgStr in pkgStrings) {
    #check if pkg exists
    if ( requireNamespace(pkgStr, quietly=T) ) { 
    } else { 
      #install pkg if not present
      install.packages(pkgStr)
    }
    if (load) { #attach pkg
      suppressMessages( require(pkgStr, quietly=T, character.only=T) )
    }
  }
}

.parallel_backend <- function(on = T, CPUcluster = NULL, nCores = NULL) {
  ## function to start or end parallel backend
  if (on) {
    .eval_package( c("parallel", "doParallel") )
    if (is.null(nCores)) nCores = parallel::detectCores()
    CPUcluster = parallel::makeCluster(nCores)
    doParallel::registerDoParallel(CPUcluster)
    return( CPUcluster )
  } else {
    parallel::stopCluster(CPUcluster)
    #unregister foreach variable scope
    env = foreach:::.foreachGlobals
    rm(list=ls(name=env), pos=env) 
    invisible( gc() )
  }
}
