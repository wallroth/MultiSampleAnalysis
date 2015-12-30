# utilities: often needed helper functions

is.datacol <- function(data, tol = .Machine$double.eps^0.5) {
  ## check if columns of data are measurements (floating numbers)
  ## used to exclude info columns in a df (e.g. sample number, outcome)
  if ( is.null(dim(data)) ) data = as.matrix(data)
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
  #strip: remove non-data (measurement) columns, see is.datacol
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
  } #else... 
  if ( class(data) == "list" ) { #if a list...
    #find out if trial list or slice list:
    if ( .is.sliced(data) ) {
      data = .data.unslice(data)
    } else {
      #retransform notice only if it was a trial list:
      .transformed <<- TRUE
      data = data.table::rbindlist(data) #transform list to df
    }
  }
  #df, dt or matrix/vector to df
  data = as.data.frame(data)
  if (strip) return( data[, is.datacol(data)] )
  return(data)  
}

data.trials.normalize <- function(data, start=NULL, end=NULL, scale=F) {
  ## remove the mean of a range of data points from the other samples of a trial
  ## generally used to subtract the average of a baseline period, i.e. 1:baseline_end
  #INPUT ---
  #data: continuous df or list of trials, returns data in the same format
  #start: first sample of the range to subtract, defaults to first sample of the trial
  #end: last sample of the range, defaults to last sample of the trial
  #scale: if True, data is both centered and scaled, otherwise only centered
  data = data.check(data, aslist=T, strip=F) #transform to list
  nsamples = data.get_samplenum(data[[1]])
  if ( is.null(start) || start < 1 ) start = 1
  if ( is.null(end) || end > nsamples ) end = nsamples
  dcol = is.datacol(data[[1]])
  normdata = lapply(data, function(trial) {
    normalized = scale(trial[start:end, dcol], scale=scale)
    if ( !identical(end, nsamples) ) {
      #normalize only via a baseline
      if (scale) {
        temp = scale(trial[ c( 0:(start-1), (end+1):nrow(trial) ), dcol], 
                     attr(normalized,"scaled:center"), attr(normalized,"scaled:scale"))
      } else {
        temp = scale(trial[ c( 0:(start-1), (end+1):nrow(trial) ), dcol], 
                     attr(normalized,"scaled:center"), scale=F)
      }
      normalized =  rbind( temp[0:(start-1), , drop=F], #pre-baseline samples
                           normalized, #baseline samples
                           temp[start:nrow(temp), , drop=F] ) #post-baseline samples
    }
    setNames( data.frame( trial[, !dcol], normalized ), names(trial) ) #return to lapply
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
  ## uses zero padding so that artifacts won't be carried into the data
  #INPUT ---
  #data: continuous df or trial list
  #old.srate: current sampling rate
  #new.srate: desired sampling rate
  .eval_package("MASS")
  data = data.check(data, aslist=T, strip=F) #transform to list
  nsamples = data.get_samplenum(data[[1]])
  frac = MASS::fractions(new.srate/old.srate)
  #calculate zero pad length:
  #note: by using fraction, rounding error accumulation is minimized
  if ( old.srate %% new.srate != 0 ) { 
    temp = strsplit( attr(frac, "fracs"), "/" )
    p = as.numeric( temp[[1]][1] )
    q = as.numeric( temp[[1]][2] )
  } else { #not a fraction
    p = new.srate/old.srate 
    q = 1
  }
  nPad = ceiling( (max(p,q) * 10)/q ) * q
  n = ceiling( (2*nPad + nsamples)*frac ) #new number of samples
  unPad = nPad * frac
  #don't interpolate across boundaries:
  resampled = lapply(data, function(trial) {
    temp = trial[, is.datacol(trial) ]
    padstart = temp[rep(1, nPad), ] #repeat first row n times
    padend = temp[rep(nrow(temp), nPad), ] #repeat last row n times
    temp = rbind(padstart, temp, padend) #zero padded trial data
    #resample
    temp = sapply( apply(temp, 2, approx, n=n), "[[", 2)
    data.frame(samples=1:(n-2*unPad), outcome=trial[1,2], temp[(unPad+1):(nrow(temp)-unPad),])
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

data.collapse_levels <- function(data, labels, new_labels=NULL, col=2, verbose=T) {
  ## function to collapse specified factor labels
  ## preferable choice over data.reduce_levels if you want to keep all samples 
  ## and some factor levels are more similar than others
  #INPUT ---
  #data: continuous df or list of trials
  #labels: multi element vector or list with multi-element vectors
  #        corresponding to the class labels you wish to merge,
  #        e.g. c(1,2) changes the labels to 1
  #        e.g. list( c(1,2), c(3,4) ) changes labels to 1 and 3 (4 labels->2 labels)
  #new_labels: a vector defining new labels 
  #col: column with the outcome (class labels), defaults to 2nd column
  #verbose: print information about collapsed levels and new names
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
  temp = as.character(outcome)
  for ( i in seq_along(labels) ) {
    lab = labels[[i]]
    if ( !is.null(new_labels) && i <= length(labels) ) {
      newlab = new_labels[i]
    } else {
      newlab = lab[1]
    }
    temp[ outcome %in% lab ] = newlab
    if (verbose) {
      print( paste("Collapsed factor levels", 
                   paste(as.character(lab),collapse=" & "), "to", newlab) )
    }
  }
  #change labels in data
  data[, col] = as.factor(temp)
  try( data[,col] <- droplevels(data[,col]), silent=T )
  return( data.check(data, aslist=T, strip=F, transform=.transformed) ) #transform to list if needed
}  

data.remove_outliers <- function(data, Q=50, IQR=90, C=3, plot=F) {
  ## remove trials that are classified as outliers based on variance
  ## the variance is computed within each trial for each channel,
  ## then the channel variances are averaged, yielding one value per trial
  ## finally the trials with global values that exceed the threshold are excluded
  ## the threshold is defined as: Q + C * IQR
  if ( any( findInterval( c(Q,IQR), c(0,100) ) != 1 ) ) {
    stop( "Q and IQR have to be in the range 0 to 100." )
  }
  trialdata = data.check(data, aslist=T, strip=T)
  # compute global variance values
  trialvar = colMeans( sapply( trialdata, function(trial) apply(trial, 2, var, na.rm=T) ) )
  # define threshold
  IQR = c( 0+(100-IQR)/200, 1-(100-IQR)/200 ) #convert IQR values to 0-1 range
  lims = quantile( trialvar, c(Q/100, IQR) ) #get threshold values
  threshold = unname( lims[1] + C * (lims[3]-lims[2]) )
  outliers = which( unname(trialvar > threshold) ) #find trials that exceed global variance threshold
  print( paste0("Removed ", length(outliers), " outliers (", 
               length(outliers)/length(trialvar)*100, "% of trials).") )
  if (plot) {
    #visualize the outliers
    ylims = c( min( min(trialvar), threshold ), max( max(trialvar), threshold ) )
    plot(trialvar, xlab="Trial", ylab="Variance", las=1, ylim = ylims)
    abline(h = threshold, col="red", lty=2)
    if ( length(outliers) > 0 ) {
      outlier_labels = rep("", length(trialvar))
      outlier_labels[outliers] = as.character(outliers)
      text(trialvar, labels=outlier_labels, cex=.7, pos=4)  
    }
  }
  if ( length(outliers) < 1 ) return( list(data = data, removed = 0, variance=trialvar, threshold=threshold) )
  # prepare output:
  # return data without the outlier trials and info which trials were excluded
  trialdata = data.check(data, aslist=T, strip=F) #get trialdata with info columns
  return( list(data = data.check( trialdata[-outliers], aslist=F, transform=.transformed ), 
               removed = outliers, variance=trialvar, threshold=threshold) )
}



.eval_ellipsis <- function(funStr, ...) {
  ## helper to evaluate input to ellipsis (...) before handing it to a function
  require(utils, quietly=T)
  args.in = list(...)
  fun.args = formals( funStr ) #function arguments
  use.args = modifyList(fun.args, args.in, keep.null=T) #overwrite matching arguments 
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
