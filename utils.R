## UTILITIES: often needed helper functions

is.datacol <- function(data, tol = .Machine$double.eps^0.5) {
  ## check if columns of data are measurements (floating numbers)
  ## used to exclude info columns in a df (e.g. sample number, outcome)
  if ( is.null(dim(data)) ) data = as.matrix(data)
  sapply( 1:ncol(data), function(i) {
    if ( is.numeric( data[,i] ) ) {
      #are there non-integer values in the column?
      return( sum( abs( data[,i] - round(data[,i]) ) < tol, na.rm=T ) < nrow(data) ) 
    }
    return( FALSE ) #non-numeric
  })
}

data.set_type <- function(data) {
  ## automatic type detection for lists. non-lists are returned without action
  ## if type is already set, returns immediately
  ## distinction can be problematic for lists of slices and subjects
  ## current decision is based on the following criteria:
  # trial list: 
  #   - if the unique sample numbers match nrow in the 1st list element
  # subject list: 
  #   - if overlapping samples between 1st and 2nd list element correspond to nrow
  #   - or if the outcome is in non-identical order (highly likely)
  # slice list: 
  #   - overlapping samples must be smaller than sample number,
  #   - outcome should be in identical order between slices
  ## the outcome is checked as well because sample numbering might vary across
  ## subject data sets (for whatever reason); in the very unusual scenario
  ## of an identical outcome order between the 1st and 2nd subject, the list
  ## will be fasely identified as slices and should instead be set manually
  ## likewise, a slice list which was created by hand should have constant outcome order
  ## otherwise it will be falsely identified as subjects
  if ( any( c("trials", "slices", "subjects") %in% attr(data, "type") ) ) {
    return(data)
  }
  if ( class(data) == "list" ) {
    if ( class(data[[1]]) == "list" ) { #nested list, perhaps subject data
      data = lapply( data, data.check, aslist=F )
    }
    #get sample number of 1st and 2nd list element
    nsamples = c( data.samplenum(data[[1]]), data.samplenum(data[[2]]) )
    if ( identical( nrow(data[[1]]), nsamples[1] ) ) { #trial list
      #unique sample numbers match row number
      attr(data, "type") = "trials" 
    } else { #slice or subject list
      #if slice list, "overlap" cannot be identical to sample number
      #because subject data could have non-overlapping sample numbers, 
      #also check if outcome has the same order in 1st and 2nd element (unlikely for subject data)
      overlap = sum( unique( data[[1]][,1] ) %in% unique( data[[2]][,1] ) )
      if ( overlap %in% nsamples || !identical( data[[1]][seq( 1, nrow(data[[1]]), by=nsamples[1] ), 2], 
                                                data[[2]][seq( 1, nrow(data[[2]]), by=nsamples[2] ), 2] ) ) {
        attr(data, "type") = "subjects"
      } else { #slice list
        attr(data, "type") = "slices"
      }
    }
  }
  return(data)
}

data.set_info <- function(data, trials=NULL, outcome=NULL) {
  ## function to append required info columns (sample numbers and/or outcome) to raw data
  ## which only contains measurements (or at least not the required info cols)
  #INPUT ---
  #data: df or list of trials, slices, subjects
  #trials: int, if specified sample numbers will be added to 1st col accordingly
  #outcome: vector with 1 value per trial, if specified will be added to 2nd col
  #NOTES ---
  #if only one of samples or outcome should be added, leave the other at NULL
  #for outcome this requires the presence of sample numbers in the 1st col
  data = data.set_type(data)
  if ( "subjects" %in% attr(data, "type") ) {
    if ( !is.null(outcome) ) {
      warning( "You are appending identical outcome info to every subject data set in the list." )
    }
    data = lapply( data, data.set_info, trials=trials, outcome=outcome )
    attr(data, "type") = "subjects"
    return(data)
  }
  data = data.check(data, aslist=F)
  if ( !is.null(outcome) ) {
    if ( !is.null(trials) ) { #trials specified
      if ( (nrow(data)/trials) %% 1 > 0 ) {
        stop( "Either your specified trial number is incorrect or your trials differ in length." )
      }
      if ( length(outcome) != trials ) {
        stop( "Length of the outcome must match number of trials." )
      }
      outcome = as.factor( rep( outcome, each=nrow(data)/trials ) )
      if ( length(outcome) < nrow(data) ) {
        stop( "Number of trials does not match the length of the data." )
      }
      data = cbind(outcome, data)
    } else if ( is.datacol(data[,1]) ) { #sample number is not present in df and no trials specified
      stop( "No sample numbers found in 1st column and no trials specified.\n", 
            "Without proper trial info, outcome cannot be appended." )
    } else { #sample number is present in df
      nsamples = data.samplenum(data)
      outcome = as.factor( rep( outcome, each=nsamples ) )
      if ( length(outcome) < nrow(data) ) {
        stop( "Length of the outcome must match number of trials." )
      }
      data = cbind(samples=data[,1], outcome, data[,-1]) #append to 2nd col
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
    data = cbind(samples, data)
  }
  return( data.check(data, aslist=T, strip=F, transform=.transformed) )
}

data.samplenum <- function(data) {
  ## helper function to retrieve sample numbers per trial
  ## data is expected to have trials of equal length and ordering
  if ( is.datacol(data[,1]) ) {
    stop( "It seems your data does not have sample numbers in its first column.\n", 
          "Please provide that information. Use data.set_info at your convenience." )
  }
  samples = unique(data[,1]); nsamples = length(samples)
  if ( !isTRUE(all.equal( rep( samples, nrow(data)/nsamples ), data[,1] )) ) {
    stop( "Perhaps your first column does not contain sample numbers.\n",
          "Other possible causes may be that the order of your sample numbers varies,\n", 
          "or your trials are not of equal length. Please meet these assumptions." )
  }
  return( nsamples )
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
  data = data.set_type(data) #set the list type if necessary
  if (aslist) {
    if ( class(data) != "list" ) { #if not a list
      .transformed <<- TRUE
      return( data.split_trials(data, strip=strip) ) #transform to trial list 
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
    if ( "trials" %in% attr(data, "type") ) {
      .transformed <<- TRUE #retransform notice only if it was a trial list
      data = data.table::rbindlist(data) #transform list to df
    } else if ( "slices" %in% attr(data, "type") ) {
      data = .data.unslice(data)
    } else if ( "subjects" %in% attr(data, "type") ) {
      stop( "Loop the individual data sets of your subject data list for this function." )
    }
  }
  #df, dt or matrix/vector to df
  data = as.data.frame(data)
  if (strip) return( data[, is.datacol(data)] )
  return(data)  
}

data.split_trials <- function(data, strip=F) {
  ## transform data into a list of trials according to sample numbering (1st col)
  ## returns a list of dfs with each element corresponding to a single trial
  ## sets the attribute "type" to "trials"
  #INPUT ---
  #strip: non-numeric columns are stripped (e.g. sample numbers, outcome)
  
  ## make sure data is a df 
  # can't call data.check because it might interfere with .transformed of previous call!
  data = data.set_type(data)
  if ( "trials" %in% attr(data, "type") ) return(data)
  if ( "slices" %in% attr(data, "type") ) {
    data = .data.unslice(data)
  } else if ( "subjects" %in% attr(data, "type") ) {
    data = lapply(data, data.split_trials, strip=strip) #iterate the subjects
    attr(data, "type") = "subjects"
    return(data)
  }
  data = as.data.frame(data) #in case of matrix
  ##
  nsamples = data.samplenum(data) #1st col with sample info required
  trials = findInterval(1:nrow(data), seq(1, nrow(data), by=nsamples))  
  if (strip) {
    trialdata = split(data[, is.datacol(data)], trials)
  } else {
    trialdata = split(data, trials)
  }
  attr(trialdata, "type") = "trials"
  return( trialdata )
}

data.split_slices <- function(data, window, overlap=0, split=T) {
  ## group consecutive data points according to the specified window size into separate slices
  ## i.e. data points from different trials but same time intervals are grouped together 
  ## does several control checks to prevent unwanted results
  #INPUT ---
  #data: continuous df or list of trials, subjects
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
  data = data.set_type(data)
  if ( "subjects" %in% attr(data, "type") ) {
    data = lapply( data, data.split_slices, window=window, overlap=overlap, split=split )
    if (split)  attr(data, "type") = "subjects"
    return(data)
  }
  data = data.check(data, aslist=F) #transform to df
  nsamples = data.samplenum(data)
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
      cat( "Slices are not of equal length. Last slice contains only ", 
           round(frac*100), "% of the samples relative to the other slices.\n", sep="" )
    }    
    #slice for a single trial:
    slice = findInterval(1:nsamples, c(seq(1, nsamples, by=window), Inf))
    if (split) {
      slicedata = split(data, slice) #slice is recycled over the df
      attr(slicedata, "type") = "slices"
      return( slicedata )
    }
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
    cat( "Sliding Window: Your last", imbalance, "slice(s) contain(s)", 
         "less samples than the previous", length(endidx)-imbalance, "slices.\n" )    
  }
  if (split) {
    temp = data[slides,] #create a new df with updated length due to overlap
    slicedata = split(temp, splitidx)
    attr(slicedata, "type") = "slices"
    return( slicedata )
  }
  return( list(slides=slides, splitidx=splitidx) )
}

.data.unslice <- function(data) {
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
  attr(data, "window") = unname(window[1])
  attr(data, "overlap") = overlap
  return(data)
}

data.normalize <- function(data, start=NULL, end=NULL, scale=F) {
  ## remove the mean of a range of data points from the other samples of a trial
  ## generally used to subtract the average of a baseline period, i.e. 1:baseline_end
  #INPUT ---
  #data: continuous df or list of trials, returns data in the same format
  #start: first sample of the range to subtract, defaults to first sample of the trial
  #end: last sample of the range, defaults to last sample of the trial
  #scale: if True, data is both centered and scaled, otherwise only centered
  data = data.check(data, aslist=T, strip=F) #transform to list
  if ( "slices" %in% attr(data, "type") ) {
    data = data.split_trials( data, strip=F )
    .transformed = T #convert to df at the end
  } else if ( "subjects" %in% attr(data, "type") ) {
    data = lapply(data, data.normalize, start=start, end=end, scale=scale)
    attr(data, "type") = "subjects"
    return(data)
  }
  nsamples = data.samplenum(data[[1]])
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

data.remove_samples <- function(data, start=1, end, relabel=T) {
  ## remove datapoints from trials, specified by sample number (1st col)
  ## can be used e.g. to remove a baseline before analysis
  #INPUT ---
  #data: continuous df or list of trials, slices, subjects
  #start: first sample to remove, defaults to first sample of the trial
  #end: last sample to remove (end of the range to remove)
  #start and end correspond to positions and not to the actual values in the column
  #relabel: if True, adapts 1st col to the new sample number (from 1:last sample)
  data = data.set_type(data)
  if ( "subjects" %in% attr(data, "type") ) {
    data = lapply( data, data.remove_samples, start=start, end=end, relabel=relabel )
    attr(data, "type") = "subjects"
    return(data)
  }
  data = data.check(data, aslist=F) #transform to df
  nsamples = data.samplenum(data)
  if ( nsamples <= 1 ) stop( "Only 1 sample per trial found." )
  if ( nsamples < end ) {
    warning( "The specified end exceeds the sample numbers. ",
             "Setting end to the last sample." )
    end = nsamples
  }
  trialnum = nrow(data)/nsamples
  #create idx in case sample numbering is discontinuous
  sampleidx = rep( 1:nsamples, trialnum )
  data = data[sampleidx<start | sampleidx>end,] #remove range
  if (relabel) { #fix sample numbers to be continuous again
    data[,1] = rep(1:( nsamples-end+start-1 ), trialnum) 
  }
  return( data.check(data, aslist=T, strip=F, transform=.transformed) ) #transform to list if needed
}

data.remove_classes <- function(data, labels, col=2) {
  ## function to reduce the outcome factor to fewer levels
  ## subset multiclass data into binary data according to specified labels for 1 vs 1 classification
  #INPUT ---
  #data: continuous df or list of trials, slices, subjects
  #labels: multi element vector corresponding to the class labels you wish to retain
  #col: column with the outcome (class labels), defaults to 2nd column
  #RETURNS ---
  #a subset of the original data with only the specified class labels
  data = data.set_type(data)
  if ( "subjects" %in% attr(data, "type") ) {
    data = lapply( data, data.remove_classes, labels=labels, col=col )
    attr(data, "type") = "subjects"
    return(data)
  }
  data = data.check(data, aslist=F) #transform to df
  #subset data according to specified labels:
  data = data[ data[,col] %in% labels, ]
  #in case the outcome was a factor:
  try( data[,col] <- droplevels(data[,col]), silent=T ) 
  return( data.check(data, aslist=T, strip=F, transform=.transformed) ) #transform to list if needed
}

data.merge_classes <- function(data, labels, new_labels=NULL, col=2, verbose=T) {
  ## function to collapse specified factor labels
  ## preferable choice over data.remove_classes if you want to keep all samples 
  ## and some factor levels are more similar than others
  #INPUT ---
  #data: continuous df or list of trials, slices, subjects
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
  data = data.set_type(data)
  if ( "subjects" %in% attr(data, "type") ) {
    data = lapply( data, data.merge_classes, labels=labels, new_labels=new_labels, col=col, verbose=verbose )
    attr(data, "type") = "subjects"
    return(data)
  }
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
    if (verbose) cat("Merged classes", paste(as.character(lab),collapse=" & "), "to", newlab,"\n")
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
  #INPUT:
  #data: df or list of trials, slices, subjects
  #Q, IQR, C: threshold parameters
  #plot: if True, diagnostic plot is printed which plots the trial variance and the
  #      threshold (red dashed line). Removed trials are marked with their trial number
  #RETURNS:
  #data without the outliers,
  #has the attributes: threshold, removed (trials), variance (1 value per trial)
  if ( any( findInterval( c(Q,IQR), c(0,100) ) != 1 ) ) {
    stop( "Q and IQR have to be in the range 0 to 100." )
  }
  data = data.set_type(data)
  if ( "slices" %in% attr(data, "type") ) {
    trialdata = data.split_trials( data, strip=T )
    .transformed = T #convert to df at the end
  } else if ( "subjects" %in% attr(data, "type") ) {
    data = lapply(data, data.remove_outliers, Q=Q, IQR=IQR, C=C, plot=plot)
    attr(data, "type") = "subjects"
    return(data)
  } else { #trial list or df
    trialdata = data.check(data, aslist=T, strip=T)  
  }
  # compute global variance values
  trialvar = colMeans( sapply( trialdata, function(trial) apply(trial, 2, var, na.rm=T) ) )
  # define threshold
  IQR = c( 0+(100-IQR)/200, 1-(100-IQR)/200 ) #convert IQR values to 0-1 range
  lims = quantile( trialvar, probs=c(Q/100, IQR), na.rm=T ) #get threshold values
  threshold = unname( lims[1] + C * (lims[3]-lims[2]) )
  outliers = which( unname(trialvar > threshold) ) #find trials that exceed global variance threshold
  cat( "Removed ", length(outliers), " outliers (",
       round(length(outliers)/length(trialvar)*100, 2), "% of trials).\n", sep="" )
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
  # return data without the outlier trials + info as attributes
  data = data.split_trials(data, strip=F) #get trialdata with info columns
  if ( length(outliers) < 1 ) { #no outliers
    outliers = 0
    data = data.check( data, aslist=F, transform=.transformed )
  } else { #remove outliers
    data = data.check( data[-outliers], aslist=F, transform=.transformed )
  }
  attr(data, "removed") = outliers
  attr(data, "threshold") = threshold
  attr(data, "variance") = trialvar
  return(data)
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

.parallel_check <- function(required=NULL, nCores=NULL, output=NULL) {
  ## helper function to determine if a parallel backend is already running,
  ## if one should be initiated, or if the user requested a sequential backend
  #INPUT:
  #required: maximum number of required cores, 
  #          e.g. number of subjects to parallelize over, i.e. length(data)
  #nCores: user input to nCores. If NULL, automatic selection of cores
  #        between maximum required and available cores
  #        can be a cluster object (makeCluster) to use a registered backend
  #output: previously generated output by this function to shut off backend (if necessary)
  #RETURNS:
  #if output is NULL, an output list with CPUcluster and nCores is generated
  #in a sequential backend, CPUcluster will be an empty list
  .eval_package("foreach", load=T) #parallelize for subject list
  if ( is.null(output) ) {
    output = list( CPUcluster = list(), outside = F, start.time = proc.time()[3] )
    if ( is.list(nCores) && getDoParRegistered() ) { #cluster input
      output$CPUcluster = nCores
      output$nCores = getDoParWorkers() 
      output$outside = T #cluster was created outside and will not be shut down
    } else if ( !is.null(nCores) && nCores <= 1 ) { #no parallelization
      registerDoSEQ() #sequential backend for foreach
      output$nCores = 1
    } else {
      if ( is.null(nCores) ) { #default to as many cores as available/required
        output$nCores = min( parallel::detectCores(), required )
      } else { #nCores was set by user
        #sanity check user input: only as many cores as available/required
        output$nCores = min( min( nCores, parallel::detectCores() ), required )
      } 
      output$CPUcluster = .parallel_backend(on=T, nCores=output$nCores) #initiate cores
    }
    #TEMP: export functions to clusters
    if ( length(output$CPUcluster) > 0 ) {
      parallel::clusterEvalQ( output$CPUcluster, expr={ source("decode.R") } )
    }
    cat( "Distributing", required, "jobs to", output$nCores, "workers . . .\n" )
    return(output) #list with CPUcluster and nCores
  } else if ( !output$outside ) { #if cluster was created by function call...
    if ( length(output$CPUcluster) == 0 && output$nCores == 1 ) { 
      #remove sequential backend
      env = foreach:::.foreachGlobals
      rm(list=ls(name=env), pos=env) 
    } else if ( length(output$CPUcluster) > 0 ) {
      .parallel_backend(on=F, output$CPUcluster) #shut down cores
    } #else externally created cluster which is left as is
  }
  cat( "Done. Time elapsed (mins):", (proc.time()[3]-output$start.time)/60, "\n" )
}