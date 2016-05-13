library(data.table)
#modifications are done globally; use copy(dt) to make local modifications
#index a column via character with get(..)
#prevent dropping of columns by specifying in by: x[ outcome == "A" , lapply(.SD, mean), .SDcols=4:9, by=.(trial, outcome) ]
#CJ to select multiple values: dat[CJ(1:2, 1:3)]
#instead of dat[x==1 | y ==3] better to: dat[.(1,3)]
#select by multiple keys: dat[.(1:5, 5, "B"), nomatch=0]
#last row: x[.N], second last: x[.N-1]
#x[,.N] returns the number of rows
# y[, lapply(.SD, scale), .SDcols = -(1:3) ] #apply scale directly
# assign: y[, (names(y)[-c(1:3)]) := lapply(.SD, scale), .SDcols = -(1:3) ]
# scale( y[, .SD, .SDcols = -(1:3) ] ) #get matrix output with attributes
# insert column-wise: for ( col in 1:ncol(yy) ) set(y, i=NULL, j=colnames(yy)[col], yy[,col])
# select non-key columns: setdiff(names(x),key(x))
# DT[c("a","b"),sum(v),by=.EACHI]   # sum for two groups of v
# DT[,sum(v),by=.(y%%2)]   # by expression
#get the actual row index: data[, .I[sample %in% baseline] ]
#data[.(c(1,2,1),c(2,3,5))] subset for subject 1 trial 2 and 5 and for subject 2 trial 3
#data[.(unique(subject),1)] subset for every subject trial 1


read.data <- function(path=".", pattern="*.csv", nCores=NULL, ...) {
  ## function to read in data files while also setting the necessary info attributes
  #INPUT ---
  #path: string specifying the location (folder) with the data files
  #pattern: pattern common to all data files, e.g. the file ending
  #nCores: number of CPU cores to use for parallelization
  #        if NULL, automatic selection; if 1, sequential execution
  #        if an empty list, an externally registered cluster will be used
  #...: further arguments to fread and to data.setinfo
  #     important arguments are 'samples' (which has to be constant for all data sets) 
  #     and 'outcome' which can only be specified as an integer that indicates the column
  #     with the outcome values for each data file. If no such column is present you can
  #     specify NULL and add the outcome afterwards by looping the sets.
  #RETURNS ---
  #data.table with info columns, accepted by all functions (if outcome was not NULL)
  .eval_package(c("foreach","data.table"), load=T)
  args = .eval_ellipsis(c("data.setinfo", "fread"), ...)
  if ( identical(args$fread$dec, .eval_ellipsis("fread")$fread$dec) ) {
    args$fread$dec = "." #do.call on fread may fail due to default of 'dec'
  }
  files = list.files(path=path, pattern=pattern, full.names=T)
  cluster = parBackend(cores=nCores, required=length(files)) #parallelize read in
  data = foreach(f=files, .combine = list, .multicombine = T, .maxcombine=length(files)+1) %dopar%
  {
    d = do.call(fread, modifyList(args$fread, list(input=f)) ) #read file
    do.call( data.setinfo, modifyList(args$data.setinfo, list(data=d)) ) #set info attributes
  }
  parBackend(cores=cluster) #shut down cluster (if created by inside call)
  if ( length(files) > 1 ) { #list of multiple subjects
    if ( all(sapply(data, function(d) d[, unique(subject) == 0L])) ) { #subject info was not set
      data = rbindlist(lapply(data, function(d) d[, subject:=NULL]), idcol="subject", use.names=T)
    } else {
      data = rbindlist(data, use.names=T)
    }
  }
  return(data.check(data)[])
}

is.float <- function(x) if (is.numeric(x)) any(x %% 1 > 0) else FALSE

data.setinfo <- function(data, .outcome=NULL, .subject=NULL, .samples=NULL, .trials=NULL, strip=F, asFactor=T) {
  #outcome: character or number specifying a column, or vector corresponding to trial number or vector corresponding to nrow
  #subject: character or integer specifying a column, or numeric (non-integer) specifying the subject id
  #samples: character or integer specifying a column, or numeric (non-integer) specifying the number of samples per trial
  #trials: character or integer specifying a column, or numeric (non-integer) specifying the number of samples per trial
  #strip: logical, if True, non-floating number columns (excluding the info columns) are removed from data
  data = as.data.table(data) #create a copy
  #begin by checking for non-NULL input:
  if ( is.character(.subject) || is.integer(.subject) ) { #supplied as column name/number
    setnames(data, .subject, "subject") #change column name to subject
  } else if ( is.numeric(.subject) ) { #supplied as id
    data[, subject := .subject] #add id as column
  }
  if ( is.character(.samples) || is.integer(.samples) ) { #supplied as column name/number
    setnames(data, .samples, "sample") #change column name to sample
  } else if ( is.numeric(.samples) ) { #supplied as quantity
    if (is.null(.trials)) cat( "Non-integer input to .samples detected:", .samples, "samples assigned to each of", nrow(data)/.samples, "trials.\n" )
    data[, sample := rep(1:.samples, times=.N/.samples)] #repeat for each trial
  }
  if ( is.character(.trials) || is.integer(.trials) ) { #supplied as column name/number
    setnames(data, .trials, "trial") #change column name to trial
  } else if ( is.numeric(.trials) ) { #supplied as quantity
    if (is.null(.samples)) cat( "Non-integer input to .trials detected:", .trials, "trials extended for", nrow(data)/.trials, "samples each.\n" )
    data[, trial := rep(1:.trials, each=.N/.trials)] #repeat for number of samples
  }
  if ( length(.outcome) == 1 ) { #supplied as column name/number
    setnames(data, .outcome, "outcome")  #change column name to outcome
  } else if ( length(.outcome) == nrow(data) ) { #supplied as separate column (vector)
    data[, outcome := .outcome] #add column to data
  }
  #check if all columns are set:
  check = c("sample","trial","outcome","subject") %in% names(data)
  if ( !check[4] ) data[, subject := 0L] #placeholder
  if ( !all(check[1:3]) ) { #in case at least one is not set...
    #case 0|0|0
    if ( all(!check[1:3]) ) { #if none are set only outcome was supplied as vector with 1 value per trial
      .samples = nrow(data)/length(.outcome) #each outcome value corresponds to a single trial
      data[, ':=' ( trial = rep(1:(.N/.samples), each=.samples), 
                    sample = rep(1:.samples, times=.N/.samples), 
                    outcome = rep(.outcome, each=.samples)) ]
    #case 1|x|x
    } else if (check[1]) { #samples are set
      if (!check[2]) data[, trial := rep(1:(.N/uniqueN(sample)), each=uniqueN(sample)) ] #case 1|0|x, trials not set
      if (!check[3]) data[, outcome := rep(.outcome, each=uniqueN(sample)) ] #case 1|x|0, outcome not set
    #case 0|1|x
    } else if (check[2]) { #trials are set
      data[, sample := 1:.N, by=.(subject,trial)] #samples not set
      if (!check[3]) data[, outcome := .outcome[.GRP], by=.(subject,trial)] #case 0|1|0, outcome not set
    #case 0|0|1
    } else { #only outcome is set
      #estimate sample/trial number via minimum consecutive occurence of an outcome value (not fail-proof)
      .samples = data[, min(rle(as.character(outcome))$lengths)] #problematic in block designs or small data sets
      cat( "Estimating the sample number via minimum consecutive occurence of an outcome to be:", .samples, "\n" )
      data[, ':=' ( trial = rep(1:(.N/.samples), each=.samples), sample = rep(1:.samples, times=.N/.samples) ) ]
    }
  }
  if (asFactor) data[, outcome := ordered(outcome)]
  data = data.check(data) #sets key and sorts columns
  if (strip) { #remove non-float, non-info columns
    measurements = sapply( names(data)[-(1:4)], function(n) is.float(data[[n]]) )
    remove = which(!measurements)+4
    for ( n in names(data)[remove] ) set(data, j=n, value=NULL)
  }
  return(data[])
}

data.check <- function(data) {
  ## checks the data set for presence of info columns, column ordering and keys
  cn = c("subject", "trial", "sample", "outcome")
  if ( is.data.table(data) && all(cn %in% key(data)) ) return(data) #everything OK
  #else...
  if ( !all(cn %in% names(data)) ) stop( "Info columns missing. Please refer to data.setinfo." )
  if ( !is.data.table(data) ) data = as.data.table(data) #if not a data table
  if ( !all(cn %in% key(data)) ) setkeyv(data, cn) #if keys are not set
  if ( !identical(names(data)[1:length(cn)], cn) ) setcolorder( data, union(cn, names(data)) ) #order columns
  return(data) #no copy if it was a DT
}

data.normalize <- function(data, baseline, scale=F, copy=F) {
  #baseline: vector with sample numbers that represent baseline or single number setting endpoint of baseline
  #scale: perform scaling in addition to centering
  #copy: modify DT directly or create a copy. non-DT input will be copied necessarily
  if (copy) data = copy(data)
  data = data.check(data)
  if ( data[, uniqueN(sample)] <= 1 ) stop( "Normalization requires at least 2 samples per trial. ")
  if (length(baseline) == 1) baseline = 1:baseline
  vec.norm <- function(x, idx, scale) {
    x = x-mean(x[idx]) #subtract baseline mean from whole vector
    if (scale) x = x/sd(x[idx]) #divide whole vector by baseline sd
    return(x)
  }
  nm = setdiff( names(data), key(data) )
  equal = data[, .N, by=sample][, uniqueN(N) == 1] #identical distribution for all samples
  if (equal) { #generalize index of first subject/trial to all (fast)
    idx = data[.(subject[1],trial[1],baseline), which=T]
    data[, (nm) := lapply(.SD, vec.norm, idx=idx, scale=scale), .SDcols=nm, by=.(subject,trial) ]
  } else { #individual indices for every subject/trial, accessed via .I (slow, better way?)
    idxDT = data[, .(idx=sample %in% baseline), by=.(subject,trial)] #index for every subject/trial
    data[, (nm) := lapply(.SD, vec.norm, idx=idxDT[.I,idx], scale=scale), .SDcols=nm, by=.(subject,trial) ]
  }
  return(data[]) #points to same DT-object if no copying performed
}

data.subset_samples <- function(data, select, invert=F, update=T) {
  #select: samples to select or remove if invert = T
  #update: if TRUE, relabel samples to go from 1:N after subsetting
  data = data.check(data)
  if (length(select) == 1) select = 1:select
  data = data[ if (invert) !sample %in% select else sample %in% select ]
  if (update) data[, sample := seq_len(.N), by=.(subject,trial)] #count samples from 1:N
  return( data.check(data)[] )
}

data.subset_classes <- function(data, classes) {
  #classes: vector specifying the outcome values to subset
  data = data.check(data)
  if ( data[, is.factor(outcome)] ) classes = as.character(classes)
  #subset and drop factor levels if applicable
  return( data.check(data[ outcome %in% classes ][, outcome := if (is.factor(outcome)) ordered(outcome) else outcome ])[] )
}

data.merge_classes <- function(data, classes, new.labels=NULL, copy=F) {
  #classes: vector specifying the outcome values to merge or list of vectors for multiple merges
  #         e.g. c("A","B") to merge the outcomes A and B;
  #         e.g. list(c("A","B"),c("C","D")) to merge A and B, and also C and D
  #new.labels: optional vector specifying the new label(s) for the merged class(es)
  #            each position of the vector corresponds to the respective merge
  #            e.g. classes = c("A","B") and new.labels = "C" would merge A and B to C
  #            e.g. list(c("A","B"),c("C","D")) and new.labels = c("E","F") merges A and B to E; and C and D to F
  #copy: modify DT directly or create a copy. non-DT input will be copied necessarily
  if (copy) data = copy(data)
  data = data.check(data)
  if (!is.list(classes)) classes = list(classes)
  if ( data[, is.factor(outcome)] ) classes = lapply(classes, as.character)
  for (i in seq_along(classes)) {
    data[ outcome %in% classes[[i]], outcome := if (is.null(new.labels[i])) classes[[i]][1] else new.labels[i] ]
  }
  if ( data[, is.factor(outcome)] ) data[, outcome := ordered(outcome)] #drop levels
  return( data.check(data)[] )
}

data.remove_outliers <- function(data, threshold=NULL, Q=50, IQR=90, C=3, min.var=0, plot=F) {
  #detect outlier trials based on variance criterion, i.e. trial's whose averaged channel-wise variance
  #exceeds a threshold are removed from the data. threshold formula: Q+C*IQR
  #threshold: can be manually set and thus fixed for all subjects
  #Q represents a quantile of the trial variance values, e.g. Q=50=median
  #IQR represents a range between two quantiles of the variance, i.e. IQR=90=range between 5% and 95% quantiles
  #C is a constant that multiplies the inter-quantile range
  #min.var: minimum variance a trial has to exceed before it is considered as outlier
  #plot: if True, variance, thresholds and outliers are plotted per subject
  #returns data with a summary attribute which lists the number of outliers, thresholds and median variance per subject
  data = data.check(data)
  nm = setdiff( names(data), key(data) )
  #compute global variance values
  trialvar = data[, { chvar = sapply(.SD, var, na.rm=T); .(trialvar = mean(chvar)) }, .SDcols=nm, by=.(subject,trial)]
  if ( is.null(threshold) ) { #define threshold via automatic criterion
    if ( any( findInterval( c(Q,IQR), c(0,100) ) != 1 ) ) stop( "Q and IQR have to be in the range 0 to 100." )
    IQR = c( 0+(100-IQR)/200, 1-(100-IQR)/200 ) #convert IQR values to 0-1 range
    thresholds = trialvar[, { #calculate thresholds per subject based on global variance values
      lims = quantile(trialvar, probs=c(Q/100, IQR), na.rm=T);
      .(threshold = lims[1] + C * (lims[3]-lims[2])) #formula: Q + C * IQR
      }, by=subject]
  } else { #create a DT with fixed (user defined) threshold value for every subject
    thresholds = data[, .(threshold = threshold), by=subject]
  }
  #determine outlier trials which exceed the threshold
  outlier = trialvar[thresholds][, .(outlier = trialvar > max(threshold, min.var)), by=.(subject,trial)]
  cat( outlier[outlier==T, .N], "outlier trials detected.\n" )
  if (plot) { #visualize the outliers
    summary = trialvar[thresholds][outlier] #join the DTs
    for (s in thresholds[, subject]) { #iterate through subjects and plot the distribution, threshold and highlight outliers
      summary[subject == s, { 
        plot( trialvar, xlab="Trial", ylab="Variance", las=1, ylim = c(0,max(trialvar,threshold[1],min.var)) )
        abline(h=threshold[1], col="red", lty=2)
        if (min.var > 0) abline(h=min.var, col="blue", lty=2)
        text(trialvar, labels=ifelse(outlier,as.character(trial),""), cex=.7, pos=4)
        text(x=0, y=max(trialvar,threshold[1],min.var), labels=paste0("S",s))
        }]
    }
  }
  #subset non-outlier trials
  data = data[outlier][outlier==F]
  summary = outlier[outlier==T,.(outliers=.N),by=subject][thresholds][trialvar[, .(median.var=median(trialvar)), by=subject]]
  return( setattr( data.check(data[, outlier:=NULL]), "summary", summary[is.na(outliers), outliers:=0] )[] )
}

slice.idx <- function(n, window, overlap=0, imbalance="ignore") {
  #n: number of samples per trial or vector specifying min/max sample number
  #window: number of samples per slice that are contributed by each trial
  #overlap: >= 0 number of overlapping samples per slice, i.e. to have a sliding time window
  #         an overlap of 0 results in simple binning of the sample sequence where each slice
  #         contains new samples (i.e. none were part of the previous slice);
  #         any other number specifies how many of the next slice's samples were present in the previous one  
  #imabalance: a character specifying what happens when the slice parameters don't exactly match n;
  #            discard: remove samples, i.e. drop last slice, 
  #            adjust: increase overlap for the last slice to have same window size,
  #            ignore: last slice has fewer samples (i.e. smaller window) but same overlap
  #RETURNS: a matrix where column = slice number, rows = slice samples
  #NOTE: subset with data[ sample %in% slice.idx[,slicenumber] ]
  min.n = ifelse(length(n) > 1, min(n), 1)
  max.n = max(n)
  if ( window > max.n-min.n+1 ) stop( "Window must be smaller or equal the number of slices." )
  if ( overlap >= window ) stop( "Overlap must be smaller than window." )
  else if (overlap < 0) overlap = 0
  imbalance = tolower(imbalance)
  start = seq(min.n, 1+max.n-window, by=window-overlap)
  end = seq(min.n+window-1, max.n, by=window-overlap)
  idx = mapply(seq, start, end) 
  #ensure idx is a matrix and not a vector when window = 1
  if ( !is.matrix(idx) ) idx = matrix(idx, nrow=window, byrow=T)
  #check if parameters add up
  miss = max.n-end[length(end)]
  if (miss) { #not all samples were included
    if ( !imbalance %in% c("ignore","discard","adjust") ) stop( "Undefined method to handle imbalance." )
    if ( imbalance == "discard" ) { #do nothing (excludes the samples from the output)
      cat( "Last", miss, paste0("sample", ifelse(miss > 1, "s", "")), "are discarded.\n" )
    } else if ( imbalance == "ignore" ) { #add slice with smaller window and same overlap
      fix = (end[length(end)]-overlap+1):max.n
      idx = cbind(idx, c( fix, rep(NA, window-length(fix)) )) #add missing as NA
      cat( "Last slice is", window-length(fix), paste0("sample", ifelse(window-length(fix) > 1, "s", "")), "shorter.\n" )
    } else if ( imbalance == "adjust" ) { #increase window to keep overlap constant
      idx = cbind(idx, (max.n-window+1):max.n)
      o = sum( idx[, ncol(idx)-1] %in% idx[, ncol(idx)] ) - overlap
      cat( "Last slice's overlap is increased by", o, paste0("sample", ifelse(o > 1, "s", ""), "."),"\n" )
    }
  }
  return(idx) 
}

data.simulate <- function(n.subjects=1, n.outcomes=2, n.channels=60, n.trials=50, n.samples=100, SNR=.5,
                          sample.rate=50, amplitudes=NULL, frequencies=NULL, space=NULL, time=NULL) {
  ## create simulated data in typical EEG format
  #INPUT ---
  #n.subjects: number of times to repeat the data generation (simulating a subject each)
  #n.outcomes: number of outcomes, i.e. signals in the data
  #n.channels: number of measurement columns (e.g. EEG electrodes)
  #n.trials: number of trials per outcome
  #n.samples: number of samples per trial
  #SNR: signal to noise ratio, a value between 0 and 1
  #sample.rate: sampling rate (n.samples per second)
  #amplitudes: a vector specifying the signal amplitude for each outcome 
  #            if NULL, a random number between 1 and 5 per outcome
  #frequencies: a vector with 2 values per outcome specifying the frequency band
  #             if NULL, random band per outcome 
  #space: a vector with 2 values per outcome for the spatial distribution of the signals,
  #       i.e. the channels in which the signal of the outcome will be present
  #       if NULL, random channel range is chosen per outcome
  #time: a vector with 2 values per outcome for the on- and offset of the signals within a trial
  #RETURNS ---
  #data table with attribute "properties" of the mixed signals
  .eval_package("signal")
  if ( SNR < 0 || SNR > 1 ) stop( "The signal-to-noise ratio (SNR) takes a value between 0 and 1." )
  if ( any( !(c(2*length(amplitudes), length(frequencies), length(space), length(time))/n.outcomes) %in% c(0,2) ) ) {
    stop( "Signal properties need to be specified for each outcome individually." )
  }
  if (n.subjects>1) {
    args.in = lapply(as.list(match.call())[-1], eval); args.in$n.subjects=1
    data = replicate(n.subjects, {
      dat = do.call(data.simulate, args.in)
      dat[, subject:=NULL] 
    }, simplify=F)
    props = rbindlist(lapply(data, attr, "properties"), idcol="subject", use.names=T)
    data = rbindlist(data, idcol="subject", use.names=T)
    return( data.check( setattr(data, "properties", props) ) )
  }
  #set defaults: every outcome gets some attribute on the dimensions time, space, spectrum
  if ( is.null(amplitudes) ) { #randomply pick a sd
    amplitudes = sample(1:5, size=n.outcomes, replace=T)
  }
  if ( is.null(frequencies) ) { #randomly pick a spectrum
    frequencies = sample(2:round(sample.rate/4), size=n.outcomes, replace=T)
    frequencies = matrix( c(frequencies, frequencies + 1), ncol=n.outcomes, byrow=T )
  }
  frequencies = matrix(frequencies, ncol=n.outcomes) #1 column per outcome
  if ( is.null(space) ) { #spatial distribution of the signal
    space = matrix( c( sample(1:ceiling(n.channels/2), size=n.outcomes, replace=T),
                       sample(floor(n.channels/2+1):n.channels, size=n.outcomes, replace=T) ), ncol=n.outcomes, byrow=T )
  }
  space = matrix(space, ncol=n.outcomes) #1 column per outcome
  if (is.null(time) ) { #define on- and offset of signal
    time = matrix( c( sample(1:ceiling(n.samples/4), size=n.outcomes, replace=T),
                      sample(floor(3*n.samples/4+1):n.samples, size=n.outcomes, replace=T) ), ncol=n.outcomes, byrow=T )
  }
  time = matrix(time, ncol=n.outcomes) #1 column per outcome
  
  #create oscillatory signals from gaussian noise by applying a bandpass butter filter
  signals = sapply( seq_len(n.outcomes), function(s) {
    bt = signal::butter(n=5,W=c(frequencies[1,s]/sample.rate*2, frequencies[2,s]/sample.rate*2), type="pass")
    tmp = rnorm(n=n.samples, sd=amplitudes[s])
    x = signal::filtfilt(bt, tmp[ time[1,s]:time[2,s] ])
    c( rep(0, time[1,s]-1), x, rep(0, n.samples-time[2,s]) ) #append 0 before/after
  })
  if (!is.matrix(signals)) signals = matrix(signals, nrow=n.samples)
  #create smooth random source patterns
  patterns = matrix(rnorm(n.channels*n.outcomes), nrow=n.channels)
  patterns = apply( patterns, 2, signal::filtfilt, filt=signal::Ma(b=c(1,1)) ) #moving average filter
  #create prototype trials for each outcome
  trials = lapply( seq_len(n.outcomes), function(s) {
    X = signals[, s] %*% t(patterns[ space[1,s]:space[2,s], s])
    #append 0 columns to align dimensions between mixed signals
    cbind( matrix(0, nrow=n.samples, ncol=space[1,s]-1), X, 
           matrix(0, nrow=n.samples, ncol=n.channels-space[2,s]) )
  })
  #create all trials
  data = do.call(rbind, lapply(trials, function(trial) {
    #repeat the prototypical trial n times
    d = matrix( rep( t(trial), n.trials), ncol=n.channels, byrow=T )
    #create sensor noise
    noise = matrix(rnorm(n.trials*n.samples*n.channels), nrow=n.trials*n.samples, ncol=n.channels)
    #add them together
    SNR * d + (1-SNR) * noise
  }) )
  #generate outcome
  shuffleidx = sample(1:(n.outcomes*n.trials))
  outcome = rep(LETTERS[1:n.outcomes], each=n.trials)[shuffleidx]
  #split and bind in the order set by shuffleidx
  data = data.table::rbindlist( split( as.data.frame(data), rep(1:(n.outcomes*n.trials), each=n.samples) )[shuffleidx] )
  setnames( data, paste0("channel",seq_len(n.channels)) )
  data = data.setinfo(data, .outcome=outcome, .samples=n.samples, .trials=n.outcomes*n.trials)
  #add data generation information as attributes and return
  props = data.table(outcome = LETTERS[1:n.outcomes], amplitude=amplitudes, lower.freq=frequencies[1,], upper.freq=frequencies[2,],
                     channel.onset=space[1,], channel.offset=space[2,], time.onset=time[1,], time.offset=time[2,])
  return( setattr(data, "properties", props)[] )
}

parBackend <- function(cores=NULL, ...) {
  ## initialize or shut down a parallel backend
  #cores: numeric specifying the number of cores to initiate for parallel backend; NULL defaults to available cores;
  #       or a cluster object to shut down; or an empty list to take no action (to leave a cluster active)
  #...: pass the argument 'required' to specify number of jobs (internally used by functions)
  #RETURNS: the cluster object (list) or nothing (if shut down) or passes the empty list without action
  .eval_package("foreach", load=T)
  if ( !is.list(cores) ) { #initialize a backend
    start.time = proc.time()[3] #includes time spent on initialization of the cluster
    cores = min( min( cores, parallel::detectCores() ), list(...)$required ) #sanity check user input 
    #start up the CPU cluster
    cluster = parallel::makeCluster(cores) #create
    doParallel::registerDoParallel(cluster) #register
    #TEMP export functions to cluster
    parallel::clusterEvalQ( cluster, expr={ source("decode.R"); require(data.table) } )
    cat( "Distributing", list(...)$required, "jobs to", cores, "workers . . .\n" )
    return( setattr(cluster,"start.time",start.time) ) #return the cluster with time attribute
  } else if ( length(cores) > 0 && getDoParRegistered() ) { #active cluster to shut down
    parallel::stopCluster(cores)
    #remove cluster from foreach environment (otherwise throws error with no new backend)
    env = foreach:::.foreachGlobals 
    rm(list=ls(name=env), pos=env) #getDoParRegistered() evaluates to FALSE again
    invisible(gc()) #clean-up after computation
    if ( !is.null(attr(cores,"start.time")) ) {
      cat( "Done. Time elapsed (mins):", round((proc.time()[3]-attr(cores,"start.time"))/60,2), "\n" )
    }
  } else { #take no action (leave cluster active)
    if ( is.null(attr(cores,"start.time")) ) {
      cat( "Distributing", list(...)$required, "jobs to", foreach::getDoParWorkers(), "workers . . .\n" )
      return( setattr(cores,"start.time",proc.time()[3]) )
    } else {
      cat( "Done. Time elapsed (mins):", round((proc.time()[3]-attr(cores,"start.time"))/60,2), "\n" )
    }
  }
}


#### private functions ####

.eval_ellipsis <- function(funStr, ...) {
  ## helper to evaluate input to ellipsis (...) before handing it to a function
  args.in = list(...) #input
  args.out = list() #output
  for (FUN in funStr) {
    fun.args = formals(FUN) #function arguments  
    use.args = modifyList(fun.args, args.in, keep.null=T) #overwrite matching arguments 
    use.args = use.args[ names(use.args) %in% names(fun.args) ] #use only legal arguments
    args.out[[FUN]] = use.args
  }
  return(args.out)
}

.eval_package <- function(pkgStrings, load=F) {
  ## helper to check if package(s) are present
  for (pkgStr in pkgStrings) {
    if ( requireNamespace(pkgStr, quietly=T) ) { #check if pkg exists
    } else { #install pkg if not present
      install.packages(pkgStr)
    }
    if (load) { #attach pkg
      suppressMessages( require(pkgStr, quietly=T, character.only=T) )
    }
  }
}





