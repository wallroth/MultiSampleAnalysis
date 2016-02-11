## FILTER functions

data.resample <- function(data, old.srate, new.srate, nCores=NULL) {
  ## resample data trial- and column-wise using signal's resample
  ## uses zero padding so that artifacts won't be carried into the data
  #INPUT ---
  #data: continuous df or list of trials, slices, subjects
  #old.srate: current sampling rate
  #new.srate: desired sampling rate
  #nCores: if data is a subject list, number of CPU cores to use for parallelization
  #        if NULL, automatic selection; if 1, sequential execution
  #        if an empty list, an externally registered cluster will be used
  #RETURNS:
  #data with new sampling rate
  require(utils, quietly=T); .eval_package(c("MASS","signal"))
  data = data.check(data, aslist=T, strip=F) #transform to list
  if ( "slices" %in% attr(data, "type") ) {
    data = data.split_trials( data, strip=F )
    .transformed = T #convert to df at the end
  } else if ( "subjects" %in% attr(data, "type") ) { #parallelize subjects
    pcheck = .parallel_check(required=length(data), nCores=nCores)
    data = foreach(d=data, .combine=list, .multicombine=T) %dopar%
      data.resample(d, old.srate, new.srate)
    .parallel_check(output=pcheck)
    attr(data, "type") = "subjects"
    return(data)
  }
  nsamples = data.samplenum(data[[1]])
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
  nPad = ceiling( (max(p,q) * 10)/q ) * q #N=10 Matlab default
  n = ceiling( (2*nPad + nsamples)*frac ) #new number of samples
  unPad = nPad * frac
  #obtain datacol info with higher sample number by rebinding temporarily
  datacols = is.datacol( as.data.frame( data.table::rbindlist(data) ) )
  cat("Resampling:\n")
  progBar = txtProgressBar(style=3) #progress bar shown during resampling
  progBarsteps = seq( 1/length(data), 1, length.out=length(data) )
  #don't interpolate across boundaries:
  resampled = lapply( seq_along(data), function(i) {
    temp = data[[i]][, datacols ]
    padstart = temp[rep(1, nPad), ] #repeat first row n times
    padend = temp[rep(nrow(temp), nPad), ] #repeat last row n times
    temp = rbind(padstart, temp, padend) #zero padded trial data
    #use signal's resample instead of simple linear interpolation (still less refined than Matlab SPT's resample)
    # temp = sapply( apply(temp, 2, approx, n=n), "[[", 2)
    temp = apply(temp, 2, signal::resample, p=p, q=q)
    setTxtProgressBar(progBar, progBarsteps[i]) #update progress bar
    data.frame(samples=1:(n-2*unPad), outcome=data[[i]][1,2], temp[(unPad+1):(nrow(temp)-unPad),])
  })
  close(progBar)
  return( data.check(resampled, aslist=F, transform=.transformed) ) #transform to df if needed
}

filter.coefficients <- function(frequencies, transwidth, ftype, srate, wtype="hamming") {
  ## creates windowed sinc FIR filter coefficients for the given frequency
  #INPUT ---
  #frequencies: scalar (high/lowpass) or tuple (bandpass/stop) to specify freq (range)
  #transwidth: scalar, specifies size of roll-off, freq is at the center of this value
  #            - e.g. 40Hz low-pass: freq 44, transwidth 8, roll-off between 40 and 48
  #ftype: "low", "high", "pass", "stop" filter type
  #wtype: hamming or blackman window function
  #RETURNS ---
  #b: filter coefficients for use with filter.apply or signal::filter in general
  #Notes ---
  #it is recommended to design bandpass as low-pass and high-pass filters separately (Widmann 2015)
  #-> high-pass filtering requires steep transition width whereas low-pass can be shallower
  #transition width: narrower = stricter (more precise but danger of ringing artifacts)
  #generally between 10%-25% of frequency bounds (but larger for very small frequencies)
  #ripples in the respective band are the allowed deviation from 0 (stop) / 1 (pass)
  #windowing smoothes filter kernel to minimize edge artifacts
  #zero-padding replaces necessity of two-pass filtering (filtfilt)
  if ( max(frequencies) > srate/2 ) stop( "Upper frequency limit is the Nyquist frequency corresponding to srate/2." )
  .eval_package("signal")
  #normalized transition widths from Widmann (2015) table p.8
  deltaF = ifelse( tolower(wtype) == "blackman", 5.5, 3.3 ) #blackman: 5.5, hamming: 3.3
  m = ceiling(deltaF/(transwidth/srate)/2)*2 #filter order
  #transform frequencies to range [0,1] where 1 corresponds to nyquist freq (srate/2):
  f = frequencies/(srate/2) 
  #using firws import (Widman) instead of fir1
  b = .firfilt.firws(m, f, type=ftype, window=wtype)
  return( signal::Ma(b) )
}

filter.apply <- function(data, coefficients, nCores=NULL) {
  ## function to apply the filter coefficients b from filter.coefficients to the data
  ## for FIR filter uses zero-phase one-pass filtering
  ## for IIR filter uses zero-phase two-pass filtering
  ## in either case filtering is done within trials, i.e. boundaries are respected
  #INPUT ---
  #data: df, matrix, vector or list of trials, slices, subjects
  #coefficients: filter coefficients, e.g. from filter.coefficients, signal::fir1, signal::butter
  #              either of class 'Ma' (moving average) or class 'Arma' (autoregressive Ma)
  #              if numeric, converted to class 'Ma' (FIR)
  #nCores: if data is a subject list, number of CPU cores to use for parallelization
  #        if NULL, automatic selection; if 1, sequential execution
  #        if an empty list, an externally registered cluster will be used
  #Notes ---
  #only filters numeric non-integer data
  #RETURNS ---
  #filtered data
  require(utils, quietly=T);  .eval_package("signal")
  if ( class(data) == "numeric" ) { #vector to df in a list
    data = list( as.data.frame(data) ) 
    attr(data, "type") = "trials" #camouflaged as trial list
    attr(data, "toVec") = T #retransform later
  }
  data = data.check(data, aslist=T, strip=F) #transform to list
  if ( "slices" %in% attr(data, "type") ) {
    data = data.split_trials( data, strip=F )
    .transformed = T #convert to df at the end
  } else if ( "subjects" %in% attr(data, "type") ) { #parallelize subjects
    pcheck = .parallel_check(required=length(data), nCores=nCores)
    data = foreach(d=data, .combine=list, .multicombine=T) %dopar%
      filter.apply(d, coefficients)
    .parallel_check(output=pcheck)
    attr(data, "type") = "subjects"
    return(data)
  }
  #obtain infocols with higher sample number by rebinding temporarily
  infoidx = !is.datacol( as.data.frame( data.table::rbindlist(data) ) )
  colnames = names(data[[1]]) #save variable names
  if ( class(coefficients) != "Ma" && is.numeric(coefficients) ) {
    coefficients = signal::Ma(coefficients) #assuming FIR 
  }
  cat("Filtering:\n")
  progBar = txtProgressBar(style=3) #progress bar shown during filtering
  progBarsteps = seq( 1/length(data), 1, length.out=length(data) )
  #iterate over trials:
  filtdata = lapply( seq_along(data), function(i) {
    trial = data[[i]]
    infocols = trial[, infoidx] #info
    trial = as.matrix(trial[, !infoidx]) #measurements
    if ( class(coefficients) == "Ma" ) { #FIR
      #zero-padding to circumvent delay of filter when one-pass filtering:
      groupdelay = (length(coefficients)-1)/2
      padstart = trial[rep(1, groupdelay), ] #repeat first row n times
      padend = trial[rep(nrow(trial), groupdelay), ] #repeat last row n times
      if (dim(trial)[2] == 1) { #signal is only a vector
        paddedsig = as.matrix(c(padstart, trial, padend)) #zero padded signal
      } else { #signal is a matrix with >= 2 cols
        paddedsig = rbind(padstart, trial, padend) #zero padded signal
      }    
      #apply filter:
      temp = apply(paddedsig, 2, signal::filter, filt=coefficients)
      fsignal = temp[(2*groupdelay+1):nrow(temp),] #remove padded data
    } else if ( class(coefficients) == "Arma" ) { #IIR
      #two-pass filtering forward/backward to avoid phase-delay
      fsignal = apply(trial, 2, signal::filtfilt, filt=coefficients)
    } else {
      stop( "coefficients must be of class 'Ma' or 'Arma', cf. signal package." )
    }
    setTxtProgressBar(progBar, progBarsteps[i]) #update progress bar
    setNames( cbind(infocols, fsignal), colnames )
  })
  close(progBar)
  if ( !is.null(attr(data, "toVec")) ) return( filtdata[[1]][,1] ) #return as vector
  return( data.check(filtdata, aslist=F, transform=.transformed) )
}


plot.frequencies <- function(data, srate, boundaries=F, cols=NULL, xlims=NULL, plot=T, 
                             title="", ylab="Spectral Power Density (dB/Hz)", lwd=1.5) {
  ## computes the averaged frequency response for all measurement columns
  ## values are the spectral power density (db/Hz)
  #INPUT ---
  #data: df or list of trials, slices
  #cols: columns (channels) to average over, defaults to all that contain non-integer numeric values
  #srate: sampling rate of the data in Hz
  #boundaries: should the frequency computation respect boundaries (trials) or not
  #xlims: limit the frequency range to plot, default is the full range (0:srate/2)
  #filt: if True, the averaged freq response is low-pass filtered at 200Hz for visualization
  #plot: if True, the frequency response is plotted. Else, the values are returned
  #title: optional title for the plot
  #ylab: y axis label
  #RETURNS ---
  #either nothing and the function directly plots (if plot=T) or the computed 
  #spectral power density with names corresponding to the frequency position
  data = data.check(data, aslist=F) #transform to df
  measurements = is.datacol(data)
  if ( is.null(cols) ) { #defaults to all channels if nothing specified
    cols = which(measurements)
  } else {
    cols = cols[ measurements[cols] ] #analyze only measurements 
  }
  if ( is.null(xlims) ) xlims = c(0, srate/2)
  if (!boundaries) { #across the whole data
    M = nrow(data)/2 + 1
    temp = sapply(cols, function(c) {
      X=abs(fft(data[,c])[1:M])^2 / (srate*nrow(data)) #time normalized
    })
    xAxis = 0:(M-1) * (srate/nrow(data))
  } else { #trial-wise computation
    nsamples = data.samplenum(data)
    M = floor(nsamples/2) + 1
    data = data.split_trials(data)
    temp = Reduce("+", lapply(data, function(trial) {
      sapply(cols, function(c) {
        X=abs(fft(trial[,c])[1:M])^2 / (srate*nsamples) #time normalized
      })
    }) )/length(data)
    xAxis = 0:(M-1) * (srate/nsamples)
  }
  temp = 10*log10(temp) #dB conversion
  #note: 20log10(V) is equivalent to 10log10(V^2) [V=Voltage magnitude]
  avgspec = rowMeans(temp) #averaged over all channels
  avgspec = avgspec[xAxis >= xlims[1] & xAxis <= xlims[2]]
  xAxis = xAxis[xAxis >= xlims[1] & xAxis <= xlims[2]]
  if (plot) {
    plot(xAxis[ xAxis >= xlims[1] & xAxis <= xlims[2] ], 
         avgspec[ xAxis >= xlims[1] & xAxis <= xlims[2] ], 
         type="l", las=1, col="#0072BD", xlim=xlims, lwd=lwd,
         xlab="Frequency (Hz)", ylab=ylab, main=title)
  } else { #return the values
    return( setNames(avgspec, xAxis) )
  }
}




############# copied firfilt (Widmann) Matlab functions #############

#Widman recommends 0.002-0.001 passband ripple (-27 - -30dB) & -54 to -60 dB stopband attenuation
.firfilt.firws <- function(m, freq, type, window) {
  ## firws creates the filter coefficients b
  #f as nyquist freq: range=[0 1]
  #type options: "high", "low", "pass", "stop"
  #window: automatically set to hamming, optionally blackman
  firfilt.windows <- function(m, window) {
    #compute window coefficients, copied from firfilt/windows.m
    m = seq(0, 1, by=1/(m-1))
    if ( tolower(window) == "blackman" ) {
      a = c(.42, .5, .08, 0)
      w = a[1] - a[2] * cos(2 * pi * m) + a[3] * cos(4 * pi * m) - a[4] * cos(6 * pi * m)
    } else { #hamming 
      a = 25/46
      w = a - (1 - a) * cos(2 * pi * m)
    }
    return(w)
  }
  
  firfilt.fkernel <- function(m, f, w) {
    #computing filter kernel, copied from firfilt/firws.m
    m = (-m/2) : (m/2)
    b = rep(0, length(m))
    b[m==0] = 2*pi*f
    b[m!=0] = sin(2*pi*f*m[m!=0]) / m[m!=0]
    b = b * w
    return(b / sum(b))
  }
  
  firfilt.fspecinv <- function(b) {
    #spectral inversion of b, copied from firfilt/firws.m
    b = -b
    b[(length(b)-1)/2+1] = b[(length(b)-1)/2+1] + 1
    return(b)
  }
  
  #actual firws function:
  w = firfilt.windows(m+1, window)
  f = freq/2 
  b = firfilt.fkernel(m, f[1], w)
  if (length(f) == 1 & type == "high") {
    b = firfilt.fspecinv(b)
  }
  if (length(f) == 2) {
    b = b + firfilt.fspecinv(firfilt.fkernel(m, f[2], w))
    if (type == "pass") {
      b = firfilt.fspecinv(b)
    }
  }
  return(b)
}
###


