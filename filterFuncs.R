## FILTER functions

# data.resample <- function(data, old.srate, new.srate, nCores=NULL) {
#   ## resample data trial- and column-wise using signal's resample
#   ## uses zero padding so that artifacts won't be carried into the data
#   #INPUT ---
#   #data: continuous df or list of trials, slices, subjects
#   #old.srate: current sampling rate
#   #new.srate: desired sampling rate
#   #nCores: if data is a subject list, number of CPU cores to use for parallelization
#   #        if NULL, automatic selection; if 1, sequential execution
#   #        if an empty list, an externally registered cluster will be used
#   #RETURNS:
#   #data with new sampling rate
#   require(utils, quietly=T); .eval_package(c("MASS","signal"))
#   data = data.check(data, aslist=T, strip=F) #transform to list
#   if ( "slices" %in% attr(data, "type") ) {
#     data = data.split_trials( data, strip=F )
#     .transformed = T #convert to df at the end
#   } else if ( "subjects" %in% attr(data, "type") ) { #parallelize subjects
#     pcheck = .parallel_check(required=length(data), nCores=nCores)
#     outlen = ifelse(length(data) > 100, length(data), 100) #for .maxcombine argument of foreach
#     data = foreach(d=data, .combine=list, .multicombine=T, .maxcombine=outlen) %dopar%
#       data.resample(d, old.srate, new.srate)
#     .parallel_check(output=pcheck)
#     attr(data, "type") = "subjects"
#     return(data)
#   }
#   nsamples = data.samplenum(data[[1]])
#   frac = MASS::fractions(new.srate/old.srate)
#   #calculate zero pad length:
#   #note: by using fraction, rounding error accumulation is minimized
#   if ( old.srate %% new.srate != 0 ) { 
#     temp = strsplit( attr(frac, "fracs"), "/" )
#     p = as.numeric( temp[[1]][1] )
#     q = as.numeric( temp[[1]][2] )
#   } else { #not a fraction
#     p = new.srate/old.srate 
#     q = 1
#   }
#   nPad = ceiling( (max(p,q) * 10)/q ) * q #N=10 Matlab default
#   n = ceiling( (2*nPad + nsamples)*frac ) #new number of samples
#   unPad = nPad * frac
#   #obtain datacol info with higher sample number by rebinding temporarily
#   datacols = is.datacol( as.data.frame( data.table::rbindlist(data) ) )
#   cat("Resampling:\n")
#   progBar = txtProgressBar(style=3) #progress bar shown during resampling
#   progBarsteps = seq( 1/length(data), 1, length.out=length(data) )
#   #don't interpolate across boundaries:
#   resampled = lapply( seq_along(data), function(i) {
#     temp = data[[i]][, datacols ]
#     padstart = temp[rep(1, nPad), ] #repeat first row n times
#     padend = temp[rep(nrow(temp), nPad), ] #repeat last row n times
#     temp = rbind(padstart, temp, padend) #zero padded trial data
#     #use signal's resample instead of simple linear interpolation (still less refined than Matlab SPT's resample)
#     # temp = sapply( apply(temp, 2, approx, n=n), "[[", 2)
#     temp = apply(temp, 2, signal::resample, p=p, q=q)
#     setTxtProgressBar(progBar, progBarsteps[i]) #update progress bar
#     data.frame(samples=1:(n-2*unPad), outcome=data[[i]][1,2], temp[(unPad+1):(nrow(temp)-unPad),])
#   })
#   close(progBar)
#   return( data.check(resampled, aslist=F, transform=.transformed) ) #transform to df if needed
# }

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

filter.apply <- function(data, coefficients, nCores=list()) {
  ## function to apply the filter coefficients b from filter.coefficients to the data
  ## for FIR filter uses zero-phase one-pass filtering
  ## for IIR filter uses zero-phase two-pass filtering
  ## in either case filtering is done within trials, i.e. boundaries are respected
  #INPUT ---
  #coefficients: filter coefficients, e.g. from filter.coefficients, signal::fir1, signal::butter
  #              either of class 'Ma' (moving average) or class 'Arma' (autoregressive Ma)
  #              if numeric, converted to class 'Ma' (FIR)
  #nCores: number of CPU cores to use for parallelization
  #        if NULL, automatic selection; if 1, sequential execution
  #        if an empty list, an externally registered cluster will be used
  #Notes ---
  #only filters numeric non-integer data
  #RETURNS ---
  #filtered data
  .eval_package("signal")
  if ( identical(class(data), "numeric") ) data = data.setinfo(data.table(outcome="0", .x=data), .samples=as.numeric(length(data))) #vector to DT
  else data = data.check(data)
  nm = setdiff( names(data), key(data) )
  if ( !class(coefficients) %in% c("Ma", "Arma") ) {
    if ( is.numeric(coefficients) ) coefficients = signal::Ma(coefficients) #assuming FIR 
    else stop( "coefficients must be of class 'Ma' or 'Arma', cf. signal package." )
  }
  FIR = class(coefficients) == "Ma" #moving average, otherwise IIR
  #transform data
  dataMat = unname( as.matrix(data[, .SD, .SDcols=nm]) ) #transform data to matrix (measurements only)
  dataInfo = data[, .SD, .SDcols=!nm] #info columns DT with subset information
  subIDs = data[, unique(subject)]
  trialList = do.call(c, lapply(subIDs, function(sub) { #create trial list per subject
    idx = dataInfo[.(sub), which=T] #indices for subject
    idx = split(idx, dataInfo[idx, trial]) #split for subject trials
    lapply(idx, function(tr) dataMat[tr,]) #subset data
  }))
  groupDelay = (length(coefficients)-1)/2
  padStart = rep(1, groupDelay) #repeat first row n times
  padEnd = rep(data[, max(sample)], groupDelay) #repeat last row n times
  cl = parBackend(nCores, required=length(trialList)) #initialize parallel backend (if not already)
  
  filtdata = foreach(trial=trialList, .combine=function(...) do.call(rbind, list(...)), 
                     .multicombine=T, .maxcombine=length(trialList)+1) %dopar% 
            {
              if (FIR) { #Ma
                apply(trial, 2, function(s) {
                  sf = signal::filter( filt=coefficients, x=c(s[padStart], s, s[padEnd]) ) #padded signal  
                  sf[(2*groupDelay+1):length(sf)] #remove padded data  
                })
              } else { #Arma
                apply(trial, 2, function(s) {
                  signal::filtfilt(filt=coefficients, x=s) #two-pass: forward and reverse filtering
                })
              }
            }

  parBackend(cl) #close backend (if initialized by function) and print final time
  if (identical(nm, ".x")) return( data$.x ) #return vector
  return( data.check( setnames(cbind(dataInfo, as.data.table(filtdata)), names(data)) )[] )
}


plot.frequencies <- function(data, srate, columns=NULL, xlims=NULL, plot=T, 
                             ylab="Spectral Power Density (dB/Hz)", lwd=1.5) {
  ## computes the averaged frequency response for all measurement columns
  #INPUT ---
  #srate: sampling rate of the data in Hz
  #columns: columns (channels) to average over, defaults to all non-key columns
  #xlims: limit the frequency range to plot, default is the full range (0:srate/2)
  #plot: if True, the frequency response is plotted
  #ylab: y axis label
  #lwd: line width
  #RETURNS ---
  #a DT with the frequency (Hz) and averaged power (dB/Hz) per subject
  data = data.check(data)
  if (is.null(columns)) columns = setdiff( names(data), key(data) )
  if (is.null(xlims)) xlims = c(0, srate/2)
  freqtab = data[, {
    M = floor(.N/2) + 1
    .(power = 10*log10( #dB conversion
      rowMeans( #averaged over the channels
        sapply(.SD, function(x) abs(fft(x)[1:M])^2/(srate*.N)) #time normalized frequency spectrum
      )), 
      frequency = 0:(M-1)*(srate/.N)) }, .SDcols=columns, by=.(subject,trial)][, {
        .(power = mean(power)) }, keyby=.(subject,frequency)] #average the power per frequency bin
  if (plot) { #plot the freq response per subject
    freqtab[frequency %between% xlims, {
      plot(frequency, power, type="l", las=1,  col="#0072BD", xlim=xlims, 
           lwd=lwd, xlab="Frequency (Hz)", ylab=ylab)
    }, by=subject]
  }
  return(freqtab)
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


