#make helper functions accessible
source("myFuncs.R")

filter.get_coeffs <- function(frequencies, transwidth, ftype, srate=500, wtype="hamming") {
  ## creates Hamming window sinc FIR filter coefficients for the given frequency
  #INPUT ---
  #frequencies: scalar (high/lowpass) or tuple (bandpass/stop) to specifiy freq (range)
  #transwidth: scalar, specifies size of roll-off, freq is at the center of this value
  #         - e.g. 40Hz low-pass: freq 44, transwidth 8, roll-off between 40 and 48
  #ftype: "low", "high", "pass", "stop" filter type
  #wtype: hamming or blackman window function
  #RETURNS ---
  #b: filter coefficients for use with filter.apply or signal::filter in general
  #Notes ---
  #it is recommended to design bandpass as low-pass and high-pass filters separately (Widmann 2015)
  #-> high-pass filtering requires steep transition width whereas low-pass can be shallower
  #transition width: narrower = stricter (more precise but danger of ringing artifacts)
  #generally between 10%-25% of frequency bounds (but larger for very small frequencies)
  #cutoff frequencies are the center of the transition band
  #ripples in the respective band are the allowed deviation from 0 (stop) / 1 (pass)
  
  #windowing smoothes filter kernel to minimize edge artifacts
  #zero-padding replaces necessity of two-pass filtering (filtfilt)
  #upper frequency edge cannot be higher than srate/2 (nyquist frequency)

  #normalized transition widths from Widmann (2015) table p.8
  deltaF = ifelse( tolower(wtype) == "blackman", 5.5, 3.3 ) #blackman: 5.5, hamming: 3.3
  m = ceiling(deltaF/(transwidth/srate)/2)*2 #filter order
  #transform frequencies to range [0,1] where 1 corresponds to nyquist freq (srate/2):
  f = frequencies/(srate/2) 
  
  #fir1 sets transition zones to 0 and smoothes afterwards with windowing
  # b = fir1(m, f, type=ftype, window=hamming(m+1), scale=T)
  #using firws import (Widman) instead of fir1
  b = .firfilt.firws(m, f, type=ftype, window=wtype)
  return(b)
}

filter.apply <- function(data, b) {
  ## function to apply the filter coefficients b from filter.get_coeffs to the data
  ## uses zero-phase one-pass filtering approach within trials 
  #INPUT ---
  #data: df, matrix or list
  #b: filter coefficients, e.g. from filter.get_coeffs
  #Notes ---
  #tries to not filter across boundaries (trials)
  #only filters numeric non-integer data
  require(utils, quietly=T)
  .eval_package("signal")
  if ( class(data) == "numeric" ) {
    data = list( as.data.frame(data) ) #vector to df
    .transformed = TRUE
  }
  trialdata = data.check(data, aslist=T, strip=F) #transform to list
  infoidx = !is.datacol(trialdata[[1]]) #non-measurement data
  colnames = names(trialdata[[1]]) #save variable names
  print("Filtering:")
  progBar = txtProgressBar(style=3) #progress bar shown during filtering
  progBarsteps = seq( 1/length(trialdata), 1, length.out=length(trialdata) )
  #iterate over trials:
  filtdata = lapply(1:length(trialdata), function(i) {
    data = trialdata[[i]]
    infocols = data[,infoidx] #info
    data = as.matrix(data[,!infoidx]) #measurements
    #zero-padding to circumvent delay of filter when one-pass filtering:
    groupdelay = (length(b)-1)/2
    padstart = data[rep(1, groupdelay), ] #repeat first row n times
    padend = data[rep(nrow(data), groupdelay), ] #repeat last row n times
    if (dim(data)[2] == 1) { #signal is only a vector
      paddedsig = as.matrix(c(padstart, data, padend)) #zero padded signal
    } else { #signal is a matrix with >= 2 cols
      paddedsig = rbind(padstart, data, padend) #zero padded signal
    }    
    #apply filter:
    temp = signal::filter(signal::Ma(b), paddedsig) #time-series with dims: [ 1, nrow*ncol ]
    temp = matrix(temp, ncol=ncol(data)) #retransform into matrix
    fsignal = temp[(2*groupdelay+1):nrow(temp),] #remove padded data
    setTxtProgressBar(progBar, progBarsteps[i]) #update progress bar
    setNames( cbind(infocols, fsignal), colnames )
  })
  close(progBar)
  return( data.check(filtdata, aslist=F, transform=.transformed) )
}


data.freq_response <- function(data, cols=NULL, srate=500, 
                               xlims=NULL, filt=F, title="") {
  ## plots the averaged frequency response for all measurement columns
  #INPUT ---
  #data: df or list of trial data
  #cols: columns (channels) to average over, defaults to all that contain non-integer numeric values
  #srate: sampling rate of the data in Hz
  #xlims: limit the frequency range to plot, default is c(1, srate/2)
  #filt: if True, the averaged freq response is low-pass filtered at 200Hz for visualization
  #title: optional title for the plot
  
  data = data.check(data, aslist=F) #transform to df
  if (is.null(cols)) { #defaults to all channels if nothing specified
    data = data[, is.datacol(data)] #analzye only measurements 
    cols = 1:ncol(data)
  } else {
      cols = cols[ is.datacol(data[,cols]) ] #analzye only measurements 
  }
  if (is.null(xlims)) { 
    xlims = c(1, srate/2) 
  }
  M = nrow(data/2 + 1)
  temp = sapply(cols, function(c) {
    X=abs(fft(data[,c])[1:M])^2 / (srate*nrow(data)) #time normalized
  })
  xAxis = 0:(M-1) * (srate/nrow(data))
  temp = 10*log10(temp) #dB conversion
  #note: 20log10(V) is equivalent to 10log10(V^2) [V=Voltage magnitude]
  avgspec = rowMeans(temp) #averaged over all channels
  if (filt) { #apply 200Hz lowpass for visualization:
    b = filter.get_coeffs(225,50,"low",nrow(data))
    avgspec = filter.apply(avgspec, b)[,1] #unroll matrix
    #note: don't filter if you want to see line noise
  }
  avgspec = avgspec[xAxis >= xlims[1] & xAxis <= xlims[2]]
  xAxis = xAxis[xAxis >= xlims[1] & xAxis <= xlims[2]]
  plot(xAxis, avgspec, type="l", las=1,
       col=rgb(0,0.4470,0.7410), xlim=xlims, lwd=1.5,
       xlab="Frequency (Hz)", ylab="Spectral power density (dB/Hz)", main=title)
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


