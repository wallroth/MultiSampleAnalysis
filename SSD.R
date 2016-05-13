## Spatio-Spectral Decomposition

## More general terminology explained (cf. Haufe et al. 2014)
# brain signals are measured in forms of activation patterns at channels + noise
# A = activation patterns; s = signals (isolated responses to stimulus), e = noise, X = measurement/data
# x = A*s(n)' + e(n) = forward (generative) model; s and e have length N (samples)
# each column of A represents the activation pattern for a certain signal
# i.e. a column vector of A (a) represents the strength of the signal in each channel
# A therefore has dimensions: channels x signals
# X is then a combination of activation patterns of superimposed signals (+ noise),
# e.g. for 1 signal: A(1:64,1) * S(1:N,1)' = [1:64,1:N] + e(1:64,1:N)
# for 5 signals: A(1:64,1:5) * S(1:N,1:5)' = [1:64,1:N] + e(1:64,1:N)
# -> values in A reflect the signal strength at different channels

# backward modeling tries to extract signals from measurements (reverse direction)
# W'*x(n) = s(n) = backward (discriminative) model
# W = transformation matrix/weights/spatial filters with dimensions: channels x signals
# usually we can assume signals < channels
# by projecting x(n) onto a column of W, each column (w) extracts one latent signal
# e.g. X(1:N,1:64) * W(1:64,1:5) = [1:N,1:5]
# the resulting columns are components which maximize similarity to a target (signal)
# and are a more informative lower-dimensional representation of the overall data
# filters (W) should amplify the signal of interest and suppress signals of no interest
# however W cannot be directly interpreted as it is a function of signal + noise
# -> transformation to A required!
# if W is a square matrix (64x64), A = the transpose of the inverse of W
# if W is not square because signals < channels, W is no longer invertible
# X must have full rank (signals are independent) & we assume noise to be uncorrelated to signal
# A = Cx * W * Cs^-1
# Cx = covariance of data; Cs = covariance of signals
# Cs is calculated with W' * Cx * W and has dimensions signals x signals
# -> W[1:64,1:5]' * Cx[1:64,1:64] * W[1:64,1:5] = Cs[1:5,1:5]
# -> Cx[1:64,1:64] * W[1:64,1:5] * Cs[1:5,1:5]^-1 = A[1:64,1:5]
# note: given A, we can also compute W = Cx^-1 * A * Cs
# if signals/components are uncorrelated (orthogonal): A = Cx * W

## Steps in the SSD framework:
#1. narrow bandpass filtering of target frequency, e.g. freq = 10-12; width=2
# this yields matrix Xs (signal) but includes noise
# width should be small (1-2 Hz) for better linear estimation of SNR
#2. flanking frequencies are min-width and max+width -> 8-10 and 12-14 Hz
# bandpass filtering from min-width : max+width and subtracting Xs yields Xn (noise)
# to get Xn one first bandpass filters the signal (e.g. 8-14Hz) and then uses 
# bandstop filtering with +/- 1 Hz transitions with respect to the bandpass (e.g. 9-13 Hz bandstop)
#3. estimation of time-averaged covariance matrices for Xs and Xn
# Cs = Xs' * Xs / N
# Cn = Xn' * Xn / N
# = covariance matrix with biased normalization factor (N instead of N-1)
# X = matrix(c(1,3,-7,3,9,2,-5,4,6), ncol=3, byrow=T)
# X0 = scale(X, scale=F) #sweep(X, 2, colMeans(X), "-") #demeaned X
# t(X0) %*% X0 / nrow(X0) #covariance with normalization factor N
# cov.wt(X, method="ML") #alternative implementation, or cov(X,1) in matlab
#4. find spatial filters with high variance at target frequency and low variance in noise freq
# generalized eigenvalue decomposition: Cs*w = λCn*w
#5. obtain SSD components by projecting data onto the weights. Dimensionality can be reduced 
# by choosing only the first few components (ordered like in PCA)

SSD.coefficients <- function(frequencies, noise.width=2, signal.position="center", ...) {
  ## get SSD filter coefficients to be supplied to SSD.filter
  #INPUT ---
  #frequencies: 4/6 (or 2) values in this order:
  #   - Signal bandpass: low, high; e.g. 10-12Hz 
  #   - Noise bandpass: low, high, surrounding the signal, e.g. 8-14Hz 
  #   - Noise bandstop: low, high; between signal and noise bp, e.g. 9-13Hz
  #     e.g. frequencies = c(10,12,8,14,9,13)
  #   if only 2 values, they are assumed as signal bandpass specification
  #   if 4 values, no bandstop is performed but instead the signal is subtracted
  #   from the noise (no transitions but likely to result in spectral leakage)
  #   at default, 6 frequencies will be set
  #noise.width: defines the width of the unattenuated noise range
  #             is used to calculate the noise band if only signal band is 
  #             specified. Otherwise the value is ignored.
  #signal.position: center, left, right; specifies if the frequency spectrum is closed to both sides
  #                 of the signal or open on one side. The noise bands will either surround the siganal
  #                 (center, default) or be only on the closed side of the signal spectrum
  #                 i.e. if left, the signal spectrum is open to the left and a noise band will only be
  #                 subtracted on the (closed) right side, e.g. when filtering the delta band (0-4Hz)
  #transwidth, srate, wtype: see filter.coefficients
  #transwidth will default to 1, wtype to blackman
  #RETURNS ---
  #a list with the filter coefficients for the two/three filtering steps of the SSD
  #i.e. signal bandpass, noise bandpass, noise bandstop (optionally)
  #Note: at default settings, the noise edges will directly border on the signal edges
  #      to minimize spectral leakage. The borders are determined via the transition width
  signal.position = tolower(signal.position)
  if ( !signal.position %in% c("center","left","right") ) stop( "Undefined signal position. Options are 'center', 'left', 'right'." )
  args.in = lapply( as.list(match.call())[-1], eval )
  args.coeff = do.call(.eval_ellipsis, modifyList(args.in, list(funStr="filter.coefficients")))$filter.coefficients
  if ( is.null(args.in$wtype) ) args.coeff$wtype = "blackman" #default to blackman for smaller deviation
  if ( is.null(args.in$transwidth) ) args.coeff$transwidth = 1 #default to 1
  if ( signal.position == "center" ) { #two-sided
    #check frequency specifications
    if ( length(frequencies) == 2 ) { #only signal defined
      cat( "Only signal bandpass frequencies specified. Defining noise correspondingly:\n" )
      #frequency is at the center of the band, so transition width has to be considered
      #minimized spectral leakage with the edges of the transition zones next to each other:
      noise.bs = c( frequencies[1] - args.coeff$transwidth, frequencies[2] + args.coeff$transwidth )
      noise.bp = c( noise.bs[1] - noise.width, noise.bs[2] + noise.width )
      frequencies = c( frequencies, noise.bp, noise.bs )
      cat( "Noise bandpass from",frequencies[3],"-",frequencies[4],"Hz",
           "and noise bandstop from",frequencies[5],"-",frequencies[6],"Hz.\n" )
    }
    if ( !length(frequencies) %in% c(4,6) ) stop( "Incorrect frequency specification! Set either 2, 4 or 6 frequencies." )
    #get filter coefficients (b):
    SSDcoeffs = list( signal.pass = do.call( filter.coefficients, modifyList(args.coeff, list(frequencies=frequencies[1:2], ftype="pass")) ), 
                      noise.pass = do.call( filter.coefficients, modifyList(args.coeff, list(frequencies=frequencies[3:4], ftype="pass")) ) )
    if ( length(frequencies) == 6 ) {
      SSDcoeffs$noise.stop = do.call( filter.coefficients, modifyList(args.coeff, list(frequencies=frequencies[5:6], ftype="stop")) )
    }
  } else { #left/right one-sided
    #check frequency specifications
    if ( length(frequencies) == 1 ) { #only signal defined
      cat( "Only signal frequency specified. Defining noise correspondingly:\n" )
      #frequency is at the center of the band, so transition width has to be considered
      #minimized spectral leakage with the edges of the transition zones next to each other:
      if ( signal.position == "right" ) {
        ftypes = c("high","low")  #noise only on the left side of the signal
        noise.bs = frequencies - args.coeff$transwidth
        noise.bp = noise.bs - noise.width
      } else { #left
        ftypes = c("low","high") #noise only on the right side of the signal
        noise.bs = frequencies + args.coeff$transwidth
        noise.bp = noise.bs + noise.width
      }
      frequencies = c( frequencies, noise.bp, noise.bs )
      cat( paste0("Noise ",ftypes[1],"pass at ",frequencies[2],", ",ftypes[2],"pass at ",frequencies[3]," Hz."),"\n" )
    }
    if ( !length(frequencies) %in% c(2,3) ) stop( "Incorrect frequency specification! Set either 1, 2 or 3 frequencies." )
    #get filter coefficients (b):
    SSDcoeffs = list( signal.pass = do.call( filter.coefficients, modifyList(args.coeff, list(frequencies=frequencies[1], ftype=ftypes[1])) ), 
                      noise.pass = do.call( filter.coefficients, modifyList(args.coeff, list(frequencies=frequencies[2], ftype=ftypes[1])) ) )
    if ( length(frequencies) == 3 ) {
      SSDcoeffs$noise.stop = do.call( filter.coefficients, modifyList(args.coeff, list(frequencies=frequencies[3], ftype=ftypes[2])) )
    }
  }
  return(SSDcoeffs)
}

SSD.filter <- function(data, SSDcoeffs, nCores=NULL) {
  ## create Signal and Noise by filtering supplied data with coefficients b
  ## expects a list with 2/3 different filter coefficients to use with filter.apply
  #INPUT ---
  #data: continuous df or list of trials, slices, subjects
  #SSDcoeffs: filter coefficients, see SSD.coefficients
  #nCores: if data is a subject list, number of CPU cores to use for parallelization
  #        if NULL, automatic selection; if 1, sequential execution
  #        if an empty list, an externally registered cluster will be used
  #RETURNS ---
  #list of 2 containing signal and noise, respectively
  data = data.check(data)
  cl = parBackend(nCores, required=data[, uniqueN(subject)]) #initialize parallel backend (if not already)
  cat( "Bandpass filtering signal . . .\n" )
  SSDdata = list( signal = filter.apply(data, SSDcoeffs[[1]], nCores=list()) )
  cat( "Bandpass filtering noise . . .\n" )
  SSDdata$noise = filter.apply(data, SSDcoeffs[[2]], nCores=list())
  if (length(SSDcoeffs) > 2) { #additionally apply stopband:
    cat( "Bandstop filtering noise . . .\n" )
    SSDdata$noise = filter.apply(SSDdata$noise, SSDcoeffs[[3]], nCores=list()) 
  } else { #subtract signal from noise if no bandstop supplied:
    nm = setdiff( names(SSDdata$signal), key(SSDdata$signal) )
    SSDdata$noise[, (nm) := .SD-SSDdata$signal[, .SD, .SDcols=nm], .SDcols=nm][]
  }
  parBackend(cl) #close backend (if initialized by function) and print final time
  return(SSDdata)
}

SSD.apply <- function(SSDdata, average=T, shrinkage=F, q=.1, denoise=T, patterns=F, nCores=list()) {
  ## Spatio-Spectral decomposition (SSD), cf. Nikulin et al. 2011 for  extraction of neuronal oscilations with improved SNR
  #INPUT ---
  #SSDdata: list of 2 (1: signal, 2: noise), output of SSD.filter
  #average: if True, compute the covariance within trials and average afterwards
  #         otherwise compute the covariance across all trials
  #shrinkage: if True, uses shrinkage to estimate the covariance matrices
  #           else the population covariance is computed after centering
  #q: control the number of selected components; 
  #   if < 1, the auto-selection criterion from Haufe et al., 2014 is applied:
  #   all components whose lambda is >= q*IQR(lambda) + 75thQ(lambda)
  #   if >= 1, a fixed number of components for every participant
  #denoise: if True, use low-rank factorization to denoise measurements (X*W*A')
  #         necessary if subsequent analyses should be performed in original input space
  #nCores: if data is a subject list, number of CPU cores to use for parallelization
  #        if NULL, automatic selection; if 1, sequential execution
  #        if an empty list, an externally registered cluster will be used
  #RETURNS ---
  #data in component space (if denoise = FALSE) or in original space with reduced rank (denoise = TRUE)
  if ( class(SSDdata) != "list" || length(SSDdata) != 2 ) {
    stop( "SSDdata must be a list with 2 elements: signal and noise. Refer to the documentation or use SSD.filter." )
  }
  if (shrinkage) .eval_package("corpcor")
  SSDdata = list(signal = data.check(SSDdata[[1]]), noise = data.check(SSDdata[[2]]))
  nm = setdiff( names(SSDdata$signal), key(SSDdata$signal) )
  Xs = as.matrix( SSDdata$signal[, .SD, .SDcols=nm] )
  Xn = as.matrix( SSDdata$noise[, .SD, .SDcols=nm] )
  norm.cov <- function(X) return( cov.wt(scale(X, scale=F), method="ML")$cov ) #population covariance (N) on centered data
  get.weights <- function(sub) {
    if (average) { #compute trial-averaged covariance matrices
      idx = split(1:nrow(sub$Xs), sub$trials) #split subject indices for trials
      if (shrinkage) {
        Cs = Reduce("+", lapply(idx, function(i) corpcor::cov.shrink(sub$Xs[i,], verbose=F)))/length(idx)
        Cn = Reduce("+", lapply(idx, function(i) corpcor::cov.shrink(sub$Xn[i,], verbose=F)))/length(idx)
      } else {
        Cs = Reduce("+", lapply(idx, function(i) norm.cov(sub$Xs[i,])))/length(idx)
        Cn = Reduce("+", lapply(idx, function(i) norm.cov(sub$Xn[i,])))/length(idx)
      }
    } else { #compute across trials
      if (shrinkage) {
        Cs = corpcor::cov.shrink(sub$Xs, verbose=F)
        Cn = corpcor::cov.shrink(sub$Xn, verbose=F)
      } else {
        Cs = norm.cov(sub$Xs)
        Cn = norm.cov(sub$Xn)
      }
    }
    #compute filters for whitening
    VD = eigen(Cs+Cn, symmetric=T)
    r = sum(VD$values > 10^-6*VD$values[1]) #estimate rank
    if (r < length(nm)) warning( "Not all columns are linearly independent. Computing only ",r," components." )
    M = diag(VD$values[1:r]^-0.5, nrow=r, ncol=r) %*% t(VD$vectors[,1:r])
    Cs_white = M %*% Cs %*% t(M) #whitened C signal
    WD = eigen(Cs_white, symmetric=T)
    WD$vectors = t(M) %*% WD$vectors #unwhiten the filters
    if (denoise || patterns) WD$C = Cs #add C for pattern matrix
    return(WD)
  }
  
  subIDs = SSDdata$signal[, unique(subject)]
  sublist = lapply(subIDs, function(sub) {
    isub = SSDdata$signal[.(sub), which=T] #subject indices
    list(Xs = Xs[isub,], Xn = Xn[isub,], trials=if(average)SSDdata$signal$trial[isub])
  })
  cl = parBackend(nCores, required=length(subIDs)) #initialize parallel backend (if not already)
  
  #compute the weights (eigenvectors) and lambdas (eigenvalues)
  WD = foreach(sub=sublist, .combine=list, .multicombine=T, .maxcombine=length(subIDs)+1) %dopar% get.weights(sub)
  if (length(subIDs)<2) WD = list(WD)
  r = q #fixed num components for all subjects if a constant >= 1
  rank = c() #save ranks
  if (patterns) pattern_list = list() #save patterns
  #project the data outside foreach to print feedback in case of auto-selection
  data = rbindlist( lapply(seq_along(subIDs), function(i) {
    W = WD[[i]]$vectors
    if (q < 1) { #auto-select: q*IQR(lambda) + Q75(lambda), but at least 1 component
      lambda = WD[[i]]$values
      r = max(1, sum( lambda >= q*IQR(lambda)+quantile(lambda, probs=.75) ))
      cat( "Subject ", subIDs[i], ": auto-selected ", r, " components.\n", sep="" )
    }
    X = sublist[[i]]$Xs %*% W[,1:r] #project data onto filters
    if (denoise || patterns) { #low-rank factorization
      rank <<- c(rank, r)
      A = WD[[i]]$C %*% W %*% solve(t(W) %*% WD[[i]]$C %*% W) #pattern matrix
      if (patterns) pattern_list[[ as.character(subIDs[i]) ]] <<- A[,1:r]
      if (denoise) X = X %*% t(A[,1:r]) #project X back to measurement space (with reduced rank)
    } #else: X remains in component space
    cbind(SSDdata$signal[.(subIDs[i]), .SD, .SDcols=!nm], X)
  }), fill=T ) #fill if auto-selection yielded different num components and no denoising is performed
  
  parBackend(cl) #close backend (if initialized by function) and print final time
  if (denoise) setattr(data, "ranks", rank)[]
  if (!patterns) return( data.check(data) )
  return( list(data=data.check(data), patterns=pattern_list) )
}

