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


### Steps of analysis in the SSD/SPoC framework:
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
#6. apply SPoC to SSD components

pipe.SSD <- function(data, ..., plotfreq = T, lowrank = NULL) {
  ## Pipeline for the whole SSD procedure
  #INPUT ---
  #data: continuous df or list of trials
  #frequencies, transwidth, srate: see SSD.coefficients
  #plotfreq: if True, frequency response of Signal/Noise is plotted
  #lowrank: int, if specified, SSD output is used to denoise measurements by low-rank
  #         factorization (measurements are then in original sensor space unlike components!)
  #RETURNS ---
  #same as SSD; if lowrank was specified additionally the denoised data
  
  args.in = .eval_ellipsis("SSD.coefficients", ...)
  #get the filters:
  SSDfilts = do.call(SSD.coefficients, args.in)
  #get signal and noise:
  Xs_Xn = SSD.filter(data, SSDfilts)
  if (plotfreq) { #plot frequency response for Signal and Noise
    freqs = args.in$frequencies[1:2]
    args.in = .eval_ellipsis("data.freq_response", ...)
    par(mfrow=c(2,1))
    do.call(data.freq_response, modifyList( args.in, list(data=Xs_Xn$signal, filt=T, title="Signal") ))
    abline(v = freqs, lty=2, col="red")
    do.call(data.freq_response, modifyList( args.in, list(data=Xs_Xn$noise, filt=T, title="Noise") ))
    abline(v = freqs, lty=2, col="red")
    par(mfrow=c(1,1))
  }
  #do the SSD:
  SSDout = SSD.apply(Xs_Xn)
  if ( !is.null(lowrank) ) {
    SSDout$X_denoised = SSD.denoise(SSDout$X_SSD, SSDout$A_SSD, lowrank)
  }
  #return SSD output
  return(SSDout)
}

SSD.coefficients <- function(frequencies, transwidth=1, srate, window="blackman") {
  ## get SSD filter coefficients to be supplied to SSD.filter
  #INPUT ---
  #frequencies: 6 (or 2) values in this order:
  #   - Signal bandpass: low, high; e.g. 10-12Hz 
  #   - Noise bandpass: low, high, surrounding the signal, e.g. 8-14Hz 
  #   - Noise bandstop: low, high; between signal and noise bp, e.g. 9-13Hz
  #     e.g. frequencies = c(10,12,8,14,9,13)
  #   if only 2 values, they are assumed as signal bandpass specification
  #transwidth: transition width in Hz for the filters to control roll-off
  #srate: sampling rate of the data in Hz
  #window: the windowing sinc function, default to the more restrictive blackman
  #RETURNS ---
  #a list with the filter coefficients for the three filtering steps of the SSD
  #i.e. signal bandpass, noise bandpass, noise bandstop
  #Note: Noise defaults w.r.t. Signal to +/- 3 Hz bandpass and +/- 1 Hz bandstop
  
  #check frequency specifications
  if ( length(frequencies) == 2 ) {
    cat( "Only signal bandpass frequencies specified. Assigning the others automatically:\n" )
    frequencies[3] = frequencies[1] - 3
    frequencies[4] = frequencies[2] + 3
    cat( "Setting noise bandpass to",frequencies[3],"-",frequencies[4],"Hz.\n" )
    frequencies[5] = frequencies[1] - 1
    frequencies[6] = frequencies[2] + 1
    cat( "Setting noise bandstop to",frequencies[5],"-",frequencies[6],"Hz.\n" )
  }
  if ( length(frequencies) != 6 ) { stop( "Incorrect frequency specification! ",
                                          "Set either 2 or 6 frequencies." ) }
  #get filter coefficients (b):
  Xs_bandpass = filter.coefficients( frequencies[1:2], transwidth, "pass", srate, window )
  Xn_bandpass = filter.coefficients( frequencies[3:4], transwidth, "pass", srate, window )
  Xn_bandstop = filter.coefficients( frequencies[5:6], transwidth, "stop", srate, window )
  return( list(Xs_bandpass=Xs_bandpass, 
               Xn_bandpass=Xn_bandpass, 
               Xn_bandstop=Xn_bandstop) )
}

SSD.filter <- function(data, SSDcoeffs, nCores=NULL) {
  ## create Signal and Noise by filtering supplied data with coefficients b
  ## expects a list with 3 different filter coefficients to use with filter.apply
  #INPUT ---
  #data: continuous df or list of trials, slices, subjects
  #SSDcoeffs: filter coefficients, see SSD.coefficients
  #nCores: if data is a subject list, number of CPU cores to use for parallelization
  #        if NULL, automatic selection; if 1, sequential execution
  #        can also be an already registered cluster object
  #RETURNS ---
  #list of 2 with dfs containing signal and noise, respectively
  data = data.set_type(data)
  if ( "subjects" %in% attr(data, "type") ) { #parallelize subjects
    pcheck = .parallel_check(required=length(data), nCores=nCores)
    SSDdata = foreach(d=data, .combine=list, .multicombine=T) %dopar%
      SSD.filter(d, SSDcoeffs)
    .parallel_check(output=pcheck)
    SSDdata = setNames( SSDdata, paste0("subject", seq_along(data)) )
    attr(SSDdata, "type") = "subjects"
    return(SSDdata)
  }
  #filter the signal for the target frequencies:
  SSDdata = list( signal = filter.apply( data, SSDcoeffs[[1]] ) )
  #filter noise with passband first:
  SSDdata$noise = filter.apply( data, SSDcoeffs[[2]] )
  #afterwards apply stopband:
  SSDdata$noise = filter.apply( SSDdata$noise, SSDcoeffs[[3]] )
  attr(SSDdata, "type") = "SSD"
  return(SSDdata)
}

SSD.apply <- function(SSDdata, nCores = NULL) {
  ## Spatio-Spectral decomposition (SSD), cf. Nikulin et al. 2011
  ## extraction of neuronal oscilations with improved SNR
  #NOTES ---
  # to later project components (S) back to channel space (X):
  # X = X_SSD %*% t(A_SSD)
  # SSD to denoise: project part of the components back:
  # X[,1:n] = X_SSD[,1:n] %*% t(A_SSD[,1:n])
  #INPUT ---
  #SSDdata: list of 2 (1: signal, 2: noise)
  #nCores: if data is a subject list, number of CPU cores to use for parallelization
  #        if NULL, automatic selection; if 1, sequential execution
  #        can also be an already registered cluster object  
  #RETURNS ---
  #W_SSD: spatial filters, demixing matrix W
  #A_SSD: activation patterns, mixing matrix A
  #X_SSD: SSD components, obtained by X*W (component space!)
  #D_SSD: eigenvalues of the SSD components a.k.a. the power ratio corresponding
  #       to the components (ordered like the components)
  .eval_package("geigen")
  if ( "subjects" %in% attr(SSDdata, "type") ) { #parallelize subjects
    pcheck = .parallel_check(required=length(SSDdata), nCores=nCores)
    SSDresult = foreach(d=SSDdata, .combine=list, .multicombine=T) %dopar%
      SSD.apply(d)
    .parallel_check(output=pcheck)
    SSDresult = setNames( SSDresult, paste0("subject", seq_along(SSDdata)) )
    attr(SSDresult, "type") = "subjects"
    return(SSDresult)
  }
  if ( !"SSD" %in% attr(SSDdata, "type") && length(SSDdata) != 2 ) {
    stop( "SSDdata must be alist with 2 elements: signal and noise. ",
          "Refer to the documentation or use SSD.filter." )
  }
  Xs = data.check(SSDdata[[1]], aslist=F)
  Xn = data.check(SSDdata[[2]], aslist=F)
  nsamples = data.samplenum(Xs)
  datacols = is.datacol(Xs)
  norm.cov <- function(X) {
    #population covariance (N) on centered data
    return( cov.wt(scale(X, scale=F), method="ML")$cov )
  }
  if (nsamples > 1) { #compute trial-averaged covariance matrices
    Xs.trials = data.split_trials(Xs, strip=T)
    Xn.trials = data.split_trials(Xn, strip=T)
    Cs = Reduce("+", lapply(Xs.trials, norm.cov))/length(Xs.trials) #compute within trials and average
    Cn = Reduce("+", lapply(Xn.trials, norm.cov))/length(Xn.trials) #compute within trials and average
  } else {
    Cs = norm.cov(Xs[,datacols])
    Cn = norm.cov(Xn[,datacols])
  }
  #intermediate step to check matrix rank with eigenvalue decomposition:
  #note: eigenvalue decomposition of covariance matrix is numerically more accurate than svd
  VD = eigen(Cs); V = VD$vectors; d = VD$values
  #compute an estimate of rank ('pratically' zero values):
  #indicates the number of linearly independent columns
  r = sum(d > 10^-6*d[1]) #rank = number of nonzero diagonal entries in d
  if ( r < sum(datacols) ) {
    cat( "Data does not have full rank, i.e.", sum(datacols)-r+1,
         "columns are collinear. Computing only",r,"components." ) 
    # dimensionality needs to be reduced accordingly
    # i.e. preserve the eigenvectors associated with > 0 eigenvalues
    # PCA approach:
    # prcomp(Xs)$x = Prinicipal components, alternatively: scale(Xs, scale=F) %*% V
    # the square root of the eigenvalues of the cov matrix are the standard deviations of the PCs
    # compare: prcomp(Xs)$sdev and d^0.5 -> d = explained variance per PC
    # to scale data during projection, scale with 1/sd = 1/d^0.5 = d^-0.5
  }
  #construct whitening matrix for C (diagonalization is done below)
  #V is weighted by d^-0.5 to whiten C (1 on the main diagonal)
  #https://en.wikipedia.org/wiki/Whitening_transformation
  #possibly also dim-reduced if r < ncol(C)
  M = V[,1:r] %*% diag(d[1:r]^-0.5) #note: Xwhite = scale(Xs, scale=F) %*% M (PCA)

  #new diagonalized and whitened Cs (Cn is scaled accordingly):
  #V' * C * V = this step diagonalizes C (0 off the diagonal)
  #i.e. redundant information between channels (correlations) is removed (set to 0)
  Cs.white = t(M) %*% Cs %*% M 
  Cn.white = t(M) %*% Cn %*% M
  
  #do the generalized eigenvalue decomposition:
  WD = geigen::geigen(Cs.white, Cs.white+Cn.white) 
  W = WD$vectors[,order(WD$values, decreasing=T)] #descending order
  D = sort(WD$values, decreasing=T)
  
  #project spatial filters (weight vectors) back to original (un-whitened) space:
  W = M %*% W
  #projection of data onto the SSD weights (=SSD components):
  X.SSD = as.matrix(Xs[,datacols]) %*% W #backward modeling: time-courses of r neural sources are estimated
  #note: X is now in component space (weights were applied to it)
  #patterns = Covariance of data * weights * inverse of covariance matrix of extracted sources
  A = Cs %*% W  %*% solve((t(W) %*% Cs %*% W), tol=10^-30) #pattern matrix
  #compile output
  output = list(components=as.data.frame(X.SSD), filters=W, patterns=A, lambda=D, info=Xs[,!datacols])
  attr(output, "type") = "SSD.result"
  return(output)
}

SSD.denoise <- function(SSD, rank, nCores = NULL) {
  ## low-rank factorization to denoise measurements with SSD filters: X*W*A'
  ## necessary if subsequent analyses should be performed in original input space
  ## e.g. to interpret topographies (cf. Haufe et al., 2014)
  ## similar concept as PCA but transformation has same dimensions as input (but lower rank)
  #INPUT ---
  #SSD: output of SSD.apply
  #rank: rank to approximate, < ncol(X); if = ncol(X), just a back-projection
  #      if < 1, the auto-selection criterion from Haufe et al., 2014 is applied:
  #      all components whose lambda is >= q*IQR(lambda) + 75th.percentile(lambda)
  #      where rank will substitute q
  #NOTES ---
  #X*W*A^T = X*W*W^-T = X*I
  #RETURNS ---
  #denoised data (original sensor space!)
  if ( "subjects" %in% attr(SSD, "type") ) { #parallelize subjects
    pcheck = .parallel_check(required=length(SSD), nCores=nCores)
    data = foreach(result=SSD, .combine=list, .multicombine=T) %dopar%
      SSD.denoise(result, rank)
    .parallel_check(output=pcheck)
    attr(data, "type") = "subjects"
    return(data)
  }
  if ( !"SSD.result" %in% attr(SSD, "type") ) stop( "Please supply the output of SSD.apply." )
  if ( rank > ncol(SSD$components) ) stop( "Rank cannot be larger than the measurement columns." )
  if ( rank < 1 ) {
    #auto-select q*IQR(lambda) + 75th.percentile(lambda), but at least 1 component
    rank = max(1, sum( SSD$lambda >= rank*IQR(SSD$lambda)+quantile(SSD$lambda, probs=.75) ) )
    cat( "Auto-selection criterion yielded", rank, "components.\n" )
  }
  data = as.matrix(SSD$components) %*% SSD$filters[,1:rank] %*% t(SSD$patterns[,1:rank])
  data = data.frame( SSD$info, data )
  attr(data, "rank") = rank
  return(data)
}

