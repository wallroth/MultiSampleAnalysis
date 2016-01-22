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

SSD.pipeline <- function(data, ..., plot=T, nCores=NULL) {
  ## Pipeline for the whole SSD procedure which takes care of the sequence
  ## of functions to execute; especially shines with subject data sets
  ## as the full paradigm gets parallelized. However, only the denoised
  ## data is returned, i.e. no other SSD output. For access to the full SSD
  ## output, refer to the individual functions. The most time consuming
  ## step in the paradigm is the temporal filtering of signal and noise
  ## i.e. the call to SSD.filter. If there is need to repeatedly access
  ## these two, it is computationally much more efficient to call the SSD
  ## functions separately. Such may be the case if the SSD procedure
  ## is nested within a cross-validation loop which repeatedly splits the data.
  ## As the filtering is done within trials there is no harm in pre-filtering
  ## the full data set(s) and thus executing the filter step only once.
  #INPUT ---
  #data: continuous df or list of trials, slices, subjects
  #frequencies, transwidth, srate: see SSD.coefficients
  #plot: if True, frequency response of Signal/Noise is plotted (only if not parallelized)
  #      additional arguments of plot.frequencies may be supplied
  #rank: see SSD.denoise; if not supplied, defaults to q=0.1 for the auto-selection criterion
  #nCores: if data is a subject list, number of CPU cores to use for parallelization
  #        if NULL, automatic selection; if 1, sequential execution
  #        can also be an already registered cluster object
  #RETURNS ---
  #the denoised data as a data frame. Attribute "rank" is added to indicate
  #rank after denoising (useful if auto-selection was chosen)
  data = data.set_type(data)
  if ( "subjects" %in% attr(data, "type") ) { #parallelize subjects
    pcheck = .parallel_check(required=length(data), nCores=nCores)
    data = foreach(d=data, .combine=list, .multicombine=T) %dopar%
      do.call(SSD.pipeline, utils::modifyList( list(data=d, plot=F), list(...) ))
    .parallel_check(output=pcheck)
    data = setNames( data, paste0("subject", seq_along(data)) )
    attr(data, "type") = "subjects"
    return(data)
  }
  #get the filter coefficients:
  args.in = list(...)
  SSDcoeffs = do.call(SSD.coefficients, args.in)
  #get signal and noise:
  SSDdata = SSD.filter(data, SSDcoeffs)
  if (plot) { #plot frequency response for Signal and Noise
    args.plot = .eval_ellipsis("plot.frequencies", ...)
    signal = do.call(plot.frequencies, modifyList(args.plot, list(data=SSDdata$signal, plot=F)))
    noise = do.call(plot.frequencies, modifyList(args.plot, list(data=SSDdata$noise, plot=F)))
    plot(names(signal), signal, type="l", las=1, col="#0072BD", lwd=args.plot$lwd,
         xlab="Frequency (Hz)", ylab="Spectral Power Density (dB/Hz)", 
         main="SSD: frequency response")
    lines(names(noise), noise, col="#D95319", lwd=args.plot$lwd)
    if ( is.null(args.in$transwidth) ) args.in$transwidth = 1
    transition = c( args.in$frequencies[1:2] - args.in$transwidth/2,
                    args.in$frequencies[1:2] + args.in$transwidth/2 )
    abline(v=args.in$frequencies[1:2], col="black", lty=2) #specified frequency
    abline(v=transition, col="grey", lty=3) #edges of the transition zone
    legend("topright", c("signal","noise"), 
           col=c("#0072BD","#D95319"), lwd=args.plot$lwd, bty="n")
  }
  #do the SSD:
  SSD.out = SSD.apply(SSDdata)
  #denoise data:
  data = SSD.denoise(SSD.out, args.in$rank)
  return(data)
}

SSD.coefficients <- function(frequencies, noise.width=2, ...) {
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
  #transwidth, srate, wtype: see filter.coefficients
  #transwidth will default to 1, wtype to blackman
  #RETURNS ---
  #a list with the filter coefficients for the two/three filtering steps of the SSD
  #i.e. signal bandpass, noise bandpass, noise bandstop (optionally)
  #Note: at default settings, the noise edges will directly border on the signal edges
  #      to minimize spectral leakage. The borders are determined via the transition width
  args.coeff = .eval_ellipsis("filter.coefficients", ...)
  if ( is.null(list(...)$wtype) ) { #if not specified...
    args.coeff$wtype = "blackman" #default to blackman for smaller deviation
  }
  if ( is.null(list(...)$transwidth) ) {
    args.coeff$transwidth = 1 #default to 1
  }
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
  if ( !length(frequencies) %in% c(4,6) ) { stop( "Incorrect frequency specification! ",
                                                  "Set either 2, 4 or 6 frequencies." ) }
  #get filter coefficients (b):
  SSDcoeffs = list( signal.pass = do.call( filter.coefficients, modifyList(args.coeff, list(frequencies=frequencies[1:2], ftype="pass")) ), 
                    noise.pass = do.call( filter.coefficients, modifyList(args.coeff, list(frequencies=frequencies[3:4], ftype="pass")) ) )
  if ( length(frequencies) == 6 ) {
    SSDcoeffs$noise.stop = do.call( filter.coefficients, modifyList(args.coeff, list(frequencies=frequencies[5:6], ftype="stop")) )
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
  #filter noise with passband:
  SSDdata$noise = filter.apply( data, SSDcoeffs[[2]] )
  if ( length(SSDcoeffs) > 2 ) { #additionally apply stopband:
    SSDdata$noise = filter.apply( SSDdata$noise, SSDcoeffs[[3]] )  
  } else { #subtract signal from noise if no bandstop supplied:
    SSDdata$signal = data.check(SSDdata$signal, aslist=F)
    SSDdata$noise = data.check(SSDdata$noise, aslist=F)
    measurements = is.datacol(SSDdata$signal)
    SSDdata$noise = data.frame( SSDdata$noise[, !measurements], 
                                as.matrix( SSDdata$noise[, measurements] ) - 
                                  as.matrix( SSDdata$signal[, measurements] ) )
  }
  attr(SSDdata, "type") = "SSD.data"
  return(SSDdata)
}

SSD.apply <- function(SSDdata, nCores=NULL) {
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
  #components: SSD components, obtained by X*W (component space!)
  #filters: spatial filters, demixing matrix W
  #patterns: activation patterns, mixing matrix A
  #lambda: eigenvalues of the SSD components a.k.a. the power ratio corresponding
  #        to the components (ordered like the components)
  #info: info columns which were stripped off data
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
  if ( !"SSD.data" %in% attr(SSDdata, "type") && length(SSDdata) != 2 ) {
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
  attr(output, "type") = "SSD"
  return(output)
}

SSD.denoise <- function(SSD, rank=NULL, nCores=NULL) {
  ## low-rank factorization to denoise measurements with SSD filters: X*W*A'
  ## necessary if subsequent analyses should be performed in original input space
  ## e.g. to interpret topographies (cf. Haufe et al., 2014)
  ## similar concept as PCA but transformation has same dimensions as input (but lower rank)
  #INPUT ---
  #SSD: output of SSD.apply
  #rank: rank to approximate, < ncol(X); if = ncol(X), just a back-projection
  #      if < 1, the auto-selection criterion from Haufe et al., 2014 is applied:
  #      all components whose lambda is >= q*IQR(lambda) + 75th.percentile(lambda)
  #      where rank will substitute q; if NULL, defaults to auto-select with q=0.1
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
  if ( !"SSD" %in% attr(SSD, "type") ) stop( "Please supply the output of SSD.apply." )
  if ( is.null(rank) ) {
    rank = .1 #default q=0.1 for auto-selection
  } else if ( rank > ncol(SSD$components) ) {
    stop( "Rank cannot be larger than the number of measurement columns." )
  }
  if ( rank < 1 ) {
    #auto-select q*IQR(lambda) + 75th.percentile(lambda), but at least 1 component
    rank = max(1, sum( SSD$lambda >= rank*IQR(SSD$lambda)+quantile(SSD$lambda, probs=.75) ) )
    cat( "Auto-selection criterion yielded", rank, "components.\n" )
  }
  #because components are already projected onto W, it is not needed here
  data = as.matrix(SSD$components[,1:rank]) %*% t(SSD$patterns[,1:rank])
  data = data.frame( SSD$info, data )
  attr(data, "rank") = rank
  return(data)
}

