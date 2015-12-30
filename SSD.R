### More general terminology explained (cf. Haufe et al. 2014)
# brain signals are measured in forms of activation patterns at channels + noise
# A = activation patterns; s = signals (isolated responses to stimulus), e = noise, X = measurement/data
# x = A*s(n)' + e(n) = forward (generative) model; s and e have length N (samples)
# each column of A represents the activation pattern for a certain signal
# a column vector of A (a) represents the strength of the signal in each channel
# A therefore has dimensions: channels x signals
# X is then a combination of activation patterns of superimposed signals (+ noise),
# e.g. for 1 signal: A(1:64,1) * S(1:N,1)' = [1:64,1:N] + e(1:64,1:N)
# for 5 signals: A(1:64,1:5) * S(1:N,1:5)' = [1:64,1:N] + e(1:64,1:N)
# -> values in A reflect the signal strength at different channels

# backward modeling tries to extract signals from measurements (reverse direction)
# W'*x(n) = s(n) = backward (discriminative) model
# W = transformation matrix/weights/spatial filters with dimensions: channels x signals
# usually we can assume signals < channels
# by projecting x(n) onto a column of W, each column extracts one latent signal
# e.g. W(1:64,1:5)' * X(1:N,1:64) = [1:5,1:64]
# the resulting rows are components which maximize similarity to a target (signal)
# and are a more informative low-dimensional representations of the overall data
# filters should amplify signal of interest and suppress signals of no interest
# however W cannot be directly interpreted as they are a function of signal + noise
# -> transformation to A required!
# if W is a square matrix (64x64), A = the transpose of the inverse of W
# if W is not square because signals < channels, W is no longer invertible
# W must have full rank (signals are independent) & we assume noise to be uncorrelated to signal
# A = Cx * W * Cs^-1
# Cx = covariance of data; Cs = covariance of signals
# Cs is calculated with W' * Cx * W and has dimensions signals x signals
# -> W[1:64,1:5]' * Cx[1:64,1:64] * W[1:64,1:5] = Cs[1:5,1:5]
# -> Cx[1:64,1:64] * W[1:64,1:5] * Cs[1:5,1:5]^-1 = A[1:64,1:5]
# note: given A, we can also compute W = Cx^-1 * A * Cs
# if signals/components are uncorrelated (orthogonal): A = Cx * W (true for SSD?)


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

#make helper functions accessible
source("myFuncs.R")


pipe.SSD <- function(data, ..., plotfreq = T, lowrank = NULL) {
  ## Pipeline for the whole SSD procedure
  #INPUT ---
  #data: continuous df or list of trials
  #frequencies, transwidth, srate: see SSD.get_filters
  #plotfreq: if True, frequency response of Signal/Noise is plotted
  #lowrank: int, if specified, SSD output is used to denoise measurements by low-rank
  #         factorization (measurements are then in original sensor space unlike components!)
  #RETURNS ---
  #same as SSD; if lowrank was specified additionally the denoised data
  
  args.in = .eval_ellipsis("SSD.get_filters", ...)
  #get the filters:
  SSDfilts = do.call(SSD.get_filters, args.in)
  #get signal and noise:
  Xs_Xn = SSD.get_SigNoise(data, SSDfilts)
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

SSD.get_filters <- function(frequencies, transwidth=1, srate=500, window="blackman") {
  ## get SSD filter coefficients to be supplied to SSD.get_SigNoise
  #INPUT ---
  #frequencies: 6 (or 2) values in this order:
  #   - Signal bandpass: low, high; e.g. 10-12Hz 
  #   - Noise bandpass: low, high, surrounding the signal, e.g. 8-14Hz 
  #   - Noise bandstop: low, high; between signal and noise bp, e.g. 9-13Hz
  #     e.g. frequencies = c(10,12,8,14,9,13)
  #   if only 2 values, they are assumed as signal bandpass specification
  #transwidth: transition width in Hz for the filters to control roll-off
  #srate: sampling rate of the data in Hz
  .loadFuncs(filt = T)
  # check frequency specifications:
  # Noise defaults w.r.t. Signal to +/- 3 Hz bandpass and +/- 1 Hz bandstop
  if (length(frequencies) == 2) {
    print("Only signal bandpass frequencies specified. Assigning the others automatically...")
    frequencies[3] = frequencies[1] - 3
    frequencies[4] = frequencies[2] + 3
    print(paste("Setting noise bandpass to",frequencies[3],"-",frequencies[4],"Hz."))
    frequencies[5] = frequencies[1] - 1
    frequencies[6] = frequencies[2] + 1
    print(paste("Setting noise bandstop to",frequencies[5],"-",frequencies[6],"Hz."))
  }
  if (length(frequencies) != 6) { stop( "Incorrect frequency specification!" ) }
  # get filter coefficients (b):
  Xs_bandpass = filter.get_coeffs(frequencies[1:2], transwidth, ftype="pass", srate=srate, wtype=window)
  Xn_bandpass = filter.get_coeffs(frequencies[3:4], transwidth, ftype="pass", srate=srate, wtype=window)
  Xn_bandstop = filter.get_coeffs(frequencies[5:6], transwidth, ftype="stop", srate=srate, wtype=window)
  return(list(Xs_bandpass=Xs_bandpass, 
              Xn_bandpass=Xn_bandpass, 
              Xn_bandstop=Xn_bandstop))
}

SSD.get_SigNoise <- function(data, SSDfilters) {
  ## create Signal and Noise with supplied data and filter coefficients b
  #INPUT ---
  #data: continuous df or list of trials
  #SSDfilters: filter coefficients, see SSD.get_filters
  #RETURNS ---
  #list of 2 with dfs containing signal and noise, respectively
  source("myFuncs.R")
  .loadFuncs(filt = T)
  #filter the signal for the target frequencies:
  Xs = data.check( filter.apply(data, SSDfilters[[1]]), aslist=F )
  #filter noise with passband first:
  Xn = filter.apply(data, SSDfilters[[2]])
  #afterwards apply stopband:
  Xn = data.check( filter.apply(Xn, SSDfilters[[3]]), aslist=F )
  return(list(signal=Xs, noise=Xn))
}

SSD.apply <- function(SSDdata) {
  ## Spatio-Spectral decomposition (SSD), cf. Nikulin et al. 2011
  ## extraction of neuronal oscilations with improved SNR
  #NOTES ---
  # to later project components (S) back to channel space (X):
  # X = X_SSD %*% t(A_SSD)
  # SSD to denoise: project part of the components back:
  # X[,1:n] = X_SSD[,1:n] %*% t(A_SSD[,1:n])
  #INPUT ---
  #SSDdata: list of 2 (signal + noise), containing continuous DFs
  #RETURNS ---
  #W_SSD: spatial filters, demixing matrix W
  #A_SSD: activation patterns, mixing matrix A
  #X_SSD: SSD components, obtained by X*W (component space!)
  #D_SSD: eigenvalues of the SSD components a.k.a. the power ratio corresponding
  #       to the components (ordered like the components)
  source("myFuncs.R")
  .eval_package("geigen")
  Xs = data.check(SSDdata[[1]], aslist=F, strip=T)
  Xn = data.check(SSDdata[[2]], aslist=F, strip=T)
  #compute covariance matrices
  norm.cov <- function(X) {
    #population covariance (N) on centered data
    return( cov.wt(scale(X, scale=F), method="ML")$cov )
  }
  Cs = norm.cov(Xs)
  Cn = norm.cov(Xn)
  C = Cs
  
  #intermediate step to check matrix rank with eigenvalue decomposition:
  #note: eigenvalue decomposition of covariance matrix is numerically more accurate than svd
  VD = eigen(C); V = VD$vectors; d = VD$values
  #compute an estimate of rank ('pratically' zero values):
  #indicates the number of linearly independent columns
  r = sum(d > 10^-6*d[1]) #rank = number of nonzero diagonal entries in d
  if (r < ncol(Xs)) {
    warning( paste("Matrix does not have full rank, i.e.", ncol(Xs)-r+1 ,"columns are collinear.",
                   "Computing only",r,"components.") ) 
    #covariance matrix C's dimensionality needs to be reduced accordingly
    #i.e. preserve the eigenvectors associated with > 0 eigenvalues
    # = low-rank factorization (same idea as PCA)
    ## PCA approach:
    ## prcomp(Xs)$x = Prinicipal components, alternatively: scale(Xs, scale=F) %*% V
    ## the square root of the eigenvalues of the cov matrix are the standard deviations of the PCs
    ## compare: prcomp(Xs)$sdev and d^0.5 -> d = explained variance per PC
    ## to scale data during projection, scale with 1/sd = 1/d^0.5 = d^-0.5
  }
  #construct whitening matrix for C (diagonalization is done below)
  #V is weighted by λ^-0.5 to whiten C (1 on the main diagonal)
  #https://en.wikipedia.org/wiki/Whitening_transformation
  #possibly also dim-reduced if r < ncol(C)
  M = V[,1:r] %*% diag(d[1:r]^-0.5) #note: Xwhite = scale(Xs, scale=F) %*% M (PCA)

  #new diagonalized and whitened Cs (Cn is scaled accordingly):
  #V' * C * V = this step diagonalizes C (0 off the diagonal)
  #i.e. redundant information between channels (correlations) is removed (set to 0)
  Csr = t(M) %*% Cs %*% M 
  Cnr = t(M) %*% Cn %*% M
  
  #do the generalized eigenvalue decomposition:
  WD = geigen::geigen(Csr, Csr+Cnr) 
  W = WD$vectors[,order(WD$values, decreasing=T)] #descending order
  D = sort(WD$values, decreasing=T)
  
  #project spatial filters (weight vectors) back to original (un-whitened) space:
  W = M %*% W
  #projection of data onto the SSD weights (=SSD components):
  X_SSD = as.matrix(Xs) %*% W #backward modeling: time-courses of r neural sources are estimated
  #note: X is now in component space (weights were applied to it) 

  #patterns = Covariance of data * weights * inverse of covariance matrix of extracted sources
  A = C %*% W  %*% solve((t(W) %*% C %*% W), tol=10^-30) #pattern matrix
  return( list(W_SSD=W, A_SSD=A, X_SSD=as.data.frame(X_SSD), D_SSD=D) )
}

SSD.denoise <- function(Xs, A, k=ncol(A), W=diag( ncol(A) )) {
  ## low-rank factorization to denoise measurements with SSD filters: X*W*A'
  ## necessary if subsequent analyses should be performed in original input space
  ## e.g. to interpret topographies (cf. Haufe et al., 2014)
  ## similar concept as PCA but transformation has same dimensions as input
  #INPUT ---
  #Xs: bandpass-filtered signal, e.g. from SSDdata$signal or SSD components (X_SSD)
  #   i.e. a signal which was already projected by W
  #W: filter matrix, e.g. as obtained from SSD, defaults to identitity matrix.
  #   only required if Xs is not components but some bandpass filtered signal
  #A: pattern matrix, e.g. as obtained from SSD
  #k: rank to approximate, < ncol(X); if = ncol(X), just a projection
  #NOTES ---
  #X*W*A^T = X*W*W^-T = X*I
  #RETURNS ---
  #denoised data (original sensor space!)
  Xs = data.check(Xs, aslist=F) #in case Xs did not come from SSD.get_SigNoise
  X = as.matrix( Xs[,is.datacol(Xs)] )
  if (k > ncol(X)) {
    stop( "Rank k cannot be larger than the columns containing measurements in X." )
  }  
  X_denoised = X %*% W[,1:k] %*% t(A[,1:k])
  return( data.frame( Xs[,!is.datacol(Xs)], X_denoised ) )
}


SPoC.apply <- function(data, outcome, npattern=NULL, k=NULL, baseline=NULL, ...) {
  ## Source Power Correlation Analysis, cf. Dähne et al. 2014
  ## Optimizes SSD filters so that covariance with target variable is maximal
  ## analogous to CSP but with continuous outcome
  ## supervised learning approach that takes the DV into account unlike ICA
  ## general idea: find a spatial filter that extracts an oscillatory signal whose
  ## power correlates with a given (continuous) target variable
  #INPUT ---
  #data: SSD output, continuous data or list of trials, 1st col with sample number
  #outcome: continuous outcome vector with elements corresponding to the number of trials
  #npattern: the first/last n SPoC filters to use for feature generation
  #k: if SSD output is supplied, k is the number of components used for low-rank factorization 
  #baseline: define last sample of baseline (if any) so that it will be excluded 
  #          for the SPoC procedure. If Null, no samples are removed.  
  #RETURNS ---
  #W_SPoC: spatial filters, demixing matrix W
  #A_SPoC: activation patterns, mixing matrix A
  #X_SPoC: SPoC components, obtained by X*W
  #D_SPoC: eigenvalues (lambda) of SPoC components, the power correlations
  #features: trial-wise (log)variance of SPoC-components
  if ( !is.numeric(outcome) ) {
    stop( "Outcome must be numeric." )
  }
  if ( sum( grepl("SSD", names(data)) ) >= 3 ) { #SSD output with at least X, W, A
    SSD = list( W=data$W_SSD, A=data$A_SSD ) #save SSD filters and patterns
    if (is.null(k)) k = ncol(SSD$A) #no denoising
    data = SSD.denoise(data$X_SSD, k=k, A=SSD$A) #project to original measurement space
    #infer sample numbers from outcome length (num trials)
    data = data.append_info(data, trials=length(outcome))
  }
  data = data.check(data, aslist=F)
  if ( (nrow(data)/length(outcome)) %% 1 > 0 ) {
    warning( "Number of outcomes does not seem to match the number of trials or your trials vary in length." )
  }
  #remove baseline for SPoC procedure
  if ( !is.null(baseline) ) {
    data = data.trials.remove_samples(data, end=baseline)
  }
  nsamples = data.get_samplenum(data)
  #get the outcome z, scaled to zero mean and unit variance (sensible for continuous variables)
  z = scale( outcome ) #z same length as trials
  trialdata = data.trials.split(data, strip=T) #transform to list
  #z is approximated in each trial by the variance of X
  #which is equal to W' * C(t) * W
  #to have C(t) mean-free: C(t) - C [the average covariance across all trials]
  C = Reduce("+", lapply(trialdata, cov)) / length(trialdata) #averaged C
  C_trials = lapply(trialdata, function(tr) cov(tr) - C) #mean-free trial-wise C
  
  #obtain a z-weighted covariance matrix:
  #Cz expresses the covariance between z and its approximation
  Ct_vec = sapply(C_trials, matrix, nrow=nrow(C)*ncol(C)) #vectorize the trial-wise cov
  Cz = matrix(Ct_vec %*% z, nrow(C), ncol(C)) / length(trialdata)
  
  #eigenvalues directly express that covariance
  #Cz needs to be whitened to make the generalized eigenvalue problem into an ordinary one
  VD = eigen(C); V = VD$vectors; d = VD$values
  r = sum(d > 10^-6*d[1]) #rank
  if (r < ncol(C) & is.null(k)) {
    warning( paste("Matrix does not have full rank, i.e.", ncol(C)-r+1 ,"columns are collinear.",
                   "Computing only",r,"components.") )   
  }
  M = V[,1:r]  %*% diag(d[1:r]^-0.5) #might not be full-rank
  Cz_white = t(M) %*% Cz %*% M
  #now the ordinary eigenvalue decomposition:
  WD = eigen(Cz_white); W = WD$vectors; d = WD$values
  W = M %*% W #project back to original (un-whitened) channel space
  
  #scale eigenvectors to unit variance since eigenvector scaling is arbitrary and not unique:
  W = apply(W, 2, function(w) { w/sqrt(t(w) %*% C %*% w) }) #does nothing if already unit variance
  A = C %*% W  %*% solve((t(W) %*% C %*% W), tol=10^-30) #pattern matrix
  
#   #if data were SSD components, change W and A to be applicable to original measurements:
#   if ( !is.null(W_SSD) ) {
#     #map filters from SSD space to original space of data
#     W = SSD$W %*% W 
#     A = SSD$A %*% A
#   } 
  if (is.null(npattern)) {
    if (is.null(k)) {
      npattern = floor( ncol(W)/4 )
    } else {
      npattern = floor(k/2)
    }
  }
  filters = W[,c(1:npattern, (ncol(W)-npattern+1):ncol(W))] #first n and last n cols
  patterns = A[,c(1:npattern, (ncol(A)-npattern+1):ncol(A))] #for visualization
  lambda = d[c(1:npattern, (length(d)-npattern+1):length(d))]
  #compute SPoC features/components
  args.in = .eval_ellipsis("SPoC.get_features", ...)
  features = do.call(SPoC.get_features, modifyList( args.in, list(data=trialdata, filters=filters) ))
  return(list(filters=filters, patterns=patterns, lambda=lambda, features=features))
}

SPoC.get_features <- function(data, filters, approximate=T, logtransform=T) {
  ## backward model W'*x(n) = s(n) -> project data onto spatial filters and extract variance
  # the idea here is to weight/filter channels that contain signal+noise so that
  # hopefully only the signal of interest will be left, cf. Haufe et al. 2014, p.98
  # the resulting features are an approximation of the "true" underlying signal,
  # i.e. with respect to an external target variable
  #INPUT ---
  #data: can be continuous or trial format
  #filters: spatial filters obtained from SPoC
  #approximate: if True, the variance per trial is computed to
  #             approximate band power
  #             otherwise, the components are returned  
  #logtransform: if True, variance is logtransformed
  trialdata = data.check(data, aslist=T, strip=T)
  if (approximate) { 
    #variance per trial
    features = t(sapply(trialdata, function(td) apply(as.matrix(td) %*% filters, 2, var)))
    if (logtransform) { 
      features = log(features)
    }
    return(features)
  }
  #components
  features = plyr::rbind.fill.matrix( lapply(trialdata, function(td) as.matrix(td) %*% filters) )
  return(features)
}



