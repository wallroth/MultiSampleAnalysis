## recommended to do temporal filterering first to improve spatial filtering results

# ToDo: implement regularization for CSP ? cf. Lu et al. 2010

data.project <- function(data, weights, approximate=T, logtransform=T) {
  ## project measurements onto weights (spatial CSP/SpecCSP/SPoC filters)
  ## the resulting features are an approximation of the "true" underlying signal,
  ## i.e. with respect to an external target variable
  ## the idea is to weight/filter channels that contain signal+noise in such a way
  ## that hopefully only the signal of interest will be left, cf. Haufe et al. 2014, p.98
  ## each projection captures a different spatial localization
  #INPUT ---
  #data: df or list of trials, slices
  #weights: spatial filters, e.g. output by decompose.CSP
  #approximate: if True, the variance per trial is computed to approximate band power
  #             otherwise, the components are returned (same number of rows as data)
  #logtransform: if True, variance is logtransformed (a form of scaling variance)
  data = data.set_type(data)
  if ( "slices" %in% attr(data, "type") ) {
    data = data.split_trials( data, strip=T )
  } else {
    data = data.check(data, aslist=T, strip=T)
  }
  if ( "filters" %in% names(weights) ) { #input is a list from CSP/SPoC
    weights = weights$filters
  } 
  projection = lapply(data, function(trial) as.matrix(trial) %*% weights)
  if (approximate) { #variance per trial
    features = t( sapply(projection, apply, 2, var) ) #column-wise
    if (logtransform) { 
      if ( identical(min(features), 0) ) { #prevent -Inf with log(0)
        #replace 0 with the minimum of 10^-30 or smallest features value after exclusion of 0
        features[ which(features<=0) ] = min( 10^-30, min( features[ which(features>0) ] ) )
      }
      features = log(features)
    }
  } else { #components
    features = plyr::rbind.fill.matrix( projection ) #rebind list to matrix
  }
  return(features)
}

decompose.CSP <- function(data, npattern=3, baseline=NULL, features=T, nCores=NULL, ...) {
  ## Common spatial pattern algorithm, cf. Lemm et al. 2005, Blankertz et al. 2008
  ## CSP paradigm: create spatial filters (weights) which can be used
  ## to extract signals of interest (via projection) from the data measurements
  ## resulting features are maximally informative with respect to the contrasted binary classes
  #INPUT ---
  #data: df or list of trials, slices, subjects 
  #npattern: patterns to extract * 2 (for each condition), 
  # -> if too low the features may not contain sufficient information for the classifiers
  # -> if too high, risk of overfitting
  #baseline: define last sample of baseline (if any) so that it will be excluded 
  #          for the CSP procedure. If Null, no samples are removed.
  #features: if T, features are computed and returned
  #nCores: if data is a subject list, number of CPU cores to use for parallelization
  #        if NULL, automatic selection; if 1, sequential execution
  #        if an empty list, an externally registered cluster will be used
  #RETURNS ---
  #list with 3 or 5 elements (if features = T)
  #filters, patterns, lambda, (features, outcome)
  #filters are the weights to be used for projection, patterns are the 
  #corresponding spatial activation pattern (inverse of the filters)
  #lambda are the eigenvalues (effect sizes) corresponding to the filters
  #NOTE: if outcome is not binary, one vs. the rest approach is taken
  data = data.set_type(data)
  if ( "subjects" %in% attr(data, "type") ) { #parallelize subjects
    pcheck = .parallel_check(required=length(data), nCores=nCores)
    args.in = as.list( match.call() )[-c(1,2)]
    #make sure variable names are evaluated before they are passed on
    args.in = lapply(args.in, function(x) if (class(x)=="name") eval(x) else x)
    CSPresult = foreach(d=data, .combine=list, .multicombine=T) %dopar%
      do.call(decompose.CSP, utils::modifyList( args.in, list(data=d) ))
    .parallel_check(output=pcheck)
    CSPresult = setNames( CSPresult, paste0("subject", seq_along(data)) )
    attr(CSPresult, "type") = "subjects"
    return(CSPresult)
  }
  data = data.check(data, aslist=F)
  #remove baseline for CSP procedure
  if ( !is.null(baseline) ) data = data.remove_samples(data, end=baseline)
  k = unique(data[,2]) #check outcome
  if ( is.null(list(...)$silent) && length(k) > 2 ) { #one vs. the rest (OVR)
    cat( "More than two classes; switching to OVR CSP for multi-class solution.\n" )
    #sort k because otherwise its order will conform to the first appearance in data
    OVR = lapply(sort(k), function(class) { #collapse all other classes to one level ("rest")
      d = data.merge_classes(data, labels = k[!k %in% class], verbose=F)
      decompose.CSP(data=d, npattern=npattern, features=F) #... only required for features
    })
    #column bind the per class OVR results
    filters = do.call(cbind, lapply(OVR, "[[", "filters"))
    patterns = do.call(cbind, lapply(OVR, "[[", "patterns"))
    lambda = do.call(c, lapply(OVR, "[[", "lambda"))
  }
  #get data into trial format
  nsamples = data.samplenum(data)
  target = data[seq(1, nrow(data), nsamples), 2] #2nd col with DV
  data = data.split_trials(data, strip=T)
  if ( length(k) == 2 ) { #binary
    #get trial-averaged Cov-matrix per condition -> C1, C2
    #ToDo: currently just centered. (x'*X / trace) (Lu et al.)?
    norm.cov <- function(X) cov(scale(X, scale=F)) #sample covariance (N-1) on centered data
    if (nsamples > 1) { #compute within trials
      C1 = Reduce("+", lapply(data[target==k[1]], norm.cov)) / length(data[target==k[1]]) #averaged C
      C2 = Reduce("+", lapply(data[target==k[2]], norm.cov)) / length(data[target==k[2]]) #averaged C
      #ToDo: add regularization
    } else { #compute across trials
      C1 = norm.cov( plyr::rbind.fill( data[target==k[1]] ) )
      C2 = norm.cov( plyr::rbind.fill( data[target==k[2]] ) )
    }
    #do the eigenvalue decomposeosition on C1+C2
    VD = eigen(C1+C2, symmetric=T); V = VD$vectors; d = VD$values
    r = sum(d > 10^-6*d[1]) #estimate rank
    if ( is.null(list(...)$silent) && r < ncol(V)) { 
      cat( "Data does not have full rank, i.e.", ncol(V)-r+1 ,
           "columns are collinear. Computing only",r,"components.\n" )
    }
    if (r < 2*npattern) {
      warning( "Too few data columns to calculate ", 2*npattern, " filters. ",
               "Computing ", 2*floor(r/2), " instead." )
      npattern = floor(r/2)
    }
    #whiten C1+C2 via matrix P, such that P * (C1+C2) * P' = I
    P = diag(d[1:r]^-0.5) %*% t(V[,1:r]) #note: called M in the SSD implementation
    #diag((P %*%(C1+C2) %*% t(P))) #whitened
    #whitened spatial covariance matrices (if added together that is):
    S1 = P %*% C1 %*% t(P) #P * C1 * P' 
    S2 = P %*% C2 %*% t(P) #P * C2 * P'
    #eigenvectors of S1 maximize the variance for class 1 and minimize it for class 2
    R = eigen(S1)$vectors #S1 = R * D * R' #R = B in Lu et al.
    #spatial filters:
    W = t(P) %*% R #W' = R' * P
    A = C1 %*% W  %*% solve((t(W) %*% C1 %*% W), tol=10^-30) #pattern matrix
    #obtain the npattern spatial filters: 
    #left side maximizes variance under 1st condition and minimizes for 2nd, vice versa on the right
    filters = W[,c(1:npattern, (ncol(W)-npattern+1):ncol(W))] #first n and last n cols
    patterns = A[,c(1:npattern, (ncol(A)-npattern+1):ncol(A))] #for visualization
    lambda = d[c(1:npattern, (length(d)-npattern+1):length(d))] #to estimate correlation
  }
  output = list(filters=filters, patterns=patterns, lambda=lambda)
  attr(output, "type") = "CSP"
  if (!features) return(output) #else...
  #create CSP features:
  args.in = .eval_ellipsis("data.project", ...)
  output$features = do.call(data.project, modifyList( args.in, list(data=data, weights=filters) ))
  output$outcome = target #add for convenience (directly corresponds to the features)
  return(output)
}

decompose.SPoC <- function(data, npattern=3, baseline=NULL, features=T, nCores=NULL, ...) {
  ## Source Power Correlation Analysis, cf. Dähne et al. 2014
  ## supervised learning approach that takes the DV into account
  ## analogous to CSP but with continuous outcome
  ## general idea: find a spatial filter that extracts an oscillatory signal whose
  ## power correlates with a given (continuous) target variable
  #INPUT ---
  #data: df or list of trials, slices, subjects
  #npattern: the first/last n SPoC filters to use for feature generation
  #baseline: define last sample of baseline (if any) so that it will be excluded 
  #          for the SPoC procedure. If Null, no samples are removed.  
  #features: if T, features are computed and returned
  #nCores: if data is a subject list, number of CPU cores to use for parallelization
  #        if NULL, automatic selection; if 1, sequential execution
  #        if an empty list, an externally registered cluster will be used
  #NOTE: outcome should be numeric, otherwise it will be converted to numbers
  #RETURNS ---
  #filters: spatial filters, demixing matrix W
  #patterns: activation patterns, mixing matrix A
  #lambda: eigenvalues (lambda) of SPoC components, the power correlations
  #features: trial-wise (log)variance of SPoC-components
  #outcome: as supplied (in numeric form), format as 1 value per trial
  data = data.set_type(data)
  if ( "subjects" %in% attr(data, "type") ) { #parallelize subjects
    pcheck = .parallel_check(required=length(data), nCores=nCores)
    SPoCresult = foreach(d=data, .combine=list, .multicombine=T) %dopar%
      do.call(decompose.SPoC, utils::modifyList( list(data=d, npattern=npattern, 
                            baseline=baseline, features=features), list(...) ))
    .parallel_check(output=pcheck)
    SPoCresult = setNames( SPoCresult, paste0("subject", seq_along(data)) )
    attr(SPoCresult, "type") = "subjects"
    return(SPoCresult)
  }
  data = data.check(data, aslist=F)
  if ( !is.null(baseline) ) { #remove baseline for SPoC procedure
    data = data.remove_samples(data, end=baseline)
  }
  nsamples = data.samplenum(data)
  outcome = data[seq( 1, nrow(data), by=nsamples ), 2]
  if ( !is.numeric(outcome) ) {
    warning( "Converting outcome to be numeric." )
    outcome = suppressWarnings( as.numeric(outcome) ) #suppress warnings in case of NAs
    if ( any( is.na(outcome) ) ) { #check if NAs were introduced due to non-numeric characters
      outcome = data[seq( 1, nrow(data), by=nsamples ), 2]
      outcome = as.numeric( as.factor(outcome) )
    }
  }
  #scale to zero mean and unit variance (sensible for continuous variables)
  z = scale(outcome) #z same length as trials
  data = data.split_trials(data[,-2], strip=T) #transform to list
  #z is approximated in each trial by the variance of X
  #which is equal to W' * C(t) * W
  #to have C(t) mean-free: C(t) - C [the average covariance across all trials]
  C = Reduce("+", lapply(data, cov)) / length(data) #averaged C
  C_trials = lapply(data, function(tr) cov(tr) - C) #mean-free trial-wise C
  #obtain a z-weighted covariance matrix:
  #Cz expresses the covariance between z and its approximation
  Ct_vec = sapply(C_trials, matrix, nrow=nrow(C)*ncol(C)) #vectorize the trial-wise cov
  Cz = matrix(Ct_vec %*% z, nrow(C), ncol(C)) / length(data)
  #eigenvalues directly express that covariance
  #Cz needs to be whitened to make the generalized eigenvalue problem into an ordinary one
  VD = eigen(C); V = VD$vectors; d = VD$values
  r = sum(d > 10^-6*d[1]) #rank
  if ( is.null(list(...)$silent) && r < ncol(C) ) {
    cat( "Data does not have full rank, i.e.", ncol(C)-r+1 ,
         "columns are collinear. Computing only",r,"components.\n" )
  }
  if (r < 2*npattern) {
    warning( "Too few data columns to calculate ", 2*npattern, " filters. ",
             "Computing ", 2*floor(r/2), " instead." )
    npattern = floor(r/2)
  }
  M = V[,1:r]  %*% diag(d[1:r]^-0.5) #might not be full-rank
  Cz_white = t(M) %*% Cz %*% M
  #now the ordinary eigenvalue decomposeosition:
  WD = eigen(Cz_white); W = WD$vectors; d = WD$values
  W = M %*% W #project back to original (un-whitened) channel space
  #scale eigenvectors to unit variance since eigenvector scaling is arbitrary and not unique:
  W = apply(W, 2, function(w) { w/sqrt(t(w) %*% C %*% w) }) #does nothing if already unit variance
  A = C %*% W  %*% solve((t(W) %*% C %*% W), tol=10^-30) #pattern matrix
  
  filters = W[,c(1:npattern, (ncol(W)-npattern+1):ncol(W))] #first n and last n cols
  patterns = A[,c(1:npattern, (ncol(A)-npattern+1):ncol(A))] #for visualization
  lambda = d[c(1:npattern, (length(d)-npattern+1):length(d))] #to estimate correlation
  output = list(filters=filters, patterns=patterns, lambda=lambda)
  attr(output, "type") = "SPoC"
  if (!features) return(output) #else...
  #compute SPoC features/components
  args.in = .eval_ellipsis("data.project", ...)
  output$features = do.call(data.project, modifyList( args.in, list(data=data, weights=filters) ))
  output$outcome = outcome #add for convenience (directly corresponds to the features)
  return(output)
}

decompose.SpecCSP <- function(data, srate, npattern=3, p=0, q=1, prior=NULL, iterations=3,
                              baseline=NULL, features=T, nCores=NULL, ...) {
  ##Spectrally weighted CSP, cf. Tomioka et al. 2006
  ##if frequency band is unknown, this algorithm tries to do simultaneous
  ##spatio-temporal filter optimization, i.e iterative updating of bandpass-filter and weights
  ##generally outperforms broad-band CSP
  #INPUT ---
  #data: df or list of trials, slices, subjects
  #srate: sampling rate of the data
  #npattern: number of CSP weights per class (W ncol = 2*npattern)
  #p: regularization parameter, seq(-1,1,by=0.5) = scaling exponent (see below)
  #q: regularization parameter, seq(0,4,by=0.5) = discriminability (see below)
  #prior: frequency band to search in, defaults to the full spectrum
  #iterations: number of iterations
  #baseline: define last sample of baseline (if any) so that it will be excluded 
  #          for the CSP procedure. If Null, no samples are removed.
  #features: if T, features are computed and returned
  #nCores: if data is a subject list, number of CPU cores to use for parallelization
  #        if NULL, automatic selection; if 1, sequential execution
  #        if an empty list, an externally registered cluster will be used
  ## see Tomioka paper page 16, figure 6 for some light on the parameters q and p
  #(p,q) = (0,0) = standard wide-band filtered CSP (alpha is set to 1 on each iteration)
  #(p,q) = (-1,1) = theoretical filter optimum
  #(p,q) = (1,0) = prior filter itself (see Eq. 6: alpha set to 1, beta remains unchanged)
  #(p,q) = (0,1) = elementwise product of Eq. 4 and 6 as in the paper (default here)
  #best performance in the area of p = 0:1, q = 0.5:1.5 (see Fig 6) 
  #however if prior information itself is not useful (broadband filter without specific assumption)
  #the theoretical filter optimum with p=-1 is better (sets beta to 1 each iteration), 
  #or more generally: p < 0 if little prior knowledge is at hand (note: p also depends on q!)
  #RETURNS ---
  #list with 4 or 6 elements (if features = T)
  #filters, patterns, lambda, alpha, (features, outcome)
  #filters are the weights to be used for projection, patterns are the 
  #corresponding spatial activation pattern (inverse of the filters)
  #lambda are the eigenvalues (effect sizes) corresponding to the filters
  #alpha are the spectral filters that were jointly optimized with W
  #alpha has attributes: "bands" indicates which part of the spectrum was searched
  #                      "frequencies" indicates the frequencies corresponding to bands
  #NOTE: if outcome is not binary, one vs. the rest approach is taken
  data = data.set_type(data)
  args.in = as.list( match.call() )[-c(1,2)]
  #make sure variable names are evaluated before they are passed on
  args.in = lapply(args.in, function(x) if (class(x)=="name") eval(x) else x)
  if ( "subjects" %in% attr(data, "type") ) { #parallelize subjects
    pcheck = .parallel_check(required=length(data), nCores=nCores)
    CSPresult = foreach(d=data, .combine=list, .multicombine=T) %dopar%
      do.call(decompose.SpecCSP, utils::modifyList( args.in, list(data=d) ))
    .parallel_check(output=pcheck)
    CSPresult = setNames( CSPresult, paste0("subject", seq_along(data)) )
    attr(CSPresult, "type") = "subjects"
    return(CSPresult)
  }
  data = data.check(data, aslist=F)
  if ( !is.null(baseline) ) { #remove baseline for CSP procedure
    data = data.remove_samples(data, end=baseline)
  }
  if ( is.null(prior) ) prior = c(0, srate/2) #full spectrum
  p = p + q
  k = unique(data[,2])
  if ( is.null(list(...)$silent) && length(k) > 2 ) { #one vs. the rest (OVR)
    cat( "More than two classes; switching to OVR CSP for multi-class solution.\n" )
    #sort k because otherwise its order will conform to the first appearance in data
    OVR = lapply(sort(k), function(class) { #collapse all other classes to one level ("rest")
      d = data.merge_classes(data, labels = k[!k %in% class], verbose=F)
      do.call(decompose.SpecCSP, utils::modifyList( args.in, list(data=d, features=F, baseline=NULL) ))
    })
    #column bind the per class OVR results
    filters = do.call(cbind, lapply(OVR, "[[", "filters"))
    patterns = do.call(cbind, lapply(OVR, "[[", "patterns"))
    lambda = do.call(c, lapply(OVR, "[[", "lambda"))
    alpha = do.call(cbind, lapply(OVR, "[[", "alpha"))
  }
  nsamples = data.samplenum(data) #samples per trial
  if ( nsamples <= 1 ) stop( "Cannot compute a frequency with less than 2 samples." )
  target = data[seq(1, nrow(data), nsamples), 2] #outcome
  data = data.split_trials(data, strip=T) #trial format
  class.data = list(data[ target==k[1] ], data[ target==k[2] ]) #split data for class
  freqs = (0:(nsamples-1) * (srate/nsamples)) #frequency table
  bands = which(freqs >= prior[1] & freqs <= prior[2]) #idx of frequencies to search
  if ( length(bands <= 1) ) { #adjust prior
    prior = c( freqs[ max( which(freqs <= prior[1]) ) ], freqs[ min( which(freqs >= prior[2]) ) ] )
    bands = which(freqs >= prior[1] & freqs <= prior[2]) #idx of frequencies to search
    warning( "Adjusted the prior to match the frequency spectrum." )
  }
  
  #### begin optimization of coefficients W (spatial filters) and B (temporal filter) ###
  if (length(k) == 2) {
    #note: equations in Tomioka et al. 2006 refer to single trial data
    #Sigma = alpha * V = alpha * x*x' = X*U*U'*B*B'*U*U'*X' #summed for all frequency components k
    #compute the FFT of the class data and each trial separately to obtain the cross-spectrum
    specF = lapply(class.data, function(cd) { #for each class...
      #Xfft contains U (FFT) for all trials:
      Xfft = lapply(cd, function(trial) { #for each trial...
        apply(trial, 2, fft)  #for each channel compute FFT: U
      })
      #compute full spectrum F of covariance matrices x*x' for each DFT bin k and trial...
      lapply(bands, function(k) { #for each idx in bands (the kth frequency component)...
        lapply(seq_along(Xfft), function(trial) { #for each trial...
          #take the spectrum x of all channels at frequency k and compute cross-spectrum x*x':
          2*Re( as.matrix(Xfft[[trial]][k,])%*% Conj(Xfft[[trial]][k,]) ) #U = Xfft[[trial]]
          #note: conjugate (transpose for complex values) reverses the sign of the imaginary part
        })
      })
    })  #list of 2 classes, list of k frequency components, list of n trials, ch x ch cov matrix   
    #compute weighted cross-spectrum matrix V (average over trials)
    V = lapply(1:2, function(class) {
      lapply(specF[[class]], function(sF) { #for every freq component, average:
        Reduce("+", sF) / length(sF) #averaged covariance
      }) #list of 2 classes, list of k frequency components, ch x ch cov matrix
    })
    #find spectral filter coefficients alpha (B) for each freq component and spatial filter w
    #number of filters w initialized at 1, changed to 2*npattern on 2nd iteration
    alpha = list( rep(1, length(bands)) ) #initialize alpha at 1 
    #main list corresponding to J=1 with k items
    for (step in 1:iterations) { #repeat x times
      Filters = lapply(alpha, function(alphaj) { #iterate over J filters
        #get the sensor covariance matrices for each class with alpha * V (summed over k):
        Sigma = lapply(1:2, function(class) { #for each class
          S = lapply(1:length(bands), function(k) { #for each freq component k
            alphaj[k] * V[[class]][[k]] #alpha * V
          })
          S = unname(Reduce("+", S)) #Sigma for class = sum( alpha(k) * V(k) )
          #note: if names remain in the matrix it might claim to be non-symmetric
        })
        #optimizing the spatial filter coefficients w:
        #find a decomposeosition that is common to both classes (brain states)...
        #aka a set of bases that simultaneously diagonalizes both C matrices (Eq. 2)
        eig = eigen(Sigma[[1]]+Sigma[[2]], symmetric=T)
        d = eig$values; VV = eig$vectors
        r = sum(d > 10^-6*d[1]) #estimate rank
        if ( is.null(list(...)$silent) && r < ncol(VV) ) {
          cat( "Data does not have full rank, i.e.", ncol(VV)-r+1,
               "columns are collinear. Computing only",r,"components.\n" )
        }
        if (r < 2*npattern) {
          warning( "Too few data columns to calculate ", 2*npattern, " filters. ",
                   "Computing ", 2*floor(r/2), " instead." )
          npattern = floor(r/2)
        }
        P = diag(d[1:r]^-0.5) %*% t(VV[,1:r]) #aka M
        S1 = P %*% Sigma[[1]] %*% t(P)
        S2 = P %*% Sigma[[2]] %*% t(P)
        eig = eigen(S1)
        lambda = sort(eig$values, decreasing=F, index.return=T)
        R = eig$vectors[,lambda$ix] #ascending order
        W = t(P) %*% R
        A = Sigma[[1]] %*% W  %*% solve((t(W) %*% Sigma[[1]] %*% W), tol=10^-30) #patterns
        #retain npattern eigenvectors & -values for each class
        W = list(W[,1:npattern], W[,(ncol(W)-npattern+1):ncol(W)]) #filters W
        P = list(A[,1:npattern], A[,(ncol(A)-npattern+1):ncol(A)]) #patterns P
        d = list(lambda$x[1:npattern], lambda$x[(r-npattern+1):r]) #correlation d
        #top eigenvalue per class:
        lambda = data.frame(min=lambda$x[1], max=lambda$x[r]) #min, max
        list(lambda=lambda, W=W, P=P, d=d)
      }) #end of Filters
      lambda = plyr::rbind.fill(lapply(Filters, "[[", 1)) #get lambda
      #get W and P such that lambda is minimal/maximal over j
      minJ = which.min(lambda$min); maxJ = which.max(lambda$max)
      W = cbind(Filters[[minJ]]$W[[1]], #min
                Filters[[maxJ]]$W[[2]]) #max
      P = cbind(Filters[[minJ]]$P[[1]], #min
                Filters[[maxJ]]$P[[2]]) #max 
      #optimization of alpha within each class:
      #calculate across trials mean and var of w-projected cross-spectrum component (s)
      #s(w,alpha) = alpha(k)* w'*V(k)*w see below Eq. 3
      #because we also need the variance in Eq. 3 and not only the mean (which is in V):
      #go back to specF which lets us compute the variance over the trials
      alpha = lapply(1:ncol(W), function(j) { #for every spatial filter w(j)...
        coeffs = lapply(1:2, function(class) { #for every class...
          lapply(1:length(bands), function(k) { #for every frequency band k...
            s = sapply(1:length(class.data[[class]]), function(tr) { #for every trial...
              t(W[,j]) %*% specF[[class]][[k]][[tr]] %*% W[,j] # w'*F(k,t)*w
              #this signal s(w,alpha) is spatio-temporally filtered
            })
            #compute mean and variance (over all trials) of s
            list(mu=mean(s), var=var(s)) #for class and frequency k
          })
        }) #list of: 2 classes, k frequency bins, mu & var
        #Eq. 4 for class +: alpha(k,+) = s(k,+,w) - s(k,-,w) / var(s(k,+,w)) + var(s(k,-,w)) | or 0 if <
        alpha_tmp = sapply(1:2, function(class) {
          sapply(1:length(bands), function(k) {
            #update alpha according to Eq. 4:
            alpha_opt = max( 0, ( (coeffs[[class]][[k]]$mu - coeffs[[3-class]][[k]]$mu) / 
                                    (coeffs[[class]][[k]]$var + coeffs[[3-class]][[k]]$var) ) )
            #calculate prior filter according to Eq. 6: beta(k) = (s(k,+,w) + s(k,-,w))/2
            beta_k = (coeffs[[1]][[k]]$mu + coeffs[[2]][[k]]$mu)/2
            #plug everything into Eq. 6: alpha_opt(k,c)^q * beta(k)^p
            alpha_opt^q * beta_k^p
          })
        })
        #update with the maximum for both classes
        alphamax = apply(alpha_tmp, 1, max)
        #normalize alpha coefficients so that they sum to unity (1)
        alphamax / sum(alphamax)
      })
    } #end of iteration
    #reverse W and P columns to represent the order of class1, class2:
    revorder = (2*npattern):1
    filters = W[, revorder]; patterns = P[, revorder]
    lambda = c( Filters[[maxJ]]$d[[2]][npattern:1], Filters[[minJ]]$d[[1]][npattern:1] )
    #assemble alpha for output (order still corresponds to min->max):
#     alpha = rbind(matrix(0, nrow = min(bands)-1, ncol= 2*npattern), #append 0s for left-out bands 
#                   unname(t(plyr::ldply(alpha[revorder]))), #insert alpha
#                   matrix(0, nrow = max(0,nsamples/2-max(bands)+1), ncol = 2*npattern)) #append 0s
    alpha = unname(t(plyr::ldply(alpha[revorder]))) #no 0s appended
  }
  row.names(alpha) = bands #indicate searched frequency bands
  attr(alpha, "frequencies") = freqs[bands] #indicate corresponding frequencies
  attr(alpha, "npattern") = npattern #indicate number of patterns per class in alpha
  output = list(filters=filters, patterns=patterns, lambda=lambda, alpha=alpha)
  attr(output, "type") = "SpecCSP"
  if (!features) return(output) #else...
  #compute SpecCSP features/components
  args.in = .eval_ellipsis("data.project", ...)
  output$features = do.call(data.project, modifyList( args.in, list(data=data, weights=output$filters) ))
  output$outcome = target
  return(output)
}

plot.SpecCSP <- function(Spec.result, ylims=c(0,1), title="", 
                         legendpos="topright", cols=NULL, lwd=2, lty=1) {
  ## plots the alpha coefficients of the SpecCSP output
  ## these relfect the signal power per class in the searched frequency spectrum
  #INPUT ---
  #Spec.result: output from decompose.SpecCSP or alpha of that output
  #             if no outcome info is present in the Spec.result input,
  #             classes will be simply numbered
  #ylims: y-Axis is the relative power contribution of a freq in the range 0-1
  #title: title of the plot
  #legendpos: position of the legend
  #cols: colours per class, defaults up to 5 distinct colours, 
  #      define manually if more than 5 classes
  #lwd: line width, 1 value for all classes
  #lty: either 1 value for all, or vector with 1 value per class
  if ( "SpecCSP" %in% attr(Spec.result, "type") ) {
    alpha = Spec.result$alpha #list input
  } else { #just alpha
    alpha = Spec.result
  }
  if ( "outcome" %in% names(Spec.result) ) {
    classes = unique( Spec.result$outcome )
  } else { #no outcome info, simply number the classes
    classes = 1:(ncol(alpha)/attr(alpha, "npattern"))
  }
  if (is.null(cols)) { #manually set colours
    if ( length(classes) > 5 ) stop("More than 5 distinct outcomes. Please specify colours manually.")
    cols = c("#0072BD", "#D95319", "#EDB120", "#7E2F8E", "#77AC30")
  }
  if ( length(lty) < length(classes) ) lty = rep(lty[1], length(classes))
  #plot most distinctive patterns per class:
  if ( length(classes) > 2 ) { #OVR
    patterns = seq(1, ncol(alpha), by=2*attr(alpha, "npattern")) #first of each OVR
  } else {
    patterns = c(1, 2*attr(alpha, "npattern")) #first, last
  }
  plot(attr(alpha, "frequencies"), alpha[,1], type="l", las=1, lwd=lwd, lty=lty[1],
       col=cols[1], ylim=ylims, xlab="Frequency (Hz)", ylab="Power", main=title)
  for ( class in 2:length(classes) ) {
    lines(attr(alpha, "frequencies"), alpha[, patterns[class] ], 
          col=cols[class], lwd=lwd, lty=lty[class])  
  }
  legend(legendpos, as.character(classes), col=cols[1:length(classes)],
         lty=lty, lwd=rep(lwd, length(classes)), bty="n")
}


# .CSP.write_patterns <- function(patternlist) {
#   ## function to write CSP pattern data into a matlab friendly format
#   #INPUT ---
#   #patternlist: output from decompose.CSP
#   dir.create(path="./patterns", showWarnings = F)
#   lapply(1:length(patternlist), function(i) {
#     write.table(patternlist[[i]], 
#                 file=paste0("./patterns/pattern",i,".txt"), 
#                 sep=",", row.names=F, col.names=F)
#   })
#   return("Done.")
# }

