## Common spatial pattern algorithm, cf. Lemm et al. 2005, Blankertz et al. 2008
# do temporal filterering before spatial filtering of CSP is applied 

# ToDo: implement regularization for CSP ? cf. Lu et al. 2010

#make helper functions accessible
source("myFuncs.R")

pipe.CSP <- function(data, method="", frequencies=NULL, SSD=F, ...) {
  ## function that creates CSP filters, or optionally SpecCSP filters, 
  ## and validates the result with subsequent machine learning model fit
  #INPUT ---
  #data: df or list of trials expected with sample num in 1st col, binary outcome in 2nd col
  #method: default is CSP, optionally "SpecCSP" or simply "Spec" (not case sensitive)
  #frequencies: 2 frequency values defining a band-pass to filter the data, 
  #      if unspecified, data is expected to be already bandpass filtered.
  #SSD: if frequencies is specified and SSD is True, SSD will be applied instead of 
  #     the standard bandpass filtering, calling pipe.SSD
  #     note: it is advisable to set lowrank to obtain SSD denoised measurements 
  #           e.g. lowrank = 15, so that the first 15 SSD filters will be used
  #           if lowrank is not specified, the SSD components will be used which means
  #           that subsequent analyses are performed in SSD space.
  #           warning messages from CSP can be ignored as they correspond to the reduced rank
  #IMPORTANT: specify also the srate if you want filtering/SSD!
  #npattern, logtransform: see CSP.apply
  #k, independent, shuffle: see data.folds.split
  #NOTES ---
  #further arguments correspond to the functions that are called optionally:
  #if frequencies is defined: transwidth, srate; see filter.get_coeffs
  #if method = "spec": p, q, prior, steps, srate; see SpecCSP.apply 
  #if you repeat this procedure on the same data (e.g. for cross-validation purposes),
  #it is more sensible to filter the data first and then loop over the filtered data, 
  #instead of filtering the data again on every run.
  #RETURNS ---
  #a list with 3 named elements:
  #summary: vector with the AUC value per fold
  #ML: output of the machine learning functions, see data.fit_model
  #CSP: output of the CSP method, see CSP.apply or SpecCSP.apply
  source("myFuncs.R")
  .loadFuncs(ML = T, CSP=T)
  #if filtering is requested
  if ( !is.null(frequencies) ) {
    if ( !"srate" %in% names( list(...) ) ) {
      warning( "Sampling rate (srate) is not set and will default to 500Hz. Please specify if otherwise." )
    }
    if (SSD) { # use SSD
      .loadFuncs(SSD = T)
      args.in = .eval_ellipsis("pipe.SSD", ...)
      SSDout = do.call(pipe.SSD, modifyList( list(...), list(data = data, frequencies = frequencies ) ))
      if ( is.null( args.in$lowrank ) ) { #use SSD components
        warning( "Using SSD components for subsequent analyses because lowrank is not specified.\n",
                 "Be careful when interpreting results, some inferences may not be valid (e.g. topographies)." )
        data = cbind( data[, !is.datacol(data)], SSDout$X_SSD )
      } else { #use denoised measurements (preferable choice)
        data = cbind( data[, !is.datacol(data)], SSDout$X_denoised )
      }
      
    } else if ( length(frequencies) != 2 ) {
      stop( "2 frequencies need to be specified for the lower and upper limits of the bandpass." )
    
    } else { # use standard bandpass filtering
      .loadFuncs(filt = T)
      args.in = .eval_ellipsis("filter.get_coeffs", ...)
      b = do.call(filter.get_coeffs, modifyList( args.in, list(frequencies=frequencies, ftype="pass") ))
      data = filter.apply(data, b)
    }
  }
  #fold the data for cross-validated CSP filters
  args.in = .eval_ellipsis("data.folds.split", ...)
  foldedData = do.call(data.folds.split, modifyList( args.in, list(data=data) ))
  #method to use:
  methodStr = ifelse( grepl("spec", tolower(method)), "SpecCSP.apply", "CSP.apply" )
  args.in = .eval_ellipsis(methodStr, ...); args.in[["..."]] = NULL
  args.feat = .eval_ellipsis("CSP.get_features", ...)
  args.in = modifyList(args.in, args.feat[3:4])
  #loop over fold data to create CSP features:
  CSPfolds = lapply(foldedData, function(fold) {
    .CSP.create_sets(fold$train, fold$test, args.in, args.feat, methodStr=methodStr)
  })
  #continue fitting ML model to CSP features:
  args.in = list(...)
  if ( eval( args.feat$logtransform ) & eval( args.feat$approximate ) ) {
    args.in = modifyList(args.in, list(scale=F)) #scaling done via logtransformation
  }
  CSPfits = lapply(CSPfolds, function(fold) {
    do.call(data.fit_model, modifyList( args.in, list(trainData=fold$train, 
                                                      testData=fold$test) ))
  })
  return( list(summary=sapply(CSPfits, "[[", "AUC"), 
               ML = CSPfits, CSP = lapply(CSPfolds, "[[", "CSP") ))
}

.CSP.create_sets <- function(train, test, args.csp, args.feat, methodStr="CSP.apply") {
  ## helper to prepare train and test set with CSP procedure
  #args.csp: arguments to CSP function
  #args.feat: arguments to CSP.get_features function
  #methodStr: CSP.apply or SpecCSP.apply
  #RETURNS ---
  #trainData, testData, CSP output
  train = data.check(train, aslist=F)
  test = data.check(test, aslist=F)
  if ( eval(args.feat$approximate) ) { #1 value per trial
    train.nsamples = data.get_samplenum(train) 
    test.nsamples = data.get_samplenum(test)
    train.outcome = train[ seq(1, nrow(train), train.nsamples), 2]
    test.outcome = test[ seq(1, nrow(test), test.nsamples), 2]
    train.n = 1:length(train.outcome)
    test.n = 1:length(test.outcome)
  } else { #dimensionality is the same as before
    train.n = train[,1]; test.n = test[,1] 
    train.outcome = train[,2]; test.outcome = test[,2]
  }
  csp = do.call(methodStr, modifyList( args.csp, list(data=train) ))
  trainData = data.frame(train.n, train.outcome, csp$features)
  testData = data.frame(test.n, test.outcome,
      do.call(CSP.get_features, modifyList( args.feat, list(data=test, filters=csp$filters) )))
  return( list(train=trainData, test=testData, CSP=csp))
}

#CSP paradigm
CSP.apply <- function(data, npattern=3, baseline=NULL, ...) {
  ## main function of the CSP paradigm: applies the spatial filters for feature extraction
  ## features are maximally informative with respect to the contrasted condition
  #INPUT ---
  #data: df or list of trials expected with sample num in 1st col, binary outcome in 2nd col 
  #npattern: patterns to extract * 2 (for each condition), 
  # -> if too low the features may not contain sufficient information for the classifiers
  # -> if too high, risk of overfitting
  #logtransform: should variance features be logtransformed
  #baseline: define last sample of baseline (if any) so that it will be excluded 
  #          for the CSP procedure. If Null, no samples are removed.
  data = data.check(data, aslist=F)
  k = unique(data[,2]) #binary outcome expected in 2nd column
  if (length(k) != 2) { stop( "Either outcome is not binary or not in the 2nd column." ) }
  #remove baseline for CSP procedure
  if ( !is.null(baseline) ) {
    data = data.trials.remove_samples(data, end=baseline)
  }
  #get data into trial format
  nsamples = data.get_samplenum(data)
  target = data[seq(1, nrow(data), nsamples), 2] #2nd col with DV
  trialdata = data.trials.split(data, strip=T)
  
  #get trial-averaged Cov-matrix per condition -> C1, C2
  #ToDo: currently just centered. Scale? or (x'*X / trace) (Lu et al.)?
  norm.cov <- function(X) {
    return(cov(scale(X, scale=F))) #sample covariance (N-1) on centered data
  }
  if (nsamples > 1) {
    C1 = Reduce("+", lapply(trialdata[target==k[1]], norm.cov)) / length(trialdata[target==k[1]]) #averaged C
    C2 = Reduce("+", lapply(trialdata[target==k[2]], norm.cov)) / length(trialdata[target==k[2]]) #averaged C
    #ToDo: add regularization (Lu et al. 2010)
  } else {
    C1 = norm.cov( plyr::rbind.fill( trialdata[target==k[1]] ) )
    C2 = norm.cov( plyr::rbind.fill( trialdata[target==k[2]] ) )
  }
  
  #do the eigenvalue decomposition on C1+C2
  VD = eigen(C1+C2, symmetric=T); V = VD$vectors; d = VD$values
  r = sum(d > 10^-6*d[1]) #estimate rank
  if (r < ncol(V)) { 
    warning( paste("Matrix does not have full rank, i.e.", ncol(V)-r+1 ,"columns are collinear.",
                   "Computing only",r,"components.") ) 
  }
  if (r < 2*npattern) {
    stop( paste("Too few components to calculate", 2*npattern, "filters.") )
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
  A = C1 %*% W  %*% solve((t(W) %*% C1 %*% W)) #pattern matrix
  
  #this should be equivalent (does not check for rank)
  #careful with ordering (descending vs. ascending values)
  #careful w.r. to algorithm used (qz, cholesky - cf. matlab)
  #   require(geigen)
  #   eig = geigen(C1, C1+C2) 
  #   W = eig$vectors[, order(eig$values, decreasing = T)]
  #   A = solve(W) #patterns are in the rows
  
  #obtain the npattern spatial filters: 
  #left side maximizes variance under 1st condition and minimizes for 2nd, vice versa on the right
  #filters: eigenvectors from both ends of the eigenvalue spectrum
  
  #with last n cols (e.g. 62:64)
  filters = W[,c(1:npattern, (ncol(W)-npattern+1):ncol(W))] #first n and last n cols
  patterns = A[,c(1:npattern, (ncol(A)-npattern+1):ncol(A))] #for visualization
  
  #extract log-variance features according to CSP:
  #original signal is projected via filters: linear mapping so that the variance of the... 
  #reduced feature space is maximally informative regarding the contrast
  #each projection captures a different spatial localization
  #the filtered signal is uncorrelated in both conditions
  #the features are extracted per trial
  args.in = .eval_ellipsis("CSP.get_features", ...)
  features = do.call(CSP.get_features, modifyList( args.in, list(data=trialdata, filters=filters) ))
  return(list(outcome=target,
              features=features, 
              filters=filters, 
              patterns=patterns))
}

CSP.get_features <- function(data, filters, approximate=T, logtransform=T) {
  ## project measurements onto spatial CSP filters
  #INPUT ---
  #data: df or list of trials
  #filters: spatial filters, e.g. output by CSP.apply
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

# CSP.get_sliced_features <- function(data) {
#   ## gets the log variance of sliced CSP components
#   #INPUT ---
#   #data: a slice list with CSP components
#   #RETURNS ---
#   #a slice list with 1 row per trial (the logvariance per component)
#   .loadFuncs(ML = T)
#   sliced_logvar = lapply(1:length(data), function(i) {
#     dt = data.trials.split(data[[i]], strip=T)
#     x = t( sapply(dt, function(trial) {
#       log( apply(trial, 2, var) )
#     }) )
#     outcome = data.permute_labels(data[[i]], shuffle=F)$outcome
#     out = data.frame(sample=i, outcome, x)
#   })
#   return(sliced_logvar)
# } 

.CSP.write_patterns <- function(patternlist) {
  ## function to write CSP pattern data into a matlab friendly format
  #INPUT ---
  #patternlist: output from CSP.apply
  dir.create(path="./patterns", showWarnings = F)
  lapply(1:length(patternlist), function(i) {
    write.table(patternlist[[i]], 
                file=paste0("./patterns/pattern",i,".txt"), 
                sep=",", row.names=F, col.names=F)
  })
  return("Done.")
}

SpecCSP.apply <- function(data, npattern=3, p=0, q=1, prior=c(1,srate/2), 
                          steps=3, srate=500, ...) {
  ##Spectrally weighted CSP, cf. Tomioka et al. 2006
  ##if frequency band is unknown, this algorithm tries to do simultaneous
  ##spatio-temporal filter optimization, i.e iterative updating of bandpass-filter and weights
  ##generally outperforms broad-band CSP
  #INPUT ---
  #data: df or trial list with 1st col = samples, 2nd col = outcome
  #npattern: number of CSP weights per class (W ncol = 2*npattern)
  #p: regularization parameter, seq(-1,1,by=0.5) = scaling exponent (see below)
  #q: regularization parameter, seq(0,4,by=0.5) = discriminability (see below)
  #prior: frequency band to search in
  #steps: number of iterations
  #srate: sampling rate of the data
  ## see Tomioka paper page 16, figure 6 for some light on the parameters q and p
  #(p,q) = (0,0) = standard wide-band filtered CSP (alpha is set to 1 on each iteration)
  #(p,q) = (-1,1) = theoretical filter optimum
  #(p,q) = (1,0) = prior filter itself (see Eq. 6: alpha set to 1, beta remains unchanged)
  #(p,q) = (0,1) = elementwise product of Eq. 4 and 6 as in the paper (default here)
  #best performance in the area of p = 0:1, q = 0.5:1.5 (see Fig 6) 
  #however if prior information itself is not useful (broadband filter without specific assumption)
  #the theoretical filter optimum with p=-1 is better (sets beta to 1 each iteration), 
  #or more generally: p < 0 if little prior knowledge is at hand (note: p also depends on q!)
  data = data.check(data, aslist=F)
  p = p + q
  k = unique(data[,2]) #binary outcome expected in 2nd column
  if (length(k) != 2) { stop( "Either outcome is not binary or not in the 2nd column." ) }
  nsamples = data.get_samplenum(data) #samples per trial
  target = data[seq(1, nrow(data), nsamples), 2] #outcome
  trialdata = data.trials.split(data, strip=T) #trial format
  cdata = list(trialdata[target==k[1]],
               trialdata[target==k[2]]) #split trialdata for class
  freqs = (0:(nsamples-1) * (srate/nsamples)) #frequency table
  bands = which(freqs >= prior[1] & freqs <= prior[2]) #idx of frequencies to search
  
  #note: equations in Tomioka et al. 2006 refer to single trial data
  # procedure:
  # Phi(X,w,B) = log w'*X*B*B'*X'*w (Eq. 1)
  # obtain feature vector Phi = log-power
  # train LDA classifier
  # repeat
  
  #### begin optimization of coefficients W (spatial filters) and B (temporal filter) ###
  #Sigma = alpha * V = alpha * x*x' = X*U*U'*B*B'*U*U'*X' #summed for all frequency components k
  #compute the FFT of the class data and each trial separately to obtain the cross-spectrum
  specF = lapply(cdata, function(cd) { #for each class...
    
    #Xfft contains U (FFT) for all trials:
    Xfft = lapply(cd, function(td) { #for each trial...
      apply(td, 2, fft)  #for each channel compute FFT: U
    })
    
    #compute full spectrum F of covariance matrices x*x' for each DFT bin k and trial...
    lapply(bands, function(k) { #for each idx in bands (the kth frequency component)...
      lapply(1:length(Xfft), function(tr) { #for each trial...
        #take the spectrum x of all channels at frequency k and compute cross-spectrum x*x':
        2*Re( as.matrix(Xfft[[tr]][k,])%*% Conj(Xfft[[tr]][k,]) ) #U = Xfft[[tr]]
        #note: conjugate (transpose for complex values) reverses the sign of the imaginary part
      })
    })
  })  #list of 2 classes, list of k frequency components, list of n trials, ch x ch cov matrix   
  
  #compute weighted cross-spectrum matrix V (average over trials)
  V = lapply(1:2, function(c) {
    lapply(specF[[c]], function(sF) { #for every freq component, average:
      Reduce("+", sF) / length(sF) #averaged covariance
    }) #list of 2 classes, list of k frequency components, ch x ch cov matrix
  })
  
  #finding spectral filter coefficients alpha (B) for each freq component and spatial filter w
  #number of filters w initialized at 1, changed to 2*npattern on 2nd iteration
  alpha = lapply(1, function(j) { rep(1, length(bands)) }) #initialize alpha at 1 
  #main list corresponding to J=1 with k items
  
  for (step in 1:steps) { #repeat x times
    Filters = lapply(alpha, function(alphaj) { #iterate over J filters
      
      #get the sensor covariance matrices for each class with alpha * V (summed over k):
      Sigma = lapply(1:2, function(c) { #for each class
        S = lapply(1:length(bands), function(k) { #for each freq component k
          alphaj[k] * V[[c]][[k]] #alpha * V
        })
        S = unname(Reduce("+", S)) #Sigma for class c = sum( alpha(k) * V(k) )
        #note: if names remain in the matrix it might claim to be non-symmetric
      })
      
      #optimizing the spatial filter coefficients w:
      #find a decomposition that is common to both classes (brain states)...
      #aka a set of bases that simultaneously diagonalizes both C matrices (Eq. 2)
#       #problem: doesn't handle rank deficiency
#       eig = geigen(Sigma[[1]], Sigma[[1]]+Sigma[[2]]) #cholesky (symmetric matrices) solution is unique
#       lambda = sort(eig$values, decreasing=F, index.return=T) #ascending order
#       VV = eig$vectors[,lambda$ix] #sorted eigenvectors
      eig = eigen(Sigma[[1]]+Sigma[[2]], symmetric=T)
      d = eig$values
      VV = eig$vectors
      r = sum(d > 10^-6*d[1]) #estimate rank
      if (r < ncol(VV)) { 
        warning( paste("Matrix does not have full rank, i.e.", ncol(VV)-r+1 ,"columns are collinear.",
                       "Computing only",r,"components.") ) 
      }
      if (r < 2*npattern) {
        stop( paste("Too few components to calculate", 2*npattern, "filters.") )
      }
      P = diag(d[1:r]^-0.5) %*% t(VV[,1:r]) #aka M
      S1 = P %*% Sigma[[1]] %*% t(P)
      S2 = P %*% Sigma[[2]] %*% t(P)
      eig = eigen(S1)
      lambda = sort(eig$values, decreasing=F, index.return=T)
      R = eig$vectors[,lambda$ix] #ascending order
      W = t(P) %*% R
      A = Sigma[[1]] %*% W  %*% solve((t(W) %*% Sigma[[1]] %*% W)) #patterns
      
      #top eigenvalue per class:
      lambda = c(lambda$x[1], lambda$x[r]) #min, max
      #retain npattern eigenvectors & -values for each class
      W = list(W[,1:npattern], W[,(ncol(W)-npattern+1):ncol(W)]) #filters W
      P = list(A[,1:npattern], A[,(ncol(A)-npattern+1):ncol(A)]) #patterns P
      list(lambda=lambda, W=W, P=P)
    }) 
    
    lambda = plyr::rbind.fill(lapply(Filters, "[[", 1)) #get lambda
    #get W and P such that lambda is minimal/maximal over j
    minJ = which.min(lambda[,1])
    maxJ = which.max(lambda[,2])
    W = cbind(Filters[[minJ]]$W[[1]], #min
              Filters[[maxJ]]$W[[2]]) #max
    P = cbind(Filters[[minJ]]$P[[1]], #min
              Filters[[maxJ]]$P[[2]]) #max 
    
    #optimization of alpha within each class:
    #calculate across trials mean and var of w-projected cross-spectrum component (s)
    #s(w,alpha) = alpha(k)* w'*V(k)*w see below Eq. 3
    #because we also need the variance in Eq. 3 and not only the mean (which is in V), we...
    #go back to specF which lets us compute the variance over the trials
    alpha = lapply(1:ncol(W), function(j) { #for every spatial filter w(j)...
      coeffs = lapply(1:2, function(c) { #for every class c...
        lapply(1:length(bands), function(k) { #for every frequency band k...
          s = sapply(1:length(cdata[[c]]), function(tr) { #for every trial...
            t(W[,j]) %*% specF[[c]][[k]][[tr]] %*% W[,j] # w'*F(k,t)*w
            #this signal s(w,alpha) is spatio-temporally filtered
          })
          #compute mean and variance (over all trials) of s
          list(mu=mean(s), var=var(s)) #for class c and frequency k
        })
      }) #list of: 2 classes, k frequency bins, mu & var
      
      #Eq. 4 for class +: alpha(k,+) = s(k,+,w) - s(k,-,w) / var(s(k,+,w)) + var(s(k,-,w)) | or 0 if <
      alpha_tmp = sapply(1:2, function(c) {
        sapply(1:length(bands), function(k) {
          #update alpha according to Eq. 4:
          alpha_opt = max( 0, ( (coeffs[[c]][[k]]$mu - coeffs[[3-c]][[k]]$mu) / 
                                  (coeffs[[c]][[k]]$var + coeffs[[3-c]][[k]]$var) ) )
          
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
  }
  #switch around W and P columns to represent the order of class1, class2:
  W = W[,c( (ncol(W)-npattern+1):ncol(W), 1:npattern )]
  P = P[,c( (ncol(P)-npattern+1):ncol(P), 1:npattern )]
  
  args.in = .eval_ellipsis("CSP.get_features", ...)
  features = do.call(CSP.get_features, modifyList( args.in, list(data=trialdata, filters=W) ))
  return(
    list(alpha=rbind(matrix(0,min(bands)-1,2*npattern), 
                     unname(t(plyr::ldply(alpha))), 
                     matrix(0,max(0,nsamples/2-max(bands)+1),2*npattern)),
         filters=W, patterns=P, freqs=freqs[1:(length(freqs)/2 +1)], 
         bands=bands, features=features, outcome=target)
  )
}


SpecCSP.plot <- function(specOut, priorlim = T, ylims=c(0,1), title="") {
  ## plotting the freq response (alpha coefficients) of SpecCSP
  #INPUT ---
  #SpecOut: output from SpecCSP.apply
  #priorlim: if True, x-Axis gets limited to the predefined freq range
  #ylims: y-Axis is the relative power contribution of a freq in the range 0-1
  #title: title of the plot
  
  col1 = rgb(0,0.4470,0.7410) #blue
  col2 = rgb(0.4660,0.6740,0.1880) #green
  bands = 1:length(specOut$freqs)
  if (priorlim) { bands = specOut$bands  }
  if (!is.factor(specOut$outcome)) { 
    specOut$outcome = as.factor(specOut$outcome)
  }
  
  plot(specOut$freqs[bands], specOut$alpha[bands,1], type="l", las=1,
       col=col1, lwd=2, ylim=ylims, xlab="Frequency (Hz)", ylab="Power", main=title)
  lines(specOut$freqs[bands], specOut$alpha[bands,ncol(specOut$alpha)], type="l", 
        col=col2, lwd=2, lty=1)
  #2nd class is first in specCSP:
  legend("topleft", as.character(rev(unique(specOut$outcome))), col=c(col1,col2),
         lty=c(1,1), lwd=c(2,2), bty="n")
}



