# call custom functions

.loadFuncs <- function(util = F, filt = F, CSP = F, SSD = F, ML = F, plot = F) {
  ## load source files into their own environments and attach (only if not present in namespace)
  #INPUT ---
  #util: utils.R - utilities
  #filt: filterFuncs.R - temporal filtering
  #CSP: CSP.R - common spatial patterns algorithm
  #SSD: SSD.R - spatio-spectral decomposition
  #ML: ML.R - machine learning
  
  if (util & !"decode.utils" %in% search()) {
    decode.utils = new.env()
    sys.source("utils.R", envir=decode.utils)
    attach(decode.utils)
  }
  if (filt & !"decode.filts" %in% search()) {
    decode.filts = new.env()
    sys.source("filterFuncs.R", envir=decode.filts)
    attach(decode.filts)
  }
  if (CSP & !"decode.CSP" %in% search()) {
    decode.CSP = new.env()
    sys.source("CSP.R", envir=decode.CSP)
    attach(decode.CSP)
  }
  if (SSD & !"decode.SSD" %in% search()) {
    decode.SSD = new.env()
    sys.source("SSD.R", envir=decode.SSD)
    attach(decode.SSD)
  }
  if (ML & !"decode.ML" %in% search()) {
    decode.ML = new.env()
    sys.source("ML.R", envir=decode.ML)
    attach(decode.ML)
  }
  if (plot & !"decode.plots" %in% search()) {
    decode.plots = new.env()
    sys.source("plotFuncs.R", envir=decode.plots)
    attach(decode.plots)
  }
}

.unloadFuncs <- function() {
  #unloads all custom environments/functions
  idx = which( grepl("decode\\.", search()) ) #get indices of Environments with "decode."
  sapply(rev(idx), function(i) detach(pos=i)) #unload via idx, go from behind to not screw positions
}

#default call to utils when sourced
.loadFuncs(util=T, filt=T, ML=T, CSP=T, SSD=T, plot=T)