# call custom functions

.loadFuncs <- function(util = F, filt = F, CSP = F, SSD = F, ML = F, plot = F) {
  ## load source files into their own environments and attach (only if not present in namespace)
  #INPUT ---
  #util: utils.R - utilities
  #filt: filterFuncs.R - temporal filtering
  #CSP: CSP.R - common spatial patterns algorithm
  #SSD: SSD.R - spatio-spectral decomposition
  #ML: ML.R - machine learning
  
  if (util & !"myUtils" %in% search()) {
    myUtils = new.env()
    sys.source("utils.R", envir=myUtils)
    attach(myUtils)
  }
  if (filt & !"myFilts" %in% search()) {
    myFilts = new.env()
    sys.source("filterFuncs.R", envir=myFilts)
    attach(myFilts)
  }
  if (CSP & !"myCSP" %in% search()) {
    myCSP = new.env()
    sys.source("CSP.R", envir=myCSP)
    attach(myCSP)
  }
  if (SSD & !"mySSD" %in% search()) {
    mySSD = new.env()
    sys.source("SSD.R", envir=mySSD)
    attach(mySSD)
  }
  if (ML & !"myML" %in% search()) {
    myML = new.env()
    sys.source("ML.R", envir=myML)
    attach(myML)
  }
  if (plot & !"myPlots" %in% search()) {
    myPlots = new.env()
    sys.source("plotFuncs.R", envir=myPlots)
    attach(myPlots)
  }
}

.unloadFuncs <- function() {
  #unloads all custom environments/functions
  idx = which(grepl("\\my", search())) #get indices of Environments starting with "my"
  sapply(rev(idx), function(i) detach(pos=i)) #unload via idx, go from behind to not screw positions
}

#default call to utils when sourced
.loadFuncs(util = T)