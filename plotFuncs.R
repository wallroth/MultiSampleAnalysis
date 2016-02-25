# plotting functions using ggplot

# Matlab cols: blue, red, orange, purple, green
# c("#0072BD", "#D95319", "#EDB120", "#7E2F8E", "#77AC30")

.ggTheme <- function() {
  ## generate nice gg theme
  .eval_package("ggplot2", load=T)
  myTheme <- theme_bw() + 
    theme(plot.title = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.line = element_line(size=1, colour="black"), #parallel line to axis (both x,y)
          panel.border = element_blank(), 
          panel.background = element_blank(),
          axis.title.x = element_text(face="bold", size=15), 
          axis.title.y = element_text(face="bold", size=15),# angle=0, vjust=1, hjust=4),
          axis.text.x = element_text(size=14), 
          axis.text.y = element_text(size=14), 
          legend.title = element_blank(), #element_text(face="bold", size=14), 
          legend.text = element_text(size=14),
          legend.background = element_blank(), #element_rect(linetype=1, colour="black"),
          legend.position = c(.9,.9), #"none"
          legend.key=element_blank())
  return(myTheme)
}

plot.decode <- function(result, multi.layout=NULL, average=T, type="line", SE=F, threshold=F, xvar="slice", yvar="AUC", groupvar="label", 
                        id="subject", xrange=c(-500,2000), xticks=NULL, xlabels=NULL, size=NULL, signpos=NULL, 
                        hline=0.5, ylims=NULL, xlab="Time (ms)", ylab="Decoding Accuracy (AUC)", 
                        txt="", txtpos=NULL, ggTheme=.ggTheme(), cols=c("#0072BD", "#D95319"), ...) {
  ## provide visualization of per slice averaged decoding result
  #INPUT ---
  #result: list as output by decoding functions with summary and significance
  #multi.layout: nrow, ncol per page for multiple subject plots
  #              If NULL, multi-subject data will be aggregated for a single plot
  #type: character specifying either line or box for the plot type
  #SE: indicate 95% bootstrap CI of the means via ribbons
  #threshold: indicate 95th quantile of permutation performance
  #xvar, yvar: variable on the x-/y-axis
  #groupvar: grouping variable which separates the x/y data points, can be left empty
  #id: subject column, can be left empty
  #xrange: actual time range of measurements, i.e time of 1st/last sample
  #        the timepoints are then set as the transition from one slice to the next
  #xticks: ticks to put on the x axis, defaults to all slices
  #        if only xlabels are specified, xticks will be set where xlabels are
  #xlabels: time stamps you want displayed, default is the same as xticks if not NULL,
  #         else every slice_num/10th tick is labelled (to make them readable)
  #size: line width (defaults to 1.2) or boxplot width (defaults to .5)
  #signpos: y-position of significance indication (if p.vals!=none)
  #hline: dotted horizontal line to indicate theoretical chance level, pass NULL to prevent
  #ylims: limits on y-axis, automatically determined with some white space if NULL
  #xlab, ylab: labels of the axes
  #txt: text to annotate inside the plot, appears in grey colour
  #txtpos: position of txt
  #ggTheme: ggplot theme, modify only parts by assigning to a variable first
  #         with ggTheme = .ggTheme() and then changing single values such as
  #         ggTheme + theme(legend.position = "none")
  #cols: colours according to groupvar or xvar if no groupvar
  if ( class(result) == "data.frame" ) {
    data = result
  } else { #list input
    if ( "result" %in% names(result) ) {
      data = result$result
    } else { #possibly list input of several subjects
      tmp = .output_merge(result) #try to merge
      if ( !"result" %in% names(tmp) ) stop( "Cannot find result data frame." )
      data = tmp$result
    }
  }
  .eval_package("ggplot2", load=T)
  subjects = !is.null(id) && nzchar(id) && id %in% names(data) #subject id in data
  if (!subjects) {
    id = ".tmp"
    data[[id]] = 1 #temporary column to split if no id supplied
  }
  average = average && is.null(multi.layout) #one plot 
  args.test = .eval_ellipsis("decode.test", ...)
  args.test = lapply(args.test[-1], eval) #evaluate the arguments
  args.test$dv = yvar
  if ( is.null(list(...)$p.vals) ) args.test$p.vals = "none"
  timeDomain = !is.null(groupvar) && nzchar(groupvar) #xAxis = time
  lineplot = tolower(type) == "line" && timeDomain
  if ( is.null(size) ) size = ifelse(lineplot, 1.2, .5)
  ## gather the necessary information for plotting ##
  # CIs:
  if ( SE && lineplot ) { #only relevant if lineplot
    if (args.test$bootstrap.CI) { #obtain them via bootstrap
      .eval_package("foreach", load=T)
      bootfun <- function(xdata) {
        means = replicate(args.test$iterations, {
          #sampling with replacement stratified for id
          x = as.vector( sapply(unique( xdata[[id]] ), function(sub) {
            sample( xdata[[yvar]][ xdata[[id]] == sub ] , replace=T ) }) )
          mean(x)
        })
        return( quantile(means, probs=c(.025,.975)) ) #95% CI
      }
      if (average) {
        tmp = split(data, data[,xvar]) #parallelize xvar (slices)
        pcheck = .parallel_check(required=length(tmp), nCores=args.test$nCores)
        CIs = foreach(xdata=tmp, .combine=rbind, .multicombine=T, .maxcombine=length(tmp)) %dopar%
          {
            t( sapply( unique(xdata[[groupvar]]), function(g) {
              bootfun( xdata[ xdata[[groupvar]] == g, ] ) #for each label
            }) )
          }
      } else { #CIs for each subject individually
        tmp = split(data, data[,id]) #parallelize id (subjects)
        pcheck = .parallel_check(required=length(tmp), nCores=args.test$nCores)
        CIs = foreach(sub=tmp, .combine=list, .multicombine=T, .maxcombine=length(tmp)) %dopar%
          {
            xdata = split(sub, sub[,xvar])
            plyr::rbind.fill.matrix( lapply(xdata, function(x) {
              t( sapply( unique(x[[groupvar]]), function(g) {
                bootfun( x[ x[[groupvar]] == g, ] ) #for each label
              }) )              
            }) )
          }
      }
    } else { #no bootstrap for the CI, +/- 1 SE
      if (average) {
        CIs = aggregate( as.formula(paste(yvar,"~",groupvar,"+",xvar)), data=data, 
                         FUN=function(x) { m = mean(x); se = sd(x)/sqrt(length(x)); return(c(m-1.96*se, m+1.96*se)) } )[[3]]
      } else {
        CIs = lapply( split(data, data[,id]), function(d) {
          aggregate( as.formula(paste(yvar,"~",groupvar,"+",xvar)), data=d, 
                     FUN=function(x) { m = mean(x); se = sd(x)/sqrt(length(x)); return(c(m-1.96*se, m+1.96*se)) } )[[3]] })
      }
    } 
  } # CIs done
  # p-values/95th quantile of permutations:
  signtest = NULL
  if ( timeDomain && (threshold || args.test$p.vals != "none") ) {
    if ( "summary" %in% names(result) ) { #already computed
      signtest = result$summary
      if ( id %in% names(signtest) ) { #contains result of each subject
        if ( average && "summary.overall" %in% names(result) ) { #average already computed
          signtest = result$summary.overall
        } else if (!average) {
          signtest = split(signtest, signtest[,id]) #split for subjects
        } else { #average but no summary.overall present
          signtest = do.call(decode.test, modifyList( args.test, list(result=data, id=id) ))    
        }
      } #else no id present
    } else { #no summary
      if (is.null(multi.layout)) { #also applies to single subject/across/between case
        signtest = do.call(decode.test, modifyList( args.test, list(result=data)) )
      } else { #evaluate significance for each subject individually
        tmp = split(data, data[,id])
        signtest = lapply( tmp, function(d) do.call(decode.test, modifyList( args.test, list(result=d, id="fold") )) )
      }
    }
    if ( args.test$p.vals != "none" && (args.test$adjust != "none" || !is.null(list(...)$alpha)) ) {
      #recompute significance
      if (is.null(multi.layout)) {
        p.adj = p.adjust(signtest$p, method=args.test$adjust)
        tmp = p.adj < args.test$alpha
        signtest$sign = ifelse(tmp, "*", "")
      } else {
        signtest = lapply(signtest, function(st) {
          p.adj = p.adjust(st$p, method=args.test$adjust)
          tmp = p.adj < args.test$alpha
          st$sign = ifelse(tmp, "*", "")
          st
        })
      }
    }
  } # p-values done
  if (timeDomain) { 
    cols = cols[ seq_along( unique(data[[groupvar]]) ) ] #in case more colours than labels
    # xAxis timepoints - the end of each slice's time range marks a tick
    xlen = length( unique(data[[xvar]]) ) #number of points on xAxis
    timepoints = seq( xrange[1], xrange[2], length.out = xlen+1 )[-1]
    #add time info: slices are unsorted and might not even be continuously increasing
    data$time = timepoints[ findInterval( data[[xvar]], unique(data[[xvar]]) ) ]
    if (is.null(xlabels)) { #automatic assignment, every length/10-th point
      if (is.null(xticks)) { #defaults to ticks at all time points
        xticks = timepoints
        xstep = max(1, round( xlen/10 )) #stepsize for displayed labels
        xidx <- seq(1, xlen, xstep) #xticks with labels
        xlabels <- as.character( round(timepoints) )
        xlabels[-xidx] <- ""
      } else { #no labels but ticks
        xlabels = xticks #labels for specified ticks
      }
    } else if (is.null(xticks)) { #only labels specified
      xticks = xlabels #ticks for specified labels
    }
    if (lineplot) {
      # mean performance per slice
      if (average) { #aggregate over subjects (if any)
        data = aggregate( as.formula( paste(yvar,"~",groupvar,"+",xvar,"+ time") ), data=data, mean )
        if (SE) data = data.frame(data, CI=unname(CIs)) #add CI to data as CI.1 and CI.2
        if (threshold) data$threshold = signtest[[ grep("CI",names(signtest))[1] ]][ data[[xvar]] ]
      } else {
        data = aggregate( as.formula( paste(yvar,"~",groupvar,"+",xvar,"+ time +",id) ), data=data, mean )
        if (SE) {
          CIs = plyr::rbind.fill.matrix(CIs) #bind list together
          data = data.frame(data, CI=unname(CIs)) #add CI to data as CI.1 and CI.2
        }
        if (threshold) {
          data$threshold = do.call(c, lapply(signtest, function(d) {
            d[[ grep("CI",names(sub))[1] ]][ data[[xvar]][ data[[id]] == unique(d[[id]]) ] ]
          }) )
        }
      }
    }
  }
  # set yAxis
  yscale = max( data[[yvar]] ) > 1 #units or %
  maxlim = ifelse(yscale, 100, 1)
  if ( is.null(ylims) ) { #set reasonable y-limits         
    yedge = ifelse(yscale, 10, 0.1) #white space below/above
    if (SE && lineplot) { #adjust for SE
      ylims = c( max(0, min(CIs[,1])-yedge),  min(maxlim, max(CIs[,2])+yedge) )
    } else { #lineplot without SE or boxplot
      ylims = c( max(0, min( data[[yvar]] )-yedge ), min(maxlim, max( data[[yvar]] )+yedge ) )
    }
  }
  ysteps = ifelse(yscale, 5, 0.05)
  ydecimal = ifelse(yscale, -1, 1)
  ybreaks = seq( round(ylims[1], ydecimal), round(ylims[2], ydecimal), by=ysteps )
  ylabels = as.character(ybreaks)
  if ( length(ybreaks) > 10 ) { #avoid cluttering
    ylabels[ seq(2, length(ybreaks), 2) ] = "" #remove every 2nd
  }
  if (is.null(signpos)) signpos = ylims[2]-ysteps #significance indicator ypos
  ## preparation end ##
  
  ggLineplot <- function(data, signtest=NULL, txt="") {
    p = ggplot(data=data, aes_string(x="time", y=yvar, group=groupvar, 
                                     fill=groupvar, colour=groupvar)) +
      geom_line(size=size) + scale_color_manual(values=cols) +
      geom_hline(yintercept=hline , linetype="dotted", size=0.1) + #chance level indication
      coord_cartesian(ylim=ylims) + scale_y_continuous(breaks=ybreaks, labels=ylabels) +
      scale_x_continuous(breaks=xticks, labels=xlabels, expand=c(0.01,0.01)) +
      annotate("text", x=txtpos[1], y=txtpos[2], label=txt, col="#666666", size=5) +
      labs(x=xlab, y=ylab) + ggTheme
    if ( xrange[1] < (0 - diff(timepoints[1:2])) ) {
      p = p + geom_vline(xintercept=0 , linetype="dashed", size=0.1) #time 0
    }
    if (threshold) { #add permutation threshold
      p = p + geom_line(aes(x=time, y=threshold), col=cols[length(cols)], size=size/2, linetype="dashed")
    }
    if (!average && is.null(multi.layout)) {
      tmp = aggregate( as.formula( paste(yvar,"~",xvar,"+ time") ), data=data, mean )
      p = p + geom_line(mapping=aes_string(x="time",y=yvar), size=size, col="black", data=tmp, inherit.aes=F)
    }
    if (args.test$p.vals != "none") { #add significance indication
      # p = p + annotate("text", x=timepoints, y=signpos, label=signtest$sign, col="black", size=5)
      signtime = timepoints[ as.logical(nchar(signtest$sign)) ] #significant timepoints
      consec = rle( nchar(signtest$sign) ) #consecutive significant timepoints
      consec = consec$lengths[ as.logical(consec$values) ]
      end = signtime[ cumsum(consec) ] #offsets of significance
      start = signtime[ cumsum(consec)+1-consec ] #onsets of significance
      for (i in seq_along(start)) {
        p = p + geom_segment(x=start[i], xend=end[i], y=signpos, yend=signpos, size=size, col=cols[1])
      }
    }
    if (SE) { #add standard error of the mean
      p = p + geom_ribbon(aes(ymin=CI.1, ymax=CI.2), alpha=0.3, colour=NA, show.legend=F) +
        scale_fill_manual(values=cols)
    }
    return(p)
  }
  ggBoxplot <- function(data, signtest=NULL, txt="") {
    update_geom_defaults("point", list(colour = NULL)) #for boxplot outlier colouring
    if (timeDomain) {
      grouping = interaction( data[[xvar]], data[[groupvar]] )
      p = ggplot(data=data, aes_string(x="time", y=yvar, group=grouping, fill=groupvar, colour=groupvar)) +
        scale_x_continuous(breaks=xticks, labels=xlabels, expand=c(0.01,0.01))
    } else { #no groupvar
      p = ggplot(data=data, aes_string(x=xvar, y=yvar, fill=xvar, colour=xvar))
    }
    p = p + scale_color_manual(values=cols) + scale_fill_manual(values=cols) +
      geom_hline(yintercept=hline, linetype="dotted", size=0.1) + #chance level indication
      geom_boxplot(size=size, outlier.size=size/2) +
      coord_cartesian(ylim=ylims) + scale_y_continuous(breaks=ybreaks, labels=ylabels) +
      annotate("text", x=txtpos[1], y=txtpos[2], label=txt, col="#666666", size=5) +
      labs(x=xlab, y=ylab) + ggTheme
    if ( timeDomain && xrange[1] < (0 - diff(timepoints[1:2])) ) {
      p = p + geom_vline(xintercept=0 , linetype="dashed", size=0.1) #time 0
    }
    if (timeDomain && args.test$p.vals != "none") { #add significance indication
      # p = p + annotate("text", x=timepoints, y=signpos, label=signtest$sign, col="black", size=5)
      signtime = timepoints[ as.logical(nchar(signtest$sign)) ] #significant timepoints
      consec = rle( nchar(signtest$sign) ) #consecutive significant timepoints
      consec = consec$lengths[ as.logical(consec$values) ]
      end = signtime[ cumsum(consec) ] #offsets of significance
      start = signtime[ cumsum(consec)+1-consec ] #onsets of significance
      for (i in seq_along(start)) {
        p = p + geom_segment(x=start[i], xend=end[i], y=signpos, yend=signpos, size=2*size, col=cols[1])
      }
    }
    #add median as black lines
    d = layer_data(p, i=2) #2nd layer because chance hline was added first
    p = p + geom_segment(data=d, aes(x=xmin, xend=xmax, y=middle, yend=middle), colour="black", size=1, inherit.aes=F)
    update_geom_defaults("point", list(colour = "black")) #change back
    return(p)
  }
  
  ## plot 
  plotStr = ifelse(lineplot, "ggLineplot", "ggBoxplot")
  if (average) { #single lineplot aggregated for all subjects
    p = do.call(plotStr, list(data, signtest, txt))
    return(p)
  } else { #no averaging
    if (is.null(multi.layout)) { #single plot
      if ("label" %in% names(data)) data = subset(data, label == "true")
      cols = colorRampPalette(c("gray30","gray90"), alpha=T)(length(unique(data[[id]])))
      ggTheme$legend.position = "none"
      groupvar = id
      data[[id]] = factor(data[[id]])
      p = do.call(plotStr, list(data, signtest, txt))
      return(p)
    } else { #multiple plots
      .eval_package("gridExtra", load=T)
      data = split(data, data[,id])
      if ( length(txt) != length(data) ) txt = paste0(txt, seq_along(data))
      plist = lapply( seq_along(data), function(sub) do.call(plotStr, list( data[[sub]], signtest[[sub]], txt[[sub]] )) )  
      multipage = multi.layout[1] * multi.layout[2] #plots per page
      for (i in seq(1+multipage, length(data)+multipage, multipage)) { #iterate until all plots are created
        do.call("grid.arrange", c( plist[ (i-multipage):min(i-1,length(data)) ], 
                                   nrow=multi.layout[1], ncol=multi.layout[2] ))
      }
      # return(plist)
    }
  }
}

plot.GAT <- function(GAT, contrast="none", xrange=c(-500,2000), onset=T, limits=NULL, ticks=NULL, 
                     labels=NULL, xlab="Generalization time", ylab="Training time", legend="AUC", 
                     ggTheme=.ggTheme(), colScheme=NULL, colNum=100, smooth=F, ...) {
  ## plots the GAT matrix as obtained by decoding.GAT
  #INPUT:
  #GATmatrix: the output of decode.GAT
  #contrast: character naming one of difference, significance, threshold; requires permutation results
  #          - difference: the average random performance is subtracted form the avg true perf
  #          - significance: p-values are computed via decode.test and plotted instead of AUC
  #          - threshold: signifies if the 1-alpha quantile of the permutation performance 
  #                       (either bootstrapped or directly from permutations) was exceeded (i.e. binary)
  #xrange: actual time range of measurements, i.e time of 1st/last sample
  #           the timepoints are then set as the transition from one slice to the next
  #onset: if True, highlight time 0 with dashed lines
  #limits: limits of the colour scale, corresponds to AUC values (between 0 and 1)
  #ticks: ticks to put on the x/y axis, defaults to all slices
  #       if only labels are specified, ticks will be set where labels are
  #labels: time stamps you want displayed, default is the same as ticks if not NULL,
  #        else every slice_num/10th tick is labelled (to make them readable)
  #xlab, ylab, legend: titles of the respective part
  #colScheme: colour scheme of the raster, defaults to matlab style
  #colNum: number with which a colour gradient is created from colScheme, i.e. number of distinct colours
  #smooth: if True, the raster is smoothed via linear interpolation
  #...: further arguments to decode.test (if contrast is significance or threshold)
  .eval_package("ggplot2", load=T); .eval_package("reshape2")
  contrast = tolower(contrast)
  if ( !contrast %in% c("none","difference","significance","threshold") ) {
    stop( "Undefined contrast. Options are 'difference', 'significance', 'threshold'. ")
  }
  if ( is.null(colScheme) && contrast == "none" ) {
    colScheme = c("#000088","#000FFD","#00FFFF","#80FF67","#FFFF00","#FF2000","#7F0000")
  }
  #function to transform individual matrices into a data frame
  unravel <- function(L, key) {
    if ( !any( grepl(key, names(L)) ) && class(L) == "list" ) {
      return( do.call(rbind, lapply(L, unravel, key=key)) ) #go further down
    } else if ( class( L[[ grep(key, names(L)) ]] ) == "list" ) {
      d = do.call(rbind, lapply( L[[ grep(key, names(L)) ]], function(x) reshape2::melt( unname(x) ) ))
    } else {
      d = reshape2::melt( unname(L[[ grep(key, names(L)) ]]) )
    }
    return( data.frame(label=key, d) )
  }
  #evaluate GAT input & compute contrast if applicable
  if ( is.matrix(GAT) ) { #direct input of a matrix
    if ( contrast != "none" ) stop( "Supply the GAT output with permutation results for a contrast." )
    GAT.data = reshape2::melt( unname(GAT) )
  } else if ( contrast == "none" ) { #list input, find the true label average
    if ( "average" %in% names(GAT) ) GAT = GAT$average #full output from decode.GAT
    if ( any( grepl("GAT.true", names(GAT)) ) ) GAT = GAT$GAT.true #output includes permutation average
    if ( class(GAT) == "matrix" ) GAT.data = reshape2::melt( unname(GAT) )
    else GAT.data = aggregate(value ~ Var1+Var2, data=unravel(GAT, "true"), mean) #compute average
  } else { #compute contrast
    args.test = .eval_ellipsis("decode.test", ...)
    args.test = lapply(args.test[-1], eval) #evaluate the arguments
    args.test$between = "Var1"; args.test$dv = "value"; args.test$id = ""
    if ( contrast == "difference" ) { #compute difference between true and random average
      if ( "average" %in% names(GAT) ) { #output from decode.GAT
        GAT.data = reshape2::melt( unname(GAT$average[[1]] - GAT$average[[2]]) )
      } else if ( all( grepl("true|random", names(GAT)) ) ) { #average element directly supplied
        GAT.data = reshape2::melt( unname(GAT[[1]] - GAT[[2]]) )
      } else { #compute the averages
        GAT.true = unravel(GAT, "true")
        GAT.random = unravel(GAT, "random")
        GAT.data = aggregate(value ~ Var1+Var2, data=GAT.true, mean)
        GAT.data$value = GAT.data$value - aggregate(value ~ Var1+Var2, data=GAT.random, mean)$value
      } 
      if ( legend == "AUC" ) legend = "Diff" #change the default
      if ( is.null(colScheme) ) colScheme = c("#000000","#FFFFFF") #white and black
    } else { 
      if ( "individual" %in% names(GAT) ) GAT = GAT$individual
      tmp = rbind( GAT.true = unravel(GAT, "true"), GAT.random = unravel(GAT, "random") )
      GAT.data = aggregate(value ~ Var1+Var2, data=tmp[ tmp$label=="true", ], mean)
      #initialize parallelization outside of multiple decode.test calls
      startPar = !is.list(args.test$nCores) && !(contrast == "threshold" && !args.test$bootstrap.CI)
      if (startPar) {
        pcheck = .parallel_check(required=max(GAT.data[,1]), nCores=args.test$nCores)
        args.test$nCores = list() #prevent re-initialization at every loop step
      }
      if ( contrast == "significance" ) { #calculate p-values
        args.test$bootstrap.CI = F #make sure no unnecessary workload
        #p-values for each training time (identical if computed for test time)
        ptest = do.call(rbind, lapply(unique(tmp$Var2), function(train) {
          do.call(decode.test, modifyList( args.test, list(result = tmp[ tmp$Var2 == train, ]) ))
        }))
        GAT.data$value = ptest[[ ifelse(args.test$adjust == "none", "p", "p.adj") ]] #overwrite AUC with p values
        if ( legend == "AUC" ) legend = "p" #change the default
        if ( is.null(colScheme) ) colScheme = c("#FFFFFF","#000000") #white and black
      } else if ( contrast == "threshold" ) { #indicate where 1-alpha permutation quantile was exceeded
        args.test$p.vals = "none" #make sure no unnecessary workload
        #logical value for each training time (identical if computed for test time)
        thresh = do.call(rbind, lapply(unique(tmp$Var2), function(train) {
          do.call(decode.test, modifyList( args.test, list(result = tmp[ tmp$Var2 == train, ]) ))
        }))
        #overwrite AUC with with logicals indicating where the CI was exceeded
        GAT.data$value = GAT.data$value > thresh[[ paste0("CI", 100-args.test$alpha[1]*100) ]]
        if ( legend == "AUC" ) legend = paste0(">CI", 100-args.test$alpha[1]*100) #change the default
      }
      if (startPar) .parallel_check(output=pcheck)
    }
  }
  # get axis values: the end of each slice's time range marks a tick
  len = max(GAT.data[,1]) #number of ticks on axis
  #each timepoint marks the transition from one slice to the next
  timepoints = seq( xrange[1], xrange[2], length.out = len+1 )[-1]
  GAT.data$t1 = timepoints[GAT.data$Var1] #testing time
  GAT.data$t2 = timepoints[GAT.data$Var2] #training time
  if (is.null(limits)) {
    limits = round( c( min(GAT.data$value), max(GAT.data$value) ), 1 )
    #make sure all values are mapped by the colour range
    if ( limits[1] > min(GAT.data$value) ) limits[1] = round( min(GAT.data$value)-0.1, 1 )
    if ( limits[2] < max(GAT.data$value) ) limits[2] = round( max(GAT.data$value)+0.1, 1 )
  }
  if (is.null(labels)) {
    if (is.null(ticks)) {
      ticks = timepoints
      step = round( len/10 ) #stepsize for displayed labels
      idx <- seq(1, len, step) #ticks with labels
      labels <- as.character( round(timepoints) )
      labels[-idx] <- ""
    } else {
      labels = ticks
    }
  } else if (is.null(ticks)) { #only labels specified
    ticks = labels
  } 
  if ( length( ggTheme$legend.title ) == 0 ) { #create title for legend
    ggTheme$legend.title = element_text(face="bold", size = max(ggTheme$axis.title.x$size, ggTheme$axis.title.y$size))
  }
  if ( contrast == "threshold" ) {
    if ( is.null(colScheme) ) colScheme = c("#000000", "#FFFFFF") #black and white
    colfill = scale_fill_manual(values = colScheme)
  } else {
    colfill = scale_fill_gradientn(colours=colorRampPalette(colScheme)( diff(limits) * colNum ),
                                   limits=c(limits[1],limits[2]+.001), breaks=seq(limits[1], limits[2], 0.1),
                                   guide=guide_colorbar(bardwidth=0.7, barheight=20, ticks=F, nbin=diff(limits) * colNum))
  }
  p = ggplot(data=GAT.data, aes(x=t1, y=t2)) +
    geom_raster(aes(fill=value), interpolate=smooth) + colfill +
    scale_y_continuous(breaks=ticks, labels=labels, expand=c(0,0)) +
    scale_x_continuous(breaks=ticks, labels=labels, expand=c(0,0)) +
    labs(x=xlab, y=ylab, fill=legend) +
    ggTheme + theme(legend.position="right")
  if ( onset && xrange[1] < (0 - diff(timepoints[1:2])) ) {
    p = p + geom_vline(xintercept=0 , linetype="dashed", size=0.1) +
      geom_hline(yintercept=0 , linetype="dashed", size=0.1)
  }
  return(p)
}

plot.average <- function(data, columns=NULL, cols=NULL, lwd=1) {
  ## plot measurements but aggregate over trials (separately for outcome)
  ## conceptually analogous to ERP
  #columns: measurement channels to average over
  #cols: optionally specify the colours for each outcome. Necessary if more than 5 outcomes
  #lwd: line width
  data = data.check(data, aslist=F)
  if (is.null(columns)) columns = which( is.datacol(data) )
  if ( !is.factor(data[,2]) ) data[,2] = as.factor(data[,2])
  classes = levels(data[,2])
  if (is.null(cols)) {
    if ( length(classes) > 5 ) stop("More than 5 distinct outcomes. Please specify colours manually.")
    cols = c("#0072BD", "#D95319", "#EDB120", "#7E2F8E", "#77AC30")
  }
  nsamples = data.samplenum(data) #samples per trial
  target = data[seq(1, nrow(data), nsamples), 2] #outcome
  trialdata = data.split_trials(data, strip=F) #trial format
  classdata = lapply( classes, function(cl) {
    cdat = trialdata[ target == cl ] #all trials for one outcome
    #compute trial and column average
    rowMeans( do.call(cbind, lapply(cdat, function(d) rowMeans( as.matrix(d[, columns]), na.rm=T ))) )
  })
  #calculate min max for y axis:
  temp = do.call(c, classdata)
  ylims = c( min(temp), max(temp) )
  plot( classdata[[1]], type="l", las=1, xlab="Sample", ylab="Average", 
        col=cols[1], lwd=lwd, ylim=ylims )
  for (i in 2:length(classdata) ) {
    lines( classdata[[i]], las=1, xlab="Sample", ylab="Average", col=cols[i], lwd=lwd)
  }
  legend("bottomleft", classes, col=cols[1:length(classes)],
         lwd=rep(lwd, length(classes)), bty="n")
}
