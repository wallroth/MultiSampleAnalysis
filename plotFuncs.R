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

decoding.plot_ACC_AVG <- function(result, ribbons=F, timerange=c(-500,2000), onset=T,
                                  xticks=NULL, xlabels=NULL, padjust="none", 
                                  alpha=0.01, sign="*", lwd=1.2, hline=0.5, ylims=NULL,
                                  xlab="Time (ms)", ylab="Decoding Accuracy (AUC)", 
                                  ggTheme=.ggTheme(), cols=c("#0072BD", "#D95319") ) {
  ## provide visualization of per slice averaged decoding result
  #INPUT ---
  #result: list as output by decoding functions with summary and significance
  #ribbons: indicate 95% bootstrap CI of the means via ribbons
  #timerange: actual time range of measurements, i.e time of 1st/last sample
  #           the timepoints are then set as the transition from one slice to the next
  #onset: if True, highlight time 0 with dashed line
  #xticks: ticks to put on the x axis, defaults to all slices
  #        if only xlabels are specified, xticks will be set where xlabels are
  #xlabels: time stamps you want displayed, default is the same as xticks if not NULL,
  #         else every slice_num/10th tick is labelled (to make them readable)
  #padjust: adjustment method for significance test, "none" if already applied
  #alpha: significance level to test against
  #sign: significance indicator sign, pass "" to avoid
  #lwd: line width of the lines
  #hline: dotted horizontal line to indicate chance level, pass NULL to avoid
  #ylims: limits on y-axis, automatically determined with some white space if NULL
  #xlab, ylab: labels of the axes
  #ggTheme: ggplot theme, modify only parts by assigning to a variable first
  #         with ggTheme = .ggTheme() and then changing single values such as
  #         ggTheme + theme(legend.position = "none")
  #cols: colours of the lines, 1st true, 2nd random
  if ( !"summary" %in% names(result) & class(result) != "data.frame" ) {
    stop( "Please supply the list output of the used decoding function ", 
          "with element 'summary' or directly supply the data frame. ")
  }
  if (class(result) == "list") result = result$summary
  if ( !"slice" %in% names(result) ) {
    stop( "Cannot find the column 'slice' in the output summary." )
  }
  .eval_package("ggplot2", load=T)
  # check significance
  sign.avg = decoding.signtest(result) 
  sign.avg$pval = p.adjust(sign.avg$pval, method = padjust)
  sign.avg$id = ifelse(sign.avg$pval < alpha, sign, "") #plot marker
  # average for all subjects
  subj.avg = aggregate(AUC ~ label + slice, data=result, mean)
  
  # get x-axis values: the end of each slice's time range marks a tick
  xlen = max(result$slice) #number of points on x-axis
  timepoints = seq( timerange[1], timerange[2], length.out = xlen+1 )[-1]
  subj.avg$time = timepoints[subj.avg$slice]
  
  if (is.null(xlabels)) {
    if (is.null(xticks)) {
      xticks = timepoints
      xstep = max(1, round( xlen/10 )) #stepsize for displayed labels
      xidx <- seq(1, xlen, xstep) #xticks with labels
      xlabels <- as.character( timepoints )
      xlabels[-xidx] <- ""
    } else {
      xlabels = xticks
    }
  } else if (is.null(xticks)) { #only labels specified
    xticks = xlabels
  }  
  
  if (ribbons) {
    CI = .bootstrap_ci(result) #95% 1000 bootstraps CI
    subj.avg$low = CI$low; subj.avg$high = CI$high
  }
  
  yscale = max(subj.avg$AUC) > 1 #units or %
  maxlim = ifelse(yscale, 100, 1)
  if (is.null(ylims)) {
    # set y-limits         
    yedge = ifelse(yscale, 10, 0.1) #white space below/above
    if (ribbons) {
      ylims = c( max(0, min(subj.avg$low)-yedge ),  min(maxlim, max(subj.avg$high)+yedge ) )
    } else {
      ylims = c( max(0, min(subj.avg$AUC)-yedge ), min(maxlim, max(subj.avg$AUC)+yedge ) )
    }
  }
  ysteps = ifelse(yscale, 5, 0.05)
  ydecimal = ifelse(yscale, -1, 1)
  ybreaks = seq( round(ylims[1], ydecimal), round(ylims[2], ydecimal), by = ysteps )
  ylabels = as.character(ybreaks)
  if ( length(ybreaks) > 10 ) {
    ylabels[ seq(2, length(ybreaks), 2) ] = ""
  }
  ypos = min(maxlim, max(subj.avg$AUC)+ysteps) #significance indicator ypos  
  
  ## plot average for all subjects
  p = ggplot(data=subj.avg, aes(x=time, y=AUC, group=label, fill=label, colour=label)) +
    geom_line(size=lwd) +
    scale_color_manual(values=cols) +
    geom_hline(yintercept=hline , linetype="dotted", size=0.1) +
    coord_cartesian(ylim=ylims) + scale_y_continuous(breaks=ybreaks, labels=ylabels) +
    scale_x_continuous(breaks=xticks, labels=xlabels, expand=c(0,0)) +
    annotate("text", x=timepoints, y=ypos, label=sign.avg$id, col="black", size=5) +
    labs(x=xlab, y=ylab) + 
    ggTheme
  if (ribbons) {
    p = p +  
      geom_ribbon(aes(ymin=low, ymax=high), alpha=0.3, colour=NA, show_guide=F) +
      scale_fill_manual(values=cols)
  }
  if (onset) {
    p = p + geom_vline(xintercept=0 , linetype="dashed", size=0.1)
  }
  p #plot and return if function is assigned to variable
}

decoding.plot_ACC_INDV <- function(result, ribbons=F, multi.layout=c(2,2), timerange=c(-500,2000), 
                                   onset=T, xticks=NULL, xlabels=NULL, padjust="none", 
                                   alpha=0.01, sign="*", lwd=1.2, hline=0.5, ylims=NULL,
                                   xlab="Time (ms)", ylab="Decoding Accuracy (AUC)", 
                                   ggTheme=.ggTheme(), cols=c("#0072BD", "#D95319"),
                                   subj_infopos = NULL, subj_infotxt = "Subject") { 
  ## provide visualization of per slice averaged decoding result for every subject
  #INPUT ---
  #result: list as output by decoding functions with summary and significance
  #ribbons: indicate 95% bootstrap CI of the means via ribbons
  #multi.layout: if subject_plots True, nrow, ncol per page
  #timerange: actual time range of measurements, i.e time of 1st/last sample
  #           the timepoints are then set as the transition from one slice to the next
  #onset: if True, highlight time 0 with dashed line
  #xticks: ticks to put on the x axis, defaults to all slices
  #        if only xlabels are specified, xticks will be set where xlabels are
  #xlabels: time stamps you want displayed, default is the same as xticks if not NULL,
  #         else every slice_num/10th tick is labelled (to make them readable)
  #padjust: adjustment method for significance test, "none" if already applied
  #alpha: significance level to test against
  #sign: significance indicator sign, pass "" to avoid
  #lwd: line width of the lines
  #hline: dotted horizontal line to indicate chance level, pass NULL to avoid
  #ylims: limits on y-axis, automatically determined with some white space if NULL
  #       note: same y scaling for all subject plots
  #xlab, ylab: labels of the axes
  #ggTheme: ggplot theme, modify only parts by assigning to a variable first
  #         with ggTheme = .ggTheme() and then changing single values such as
  #         ggTheme + theme(legend.position = "none")
  #cols: colours of the lines, 1st true, 2nd random
  #subj_infopos: position of the subject info, if NULL as title, 
  #              otherwise specify tuple with x,y positions (data based) 
  #              to make an annotation inside the plot
  #subj_infotxt: Prefix txt of the subject id
  if ( !"summary" %in% names(result) & class(result) != "data.frame" ) {
    stop( "Please supply the list output of the used decoding function ", 
          "with element 'summary' or directly supply the data frame. ")
  }
  if (class(result) == "list") result = result$summary
  if ( !"slice" %in% names(result) ) {
    stop( "Cannot find the column 'slice' in the output summary." )
  }
  if ( !"subject" %in% names(result) | length(unique(result$subject)) < 2 ) {
    stop( "No multiple subject information found." )
  }
  .eval_package( c("ggplot2","gridExtra"), load=T)
  LOSO = any(result$subject < 0) #LOSO Or SS? 
  if (LOSO & ribbons) {
    warning( "For a single instance of LOSO decoding there is no variation, thus no ribbons." )
    ribbons = F
  }
  # place subject info as title?
  plotTitle = length(subj_infopos) != 2
  if ( plotTitle & length( ggTheme$plot.title ) == 0 ) {
    ggTheme$plot.title = element_text(face="bold", 
                                      size = max(ggTheme$axis.title.x$size, 
                                                 ggTheme$axis.title.y$size))
  } else if ( !plotTitle ) {
    ggTheme$plot.title = element_blank()
  }
  # get x-axis values: the end of each slice's time range marks a tick
  xlen = max(result$slice) #number of points on x-axis
  timepoints = seq( timerange[1], timerange[2], length.out = xlen+1 )[-1]
  
  if (is.null(xlabels)) {
    if (is.null(xticks)) {
      xticks = timepoints
      xstep = max(1, round( xlen/10 )) #stepsize for displayed labels
      xidx <- seq(1, xlen, xstep) #xticks with labels
      xlabels <- as.character( timepoints )
      xlabels[-xidx] <- ""
    } else {
      xlabels = xticks
    }
  } else if (is.null(xticks)) { #only labels specified
    xticks = xlabels
  }  
  
  yscale = max(result$AUC) > 1 #units or %
  maxlim = ifelse(yscale, 100, 1)
  if (is.null(ylims)) {
    # set y-limits         
    yedge = ifelse(yscale, 10, 0.1) #white space below/above
    ylims = c( max(0, min(result$AUC)-yedge ), min(maxlim, max(result$AUC)+yedge ) )
  }
  ysteps = ifelse(yscale, 5, 0.05)
  ydecimal = ifelse(yscale, -1, 1)
  ybreaks = seq( round(ylims[1],ydecimal), round(ylims[2],ydecimal), by = ysteps )
  ylabels = as.character(ybreaks)
  if ( length(ybreaks) > 10 ) {
    ylabels[ seq(2, length(ybreaks), 2) ] = ""
  }
  
  #plot individual subjects
  plist = lapply( unique(result$subject), function(i) {
    #crop and avg subj data
    subj.s = aggregate(AUC ~ label + slice, data=result[result$subject==i,], mean)
    subj.s$time = timepoints[subj.s$slice] #add time info
    if (ribbons) {
      CI = .bootstrap_ci( result[result$subject==i,] ) #get CIs
      subj.s$low = CI$low; subj.s$high = CI$high
    }
    #create plot
    ps = ggplot(data=subj.s, aes(x=time, y=AUC, group=label, fill=label, colour=label)) +
      geom_line(size=lwd) +
      scale_color_manual(values=cols) +
      geom_hline(yintercept=hline , linetype="dotted", size=0.1) +
      coord_cartesian(ylim=ylims) + scale_y_continuous(breaks=ybreaks, labels=ylabels) +
      scale_x_continuous(breaks=xticks, labels=xlabels, expand=c(0.01,0.01)) +
      labs(x=xlab, y=ylab, title=paste0(subj_infotxt, i)) + 
      ggTheme
    if (!LOSO) {
      #calculate significance
      sign.s = decoding.signtest( result[result$subject==i, names(result) != "subject"] ) 
      sign.s$pval = p.adjust(sign.s$pval, method = padjust)
      sign.s$id = ifelse(sign.s$pval < alpha, sign, "") #plot marker
      ypos = min(maxlim, max(subj.s$AUC)+ysteps) #significance indicator ypos
      ps = ps + annotate("text", x=timepoints, y=ypos, label=sign.s$id, col="black", size=5)
    }
    if ( !plotTitle ) {
      ps = ps + annotate("text", x=subj_infopos[1], y=subj_infopos[2], 
                         label=paste0(subj_infotxt, i), col="#666666", size=5)
    }
    if (ribbons) {
      ps =  ps + 
        geom_ribbon(aes(ymin=low, ymax=high), alpha=0.3, colour=NA, show_guide=F) +
        scale_fill_manual(values=cols)
    }
    if (onset) {
      ps = ps + geom_vline(xintercept=0 , linetype="dashed", size=0.1)
    }
    ps #return to lapply
  })
  #present plots
  multipage = multi.layout[1] * multi.layout[2] #plots per page
  for (i in seq(1+multipage, length(plist)+multipage, multipage)) {
    do.call("grid.arrange", c(plist[(i-multipage):min( (i-1),length(plist) )], 
                              nrow=multi.layout[1], ncol=multi.layout[2]) )
  }
}

.bootstrap_ci <- function(result) {
  # helper to get stratified 95% bootstrapped CIs
  # 1000 ordinary bootstraps, normal CI
  # stratified for slice and label
  .eval_package("boot")
  mean.strat <- function(d, idx) mean(d[idx])
  
  CI = plyr::rbind.fill( lapply( unique(result$slice), function(slice) {
    plyr::rbind.fill( lapply( unique(result$label), function(label) {
      d = result[ result$slice == slice & result$label == label, "AUC" ]
      strat.boot = boot::boot(d, mean.strat, R=1000, stype="i")
      strat.ci = boot::boot.ci(strat.boot, type="norm")$normal[2:3]
      data.frame(slice, label, low=strat.ci[1], high=strat.ci[2])
    }) )
  }) )
  return( CI )
}

decoding.plot_GAT <- function(GATmatrix, timerange=c(-500,2000), onset=T, 
                              limits=NULL, ticks=NULL, labels=NULL, 
                              xlab="Generalization time", ylab="Training time", 
                              legend="AUC", title="", ggTheme=.ggTheme()) {
  ## plots the GAT matrix as obtained by decoding.GAT
  #INPUT:
  #GATmatrix: the output of decoding.GAT
  #timerange: actual time range of measurements, i.e time of 1st/last sample
  #           the timepoints are then set as the transition from one slice to the next
  #onset: if True, highlight time 0 with dashed lines
  #limits: limits of the colour scale, corresponds to AUC values (between 0 and 1)
  #ticks: ticks to put on the x/y axis, defaults to all slices
  #       if only labels are specified, ticks will be set where labels are
  #labels: time stamps you want displayed, default is the same as ticks if not NULL,
  #        else every slice_num/10th tick is labelled (to make them readable)
  .eval_package("ggplot2", load=T)
  .eval_package("reshape2")
  matlabcols <- c("#000088","#000FFD","#00FFFF","#80FF67","#FFFF00","#FF2000","#7F0000")
  GAT = reshape2::melt( unname( GATmatrix ) )
  # get axis values: the end of each slice's time range marks a tick
  len = max(GAT[,1]) #number of ticks on axis
  #each timepoint marks the transition from one slice to the next
  timepoints = seq( timerange[1], timerange[2], length.out = len+1 )[-1]
  GAT$t1 = timepoints[GAT$Var1]
  GAT$t2 = timepoints[GAT$Var2]
  
  if (is.null(limits)) {
    limits = round( c( min(GAT$value), max(GAT$value) ), 1 )
    if ( limits[1] > min(GAT$value) ) limits[1] = round( min(GAT$value)-0.1, 1 )
    if ( limits[2] < max(GAT$value) ) limits[2] = round( max(GAT$value)+0.1, 1 )
  }
  if (is.null(labels)) {
    if (is.null(ticks)) {
      ticks = timepoints
      step = round( len/10 ) #stepsize for displayed labels
      idx <- seq(1, len, step) #ticks with labels
      labels <- as.character( timepoints )
      labels[-idx] <- ""
    } else {
      labels = ticks
    }
  } else if (is.null(ticks)) { #only labels specified
    ticks = labels
  } 
  if ( length( ggTheme$legend.title ) == 0 ) {
    ggTheme$legend.title = element_text(face="bold", size = max(ggTheme$axis.title.x$size, 
                                                                ggTheme$axis.title.y$size))
  }
  p = ggplot(data=GAT, aes(x=t1, y=t2)) +
    geom_tile(aes(fill=value)) +
    scale_fill_gradientn(colours=colorRampPalette(matlabcols)( diff(limits) * 100 ),
                         limits=c(limits[1],limits[2]+.001), breaks=seq(limits[1], limits[2], 0.1),
                         guide=guide_colorbar(bardwidth=0.7, barheight=20, 
                                              ticks=F, nbin=diff(limits) * 100)) +
    scale_y_continuous(breaks=ticks, labels=labels, expand=c(0,0)) +
    scale_x_continuous(breaks=ticks, labels=labels, expand=c(0,0)) +
    labs(x=xlab, y=ylab, fill=legend, title=title) +
    ggTheme + theme(legend.position="right")
  if (onset) {
    p = p + geom_vline(xintercept=0 , linetype="dashed", size=0.1) +
      geom_hline(yintercept=0 , linetype="dashed", size=0.1)
  }
  p
}

data.plot_conditions <- function(data, columns=NULL, cols=NULL, lwd=1) {
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
  nsamples = data.get_samplenum(data) #samples per trial
  target = data[seq(1, nrow(data), nsamples), 2] #outcome
  trialdata = data.trials.split(data, strip=F) #trial format
  classdata = lapply( classes, function(cl) {
    cdat = trialdata[ target == cl ] #all trials for one outcome
    #compute trial and column average
    rowMeans( do.call(cbind, lapply(cdat, function(d) rowMeans(d[, columns], na.rm=T))) )
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
