# plotting functions using ggplot

# Matlab cols: blue, red, orange, purple, green
# c("#0072BD", "#D95319", "#EDB120", "#7E2F8E", "#77AC30")

.ggTheme <- function() { 
  ## generate nice gg theme
  .eval_package("ggplot2", load=T)
  ggTheme <- theme_bw() + 
    theme(plot.title = element_text(size=14, face="plain", hjust=.5, margin=margin(b=10)),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.line = element_line(size=.75, colour="black"), #parallel line to axis (both x,y)
          panel.border = element_blank(), 
          panel.background = element_blank(),
          axis.title.x = element_text(face="plain", size=14), 
          axis.title.y = element_text(face="plain", size=14),# angle=0, vjust=1, hjust=4),
          axis.text.x = element_text(size=12, colour="black"), 
          axis.text.y = element_text(size=12, colour="black"), 
          axis.ticks = element_line(colour="black"),
          legend.title = element_text(face="plain", size=13),
          legend.text = element_text(size=12),
          legend.background = element_blank(), #element_rect(linetype=1, colour="black"),
          legend.position = "right", #c(.9,.9), #"none"
          legend.key=element_blank(),
          strip.background = element_rect(fill = "white", colour = "white"),
          strip.text = element_text(size=13, face="plain"))
  return(ggTheme)
}

plot.decode <- function(result, multi.layout=NULL, average=T, lineplot=T, CI=F, hline=NULL, vline=NULL,
                        xvar="slice", yvar="AUC", groupvar="label", idvar="subject", xlab="Time (ms)", ylab="AUC", legend="Labels",
                        xrange=NULL, xticks=NULL, xlabels=NULL, ylims=NULL, yticks=NULL, ylabels=NULL, size=NULL, signpos=NULL, 
                        txt="", txtpos=NULL, txtcol="#666666", txtsize=5, cols=c("#0072BD", "#666666"), ggTheme=.ggTheme(), ...) {
  ## provide visualization of per slice averaged decoding result
  #INPUT ---
  #result: decoding output, or more generally any data frame/data table
  #multi.layout: nrow, ncol per page for multiple subject plots; 
  #              if NULL, multi-subject data will be aggregated for a single plot
  #              if NULL and average FALSE, subjects will all be plotted into one plot
  #average: if True, data will be aggregated over id; if False and multi.layout is not specified,
  #         all individual id data sets will be plotted into one plot (only true label data)
  #         a black line will be added to indicate the mean of the true performance,
  #         the other colours will be shades of grey if cols does not match the unique values of id
  #lineplot: if True, creates a line plot, otherwise a box plot
  #CI: indicate CI of the means via ribbons, either via bootstrap (if TRUE) or standard error of the mean
  #hline, vline: if non-NULL, horizontal/vertical line will be added to the plot at specified coordinate
  #              i.e. to indicate chance/significance level (hline) or stimulus onset time (vline)
  #xvar, yvar: variable on the x-/y-axis
  #groupvar: grouping variable which separates the x/y data points, can be left empty
  #idvar: subject column, can be left empty
  #xlab, ylab, legend: labels of the axes/legend
  #xrange: actual time range of measurements, i.e time of 1st/last sample;
  #        if non-NULL, the time points wil replace the slice numbering on the x-axis
  #        time points will be left-aligned to the slices, i.e. slice 1 will be the first time point
  #xticks, xlabels: ticks/labels to put on the x axis specified as x-axis values
  #        if only one is specified, xticks and xlabels will be the same
  #        if a single number (i) is given, every i-th value will be labelled/ticked
  #        only labels which match to the tick positions will be displayed (ggplot behavior)
  #ylims: lower/upper limit for values on y-axis
  #yticks, ylabels: ticks/labels to put on the y axis specified as y-axis values (ignored if ylims is NULL)
  #        if only one is specified, yticks and ylabels will be the same
  #        if a single number (i) is given, every i-th value will be labelled/ticked
  #        only labels which match to the tick positions will be displayed (ggplot behavior)
  #size: line width (defaults to 1.2) or boxplot width (defaults to .5)
  #signpos: y-position of significance indication (if test argument is supplied, cf. decode.test)
  #txt: text to annotate inside the plot
  #txtpos: position of txt
  #txtcol, txtsize: further txt properties
  #cols: colours according to groupvar or xvar if no groupvar
  #ggTheme: ggplot theme, possible to modify details by adding + theme(...) 
  #         e.g. ggTheme + theme(legend.position = "none")
  #...: further arguments to decode.test
  .eval_package("ggplot2", load=T)
  result = copy(result)
  if (!is.data.table(result)) setDT(result)
  # args.in = lapply( as.list(match.call())[-c(1,2)], eval.parent, n=2 )
  args.in = as.list(match.call())[-c(1,2)]
  for (e in args.in) eval(e) #evaluate variable input within function scope
  if (is.null(size)) size = ifelse(lineplot, 1.2, .5)
  args = do.call( .eval_ellipsis, modifyList(args.in, list(funStr="decode.test")) )
  args$decode.test = lapply( modifyList(args$decode.test, list(dv="y", within.var="g", between="x"))[-1], eval )
  IDs = !is.null(idvar) && nzchar(idvar) && idvar %in% names(result)
  groups = !is.null(groupvar) && nzchar(groupvar) && groupvar %in% names(result)
  ptest = !is.null(args.in$test) && args$decode.test$test != "none"
  average = (average && is.null(multi.layout)) || !IDs #one plot
  if (ptest) { #make sure everything is in order for the significance testing (if requested)
    if (!IDs) stop( "Please supply an id for the group comparison." )
    if (!groups) stop( "Please supply a group variable for which to compute the comparison." )
    if ( length(unique(result[[groupvar]])) != 2 ) stop( "The group variable must have exactly 2 levels." )
    if (!(average || (lineplot && is.null(multi.layout))) && args$decode.test$id == idvar) {
      if ("fold" %in% names(result)) args$decode.test$id = "fold"
      else stop( "Group and test id are identical with no fold information present.", 
                 "Please specifiy with 'id' the identifier for the within case comparison." )
    } 
  }
  #rename (or add if missing) columns
  setnames(result, old=c(yvar,xvar), new=c("y","x"))
  if (groups) setnames(result, old=groupvar, new="g") else result[, g:=0]
  if (IDs) setnames(result, old=idvar, new="id") else result[, id:=0]
  data = result[, .(y=mean(y)), by=.(id,x,g)] #average y per id, group and x
  signtest = NULL
  if (ptest) { #calculate p-values
    if ( average || (lineplot && is.null(multi.layout)) ) {
      signtest = do.call(decode.test, modifyList( args$decode.test, list(result=data, bootstrap=F, id="id")) )
    } else { #evaluate significance for each subject individually
      signtest = result[, do.call(decode.test, modifyList( args$decode.test, list(result=.SD, bootstrap=F)) ), by=id]
    }
  }
  CI = CI && lineplot
  if (CI) { #compute confidence intervals for the means
    data = result[, {
      if (average) {
        m = mean(y) #global mean of all subject data (within x,g)
        if (args$decode.test$bootstrap) {
          se = .SD[, .(y=mean(y)), by=id][, sd( replicate(args$decode.test$iterations, mean(sample(y, replace=T))) )]
        } else { #normalize subject-specific means to be identical  
          se = .SD[, .(yn = y - mean(y) + m), by=id][, sd(yn)/sqrt(.N)]
        }
        .(y = m, se = se) #return
      } else { #no averaging
        if (args$decode.test$bootstrap) {
          .SD[, .(y=mean(y), se=sd( replicate(args$decode.test$iterations, mean(sample(y, replace=T))) )), by=id]
        } else {
          .SD[, .(y=mean(y), se=sd(y)/sqrt(.N)), by=id]
        }
      } 
    }, by=.(x,g)]
#     if (average && !args$decode.test$bootstrap) { #apply correction
#       nWithin = data[, uniqueN(g)]
#       data[, se := se * sqrt(nWithin/(nWithin-1))]
#     }
  } else if (!CI && average && lineplot) {
    data = data[, .(y=mean(y)), by=.(x,g)]
  } else if (!lineplot) {
    data = copy(result) #unaggregated data for individual boxplots
  }
  if (lineplot || groups) { #design xAxis if lineplot or for boxplot if groups were specified in addition to xvar
    npoints = data[, uniqueN(x)] #number of points on xAxis
    if (!is.null(xrange)) { #replace slice number with time stamp
      timestep = diff(xrange)/(npoints-1)
      timepoints = seq(xrange[1], xrange[2], by=timestep) #left aligned
      data[, x := timepoints[ findInterval(x, unique(x)) ]] #add time info: slices are likely unsorted
      if (ptest) signtest[, x := timepoints[ findInterval(x, unique(x)) ]]
    }
    #create labels/ticks if any are non-NULL
    if ( !is.null(c(xlabels, xticks)) ) {
      xVals = data[, unique(x)]
      if (length(xlabels) == 1) xlabels = xVals[ seq(1,npoints,xlabels) ] #step size defined
      if (length(xticks) == 1) xticks = xVals[ seq(1,npoints,xticks) ] #step size defined
      if (is.null(xlabels)) xlabels = xticks
      if (is.null(xticks)) xticks = xlabels
      if (length(xticks) != length(xlabels)) {
        tmp = xticks
        tmp[ !xticks %in% xlabels ] = ""
        xlabels = tmp
      }
    }
  }
  if ( !is.null(ylims) && !is.null(c(ylabels, yticks)) ) { #yAxis
    if (length(ylabels) == 1) ylabels = seq(ylims[1], ylims[2], length.out=ylabels)
    if (length(yticks) == 1) yticks = seq(ylims[1], ylims[2], length.out=yticks)
  }
  if (is.null(ylabels)) ylabels = yticks
  if (is.null(yticks)) yticks = ylabels
  if (length(yticks) != length(ylabels)) {
    tmp = yticks
    tmp[ !yticks %in% ylabels ] = ""
    ylabels = tmp
  }
  if (is.null(signpos)) signpos = ifelse(is.null(ylims), data[, max(y)], ylims[2]) #significance indicator ypos
  ## preparation end ##
  ggLineplot <- function(data, signtest=NULL, txt="") {
    p = ggplot(data=data, aes(x=x, y=y, group=g, fill=g, colour=g)) +
      geom_line(size=size) + ggTheme +
      labs(x=xlab, y=ylab, group=legend, fill=legend, colour=legend)
    if (length(cols) >= data[, uniqueN(g)]) p = p + scale_color_manual(values=cols)
    if (!is.null(ylims)) {
      p = p + coord_cartesian(ylim=ylims)
      if (!is.null(ylabels)) p = p + scale_y_continuous(breaks=yticks, labels=ylabels)
    }
    if (!is.null(xlabels)) p = p + scale_x_continuous(breaks=xticks, labels=xlabels, expand=c(0,0))
    if (!is.null(hline)) p = p + geom_hline(yintercept=hline , linetype="dotted")
    if (!is.null(vline)) p = p + geom_vline(xintercept=vline , linetype="dashed")
    if (nzchar(txt)) p = p + annotate("text", x=txtpos[1], y=txtpos[2], label=txt, col=txtcol, size=txtsize)
    if (CI) p = p + geom_ribbon(aes(ymin=y-t*se, ymax=y+t*se), alpha=0.3, colour=NA, show.legend=F) + scale_fill_manual(values=cols)
    if (ptest) { #draw horizontal line with its center above significant time point(s)
      signtime = signtest[nzchar(sign), x]
      consec = signtest[, rle(nzchar(sign))] #consecutive (non-)significant timepoints
      consec = consec$lengths[consec$values] #subset significant streaks
      xstep = data[, diff(unique(x)[1:2])]/2 #to adjust center point of line
      end = signtime[ cumsum(consec) ]+xstep #offsets of significance
      start = signtime[ cumsum(consec)+1-consec ]-xstep #onsets of significance
      for (i in seq_along(start)) {
        p = p + geom_segment(x=start[i], xend=end[i], y=signpos, yend=signpos, size=size, col=ifelse(average, cols[1], "black"))
      }
    }
    return(p)
  }
  ggBoxplot <- function(data, signtest=NULL, txt="") {
    update_geom_defaults("point", list(colour = NULL)) #for boxplot outlier colouring
    if (groups) {
      data[, grouping:=interaction(x,g)]
      p = ggplot(data=data, aes(x=x, y=y, group=grouping, fill=g, colour=g))
      if (!is.null(vline)) p = p + geom_vline(xintercept=vline , linetype="dashed")
      if (!is.null(xlabels)) p = p + scale_x_continuous(breaks=xticks, labels=xlabels, expand=c(0.01,0.01))
      if (length(cols) >= data[, uniqueN(g)]) p = p + scale_color_manual(values=cols) + scale_fill_manual(values=cols)
    } else { #no groupvar
      p = ggplot(data=data, aes(x=x, y=y, fill=x, colour=x))
      if (length(cols) >= data[, uniqueN(x)]) p = p + scale_color_manual(values=cols) + scale_fill_manual(values=cols)
    }
    if (!is.null(hline)) p = p + geom_hline(yintercept=hline , linetype="dotted")
    p = p + geom_boxplot(size=size, outlier.size=size/2) + labs(x=xlab, y=ylab, fill=legend, colour=legend) + ggTheme
    if (!is.null(ylims)) {
      p = p + coord_cartesian(ylim=ylims)
      if (!is.null(ylabels)) p = p + scale_y_continuous(breaks=yticks, labels=ylabels)
    }
    if (nzchar(txt)) p = p + annotate("text", x=txtpos[1], y=txtpos[2], label=txt, col=txtcol, size=txtsize)
    if (ptest) { #draw horizontal line with its center above significant time point(s)
      signtime = signtest[nzchar(sign), x]
      consec = signtest[, rle(nzchar(sign))] #consecutive (non-)significant timepoints
      consec = consec$lengths[consec$values] #subset significant streaks
      xstep = data[, diff(unique(x)[1:2])]/2 #to adjust center point of line
      end = signtime[ cumsum(consec) ]+xstep #offsets of significance
      start = signtime[ cumsum(consec)+1-consec ]-xstep #onsets of significance
      for (i in seq_along(start)) p = p + geom_segment(x=start[i], xend=end[i], y=signpos, yend=signpos, size=size, col=cols[1])
    }
    #add median as black lines
    d = layer_data(p, i=as.integer(!is.null(hline)) + as.integer(!is.null(vline)) + 1L ) 
    p = p + geom_segment(data=d, aes(x=xmin, xend=xmax, y=middle, yend=middle), colour="black", size=1, inherit.aes=F)
    update_geom_defaults("point", list(colour = "black")) #change back
    return(p)
  }
  if (CI && is.null(args.in$t.ci)) {
    t = qt(1-args$decode.test$alpha/2, result[x==x[1], .N, by=g][1,N-1]) #t-statistic for confidence interval
  } else if (CI) {
    t = args.in$t.ci
  }
  if (!groups && lineplot) ggTheme = ggTheme + theme(legend.position="none")
  data[, g:=as.factor(g)]
  # #undo changes to result
  # setnames(result, old=c("y","x"), new=c(yvar,xvar))
  # if (!groups) result[, g:=NULL] else setnames(result, "g", groupvar)
  # if (!IDs) result[, id:=NULL] else setnames(result, "id", idvar)
  ## plot 
  plotStr = ifelse(lineplot, "ggLineplot", "ggBoxplot")
  if (average) { #single lineplot aggregated for all subjects
    p = do.call(plotStr, list(data, signtest, txt))
    return(p)
  } else { #no averaging
    n = data[, uniqueN(id)]
    if (is.null(multi.layout) && lineplot) { #single plot
      data = data[ g != args$decode.test$within.null ][, g:=NULL] #keep only true label data
      if (length(cols) < n) cols = colorRampPalette(c("gray30","gray90"), alpha=T)(n)
      setnames(data, "id", "g")[, g := as.factor(g)] #group (colour) by id
      p = ggLineplot(data, signtest, txt=txt) + 
        geom_line(mapping=aes(x=x,y=y), size=size, col="black", data=data[, .(y=mean(y)), by=x], inherit.aes=F) #average
      return(p)
    } else { #multiple plots
      if (is.null(multi.layout)) multi.layout=c(1,1) #one boxplot per page
      .eval_package("gridExtra", load=T)
      idVals = data[, unique(id)]
      if (length(txt) != n) txt = paste0(txt, idVals)
      plist = lapply(seq_len(n), function(i) do.call( plotStr, list(data[id == idVals[i]], if (ptest) signtest[id == idVals[i]], txt=txt[i]) ))
      multipage = multi.layout[1] * multi.layout[2] #plots per page
      for (i in seq(1+multipage, n+multipage, multipage)) { #iterate until all plots are created
        pl = lapply( (i-multipage):min(i-2,n-1), function(x) plist[[x]]+theme(legend.position="none") )
        do.call("grid.arrange", c( c(pl, plist[min(i-1,n)]), nrow=multi.layout[1], ncol=multi.layout[2] ))
        # do.call("grid.arrange", c( plist[ (i-multipage):min(i-1,n) ], nrow=multi.layout[1], ncol=multi.layout[2] ))
      }
      # return(plist)
    }
  }
}

plot.GAT <- function(GAT, contrast=F, range=NULL, onset=T, limits=NULL, ticks=NULL, labels=NULL, 
                     dv="AUC", groupvar="label", idvar="subject", xlab="Generalization time",
                     ylab="Training time", legend="AUC", txt="", txtpos=NULL, txtcol="#666666", txtsize=5, 
                     contrast.col="black", showBar=T, ggTheme=.ggTheme(), colScheme=NULL, colNum=100, smooth=F, ...) {
  ## plots the GAT matrix
  #INPUT:
  #GAT: the output of decode with GAT=T, more generally a data table with columns slice, test
  #contrast: if True, p values are computed for the difference in dv between groupvar;
  #          significant pvals will be plotted ontop as contour lines (colour set via contrast.col)
  #          arguments of decode.test (via ...) can be used to modify the test settings
  #range: actual time range of measurements, i.e time of 1st/last sample
  #       the timepoints are left-aligned, i.e. first slice = range[1]
  #onset: if True, highlight time 0 with dashed lines
  #limits: limits of the colour scale (corresponds to dv)
  #ticks: ticks to put on the x/y axis; if only labels are specified, ticks will be set where labels are
  #labels: time stamps to display, defaults to every tick
  #groupvar: group for contrast, can be left empty if no contrast;
  #          any group factor that is different from within.null will be aggregated for the plot
  #idvar: id variable for contrast, can be left empty if no contrast
  #xlab, ylab, legend: titles of the respective part
  #txt, txtpos, txtcol, txtsize: text to display within the plot and respective settings
  #contrast.col: colour of the contour lines if contrast = T and significant pvals found
  #showBar: show the colour bar to indicate colour coding of dv; 
  #         these settings can be modified: barwidth, barheight, nbin, barticks
  #colScheme: colour scheme of the raster, defaults to matlab style (jetcolours)
  #colNum: number with which a colour gradient is created from colScheme, i.e. number of distinct colours
  #smooth: if True, the raster is smoothed via linear interpolation
  #...: further arguments to decode.test (if contrast)
  if ( !all(c("slice","test") %in% names(GAT)) ) stop( "GAT is expected with columns 'slice' and 'test'." )
  .eval_package("ggplot2", load=T)
  args.in = as.list(match.call())[-c(1,2)]
  for (e in args.in) eval(e) #evaluate variable input within function scope
  args = lapply( do.call(.eval_ellipsis, modifyList(args.in, list(funStr="decode.test")))$decode.test[-1], eval )
  if (is.null(colScheme)) colScheme = c("#00007F","blue","#007FFF","cyan","#7FFF7F","yellow","#FF7F00","red","#7F0000")
  IDs = !is.null(idvar) && nzchar(idvar) && idvar %in% names(GAT)
  groups = !is.null(groupvar) && nzchar(groupvar) && groupvar %in% names(GAT)
  if (contrast) {
    if (!IDs) stop( "Please supply an id for the contrast." )
    if (!groups) stop( "Please supply a group variable for which to compute the contrast." )
    if ( length(unique(GAT[[groupvar]])) != 2 ) stop( "The group variable must have exactly 2 levels." )
    if ( !args$within.null %in% unique(GAT[[groupvar]]) ) stop( "Please specify the null variable for the contrast via 'within.null'. ")
    setnames(GAT, old=idvar, new="id")
  }
  setnames(GAT, old=dv, new="dv")
  #aggregate data & compute contrast between groups if applicable
  if (groups) {
    setnames(GAT, old=groupvar, new="g")
    if (contrast) { #test at each slice,test pair the group contrast for id
      GAT.data = GAT[, .(dv=mean(dv)), by=.(slice,test,g,id)]
      args$test = tolower(args$test)
      if (is.null(args.in$one.sided)) alt = "two.sided" else alt = ifelse(one.sided, "greater", "two.sided")
      #transform for speedup
      N = GAT.data[, uniqueN(slice)] #N for both slice/test
      dvar = GAT.data$dv #vector for fast subsetting
      gidx = GAT.data[.(1,1)][g!=args$within.null, which=T] #group 1 idx within a slice/test time, rest is group 2
      sidx = GAT.data[.(1), which=T] #a single training time
      ptest = lapply(seq_len(N), function(s) {
        slice = sidx+length(sidx)*(s-1) #current training time indices
        tidx = split(slice, GAT.data[slice,test]) #all test times
        sapply(tidx, function(i) {
          if (args$test == "wilcox") wilcox.test(x=dvar[i][gidx], y=dvar[i][-gidx], paired=args$paired, alternative=alt)$p.value
          else t.test(x=dvar[i][gidx], y=dvar[i][-gidx], paired=args$paired, alternative=alt)$p.value
        })
      })
      #remove null group, aggregate id and merge with ptest
      GAT.data = GAT.data[g != args$within.null, .(dv=mean(dv)), by=.(slice,test)][, p:=p.adjust(unlist(ptest), method=args$adjust)]
    } else { #no contrast but remove null group and aggregate
      GAT.data = GAT[g != args$within.null, .(dv=mean(dv)), by=.(slice,test)]
    }
  } else { #no grouping / contrast, id does not matter
    GAT.data = GAT[, .(dv=mean(dv)), by=.(slice,test)]
  }
  #undo changes to GAT
  setnames(GAT, old="dv", new=dv)
  if ("id" %in% names(GAT)) setnames(GAT, "id", idvar)
  if ("g" %in% names(GAT)) setnames(GAT, "g", groupvar)
  #create timeAxis
  npoints = GAT.data[, uniqueN(slice)] #number of points on axis
  if (!is.null(range)) { #replace slice number with time stamp
    timestep = diff(range)/(npoints-1)
    timepoints = seq(range[1], range[2], by=timestep) #left aligned
    GAT.data[, ':=' (slice = timepoints[ findInterval(slice, unique(slice)) ], 
                     test = timepoints[ findInterval(test, unique(test)) ])]
  }
  #create labels/ticks if any are non-NULL
  if ( !is.null(c(labels, ticks)) ) {
    tVals = GAT.data[, unique(slice)]
    if (length(labels) == 1) labels = tVals[ seq(1,npoints,labels) ] #step size defined
    if (length(ticks) == 1) ticks = tVals[ seq(1,npoints,ticks) ] #step size defined
    if (is.null(labels)) labels = ticks
    if (is.null(ticks)) ticks = labels
    if (length(ticks) != length(labels)) {
      tmp = ticks
      tmp[ !ticks %in% labels ] = ""
      labels = tmp
    }
  }
  #set colour bar settings
  if (showBar) {
    if ( !is.character(ggTheme$legend.position) || ggTheme$legend.position == "none" ) {
      ggTheme = ggTheme + theme(legend.position="right")
    }
    if ( is.null(args.in$barwidth) ) args.in$barwidth = .7
    if ( is.null(args.in$barheight) ) args.in$barheight = 20
    if ( is.null(args.in$nbin) ) args.in$nbin = colNum
    if ( is.null(args.in$barticks) ) args.in$barticks = F
  } else { #don't show colorbar
    ggTheme = ggTheme + theme(legend.position="none") 
  }
  p = ggplot(data=GAT.data, aes(x=test, y=slice)) + 
    geom_raster(aes(fill=dv), interpolate=smooth) + 
    annotate("text", x=txtpos[1], y=txtpos[2], label=txt, col=txtcol, size=txtsize) +
    labs(x=xlab, y=ylab, fill=legend) + ggTheme
  if (!is.null(limits)) {
    p = p + scale_fill_gradientn(colours=colorRampPalette(colScheme)(colNum), limits=c(limits[1],limits[2]),
                                 guide=guide_colorbar(bardwidth=args.in$barwidth, barheight=args.in$barheight, 
                                                      ticks=eval(args.in$barticks), nbin=args.in$nbin))
  } else { #auto limits
    p = p + scale_fill_gradientn(colours=colorRampPalette(colScheme)(colNum),
                                 guide=guide_colorbar(bardwidth=args.in$barwidth, barheight=args.in$barheight, 
                                                      ticks=eval(args.in$barticks), nbin=args.in$nbin))
  }
  if ( !is.null(c(labels, ticks)) ) {
    p = p + scale_y_continuous(breaks=ticks, labels=labels, expand=c(0,0)) +
      scale_x_continuous(breaks=ticks, labels=labels, expand=c(0,0))
  } else {
    p = p + scale_y_continuous(expand=c(0,0)) + scale_x_continuous(expand=c(0,0))
  }
  if (contrast) p = p + stat_contour(aes(z=p), col=contrast.col, breaks=args$alpha)
  if ( onset && GAT.data[, min(slice) < 0] )  {
    p = p + geom_vline(xintercept=0 , linetype="dashed", size=0.1) +
      geom_hline(yintercept=0 , linetype="dashed", size=0.1)
  }
  return(p)
}


plot.ERP <- function(data, multi.layout=NULL, butterfly=F, columns=NULL, colours=NULL, 
                     xrange=NULL, xticks=NULL, xlabels=NULL, ylims=NULL, yticks=NULL, ylabels=NULL, 
                     xlab="Time (ms)", ylab="Amplitude", legend="Class", vline=NULL, size=1.2,
                     txt="", txtpos=NULL, txtcol="#666666", txtsize=5, ggTheme=.ggTheme()) {
  ## ERP: plot amplitudes at each time point aggregated over trials (separately for outcome)
  #multi.layout: if NULL, aggregate over subjects, else a vector specifying nrow, ncol
  #butterfly: if True, plot each channel, else aggregate over channels
  #columns: measurement channels to average over, if NULL non-key columns
  #colours: if butterfly, colour for each channel, else colour for each outcome
  #rest: see plot.decode
  .eval_package("ggplot2", load=T)
  # data = data.check(data)
  if (is.null(columns)) columns = setdiff( names(data), key(data) )
  n = ifelse(butterfly, length(columns), data[, uniqueN(outcome)])
  if (is.null(colours) && n <= 5) colours = c("#0072BD", "#D95319", "#EDB120", "#7E2F8E", "#77AC30")[1:n]
  if (butterfly && legend == "Class") legend = "Channel"
  #aggregate ERP for each channel:
  if (is.null(multi.layout)) { #aggregate over subjects
    ERP = data[, lapply(.SD, mean), .SDcols=columns, by=.(sample,outcome)]
    if (!butterfly) { #aggregate the channels
      ERP = ERP[, .(erp=rowMeans(.SD)), by=.(sample,outcome)]
    }
  } else {
    ERP = data[, lapply(.SD, mean), .SDcols=columns, by=.(subject,sample,outcome)]
    if (!butterfly) { #aggregate the channels
      ERP = ERP[, .(erp=rowMeans(.SD)), by=.(subject,sample,outcome)]
    }
  }
  if (butterfly) { #each channel individually
    ERP = melt(ERP, measure.vars=columns, value.name="erp", variable.name="channel")
  }
  npoints = ERP[, uniqueN(sample)] #number of points on xAxis
  if (!is.null(xrange)) { #replace slice number with time stamp
    timestep = diff(xrange)/npoints
    timepoints = seq(xrange[1], xrange[2], by=timestep)[-(npoints+1)] #left aligned
    ERP[, sample := timepoints[ findInterval(sample, unique(sample)) ]] #add time info
  }
  #create labels/ticks if any are non-NULL
  if ( !is.null(c(xlabels, xticks)) ) {
    xVals = ERP[, unique(sample)]
    if (length(xlabels) == 1) xlabels = xVals[ seq(1,npoints,xlabels) ] #step size defined
    if (length(xticks) == 1) xticks = xVals[ seq(1,npoints,xticks) ] #step size defined
    if (is.null(xlabels)) xlabels = xticks
    if (is.null(xticks)) xticks = xlabels
    if (length(xticks) != length(xlabels)) {
      tmp = xticks
      tmp[ !xticks %in% xlabels ] = ""
      xlabels = tmp
    }
  }
  if ( !is.null(ylims) && !is.null(c(ylabels, yticks)) ) { #yAxis
    if (length(ylabels) == 1) ylabels = seq(ylims[1], ylims[2], length.out=ylabels)
    if (length(yticks) == 1) yticks = seq(ylims[1], ylims[2], length.out=yticks)
  }
  if (is.null(ylabels)) ylabels = yticks
  if (is.null(yticks)) yticks = ylabels
  if (length(yticks) != length(ylabels)) {
    tmp = yticks
    tmp[ !yticks %in% ylabels ] = ""
    ylabels = tmp
  }
  ggLineplot <- function(dat, txt) {
    if (butterfly) p = ggplot(data=dat, aes(x=sample, y=erp, group=channel, colour=channel)) + facet_grid(.~outcome)
    else p = ggplot(data=dat, aes(x=sample, y=erp, group=outcome, colour=outcome))
    p = p + geom_line(size=size) + labs(x=xlab, y=ylab, group=legend, colour=legend) + ggTheme
    if (!is.null(colours)) p = p + scale_color_manual(values=colours)
    if (!is.null(ylims)) {
      p = p + coord_cartesian(ylim=ylims)
      if (!is.null(ylabels)) p = p + scale_y_continuous(breaks=yticks, labels=ylabels)
    }
    if (!is.null(xlabels)) p = p + scale_x_continuous(breaks=xticks, labels=xlabels, expand=c(0,0))
    if (!is.null(vline)) p = p + geom_vline(xintercept=vline , linetype="dashed", size=0.1)
    if (nzchar(txt)) p = p + annotate("text", x=txtpos[1], y=txtpos[2], label=txt, col=txtcol, size=txtsize)
    return(p)
  }
  if (is.null(multi.layout)) {
    return( ggLineplot(ERP, txt) )
  } #else multi.layout specified
  .eval_package("gridExtra", load=T)
  subIDs = ERP[, unique(subject)]
  if (length(txt) != length(subIDs)) txt = paste0(txt, subIDs)
  plist = lapply(seq_along(subIDs), function(i) ggLineplot(ERP[subject == subIDs[i]], txt[i]) )
  multipage = multi.layout[1] * multi.layout[2] #plots per page
  for (i in seq(1+multipage, length(subIDs)+multipage, multipage)) { #iterate until all plots are created
    pl = lapply( (i-multipage):min(i-2,length(subIDs)-1), function(x) plist[[x]]+theme(legend.position="none") )
    do.call("grid.arrange", c( c(pl, plist[min(i-1,length(subIDs))]), nrow=multi.layout[1], ncol=multi.layout[2] ))
    # do.call("grid.arrange", c( plist[ (i-multipage):min(i-1,length(subIDs)) ], nrow=multi.layout[1], ncol=multi.layout[2] ))
  }
}




