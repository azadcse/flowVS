## ====================================================================================
## Identifying lymphocytes from a flowFrame 
## (note: there is a function in flowStats package called lymphGate performing similar task)
## arguments:
##      ff: a flowFrame from which lympocytes are identified 
##      lymph.boundary: a list denoting a rough bundary for lymphocites 
##                  1st element represents the lower and upper limit of fsc
##                  2nd element represents the lower and upper limit of ssc
##                  example: lymph.boundary <- list("FS"=c(180000, 500000),"SS"=c(0, 180000))
##      fsc/ssc: name (or numeric index) of the forward and side scatter channels
##      plot: if true then plots the rectangular and elliptical gates for the lymphocyes
## Return: A flowFrame only with lymphocytes
## ====================================================================================

lymphs <- function(ff, lymph.boundary, fsc, ssc, plot=FALSE)
{    
  allcols <- colnames(ff)
  if(is.numeric(fsc)) fsc <- allcols[fsc]
  if(is.numeric(ssc)) ssc <- allcols[ssc]
  if(is.na(match(fsc,allcols)))
    stop("Could not find forward scatter parameter. ",
         "Please set the fsc parameter", call.=FALSE)
  if(is.na(match(ssc,allcols)))
    stop("Could not find side scatter parameter. ",
         "Please set the ssc parameter", call.=FALSE)
  names(lymph.boundary) = c(fsc,ssc)
  rect.gate <- rectangleGate(filterId="lymphs", .gate=lymph.boundary )
  ff.rect <- Subset(ff, rect.gate)
  
  n2f <- norm2Filter(fsc, ssc, scale.factor=2)
  ff.lymph = Subset(ff.rect, n2f)
  
  if(plot==TRUE)
  {
    ylim = c(min(exprs(ff[,fsc])), max(exprs(ff[,fsc])))
    xlim = c(min(exprs(ff[,ssc])), max(exprs(ff[,ssc])))
    plot.formula = as.formula(sprintf("`%s` ~ `%s`",fsc, ssc))
    p = xyplot(plot.formula, main='(1) Approx. rectangular gate', nbin=128, nrpoints = 100, data=ff, filter = rect.gate, smooth=FALSE, colramp = colorRampPalette(c("blue", "green", "yellow", "red")), xlim=xlim, ylim=ylim)
    print(p,  position = c(0.0, 0, 0.5, 1), more=TRUE)
    ylim = c(min(exprs(ff.rect[,fsc])), max(exprs(ff.rect[,fsc])))
    xlim = c(min(exprs(ff.rect[,ssc])), max(exprs(ff.rect[,ssc])))
    p = xyplot(plot.formula, main='(2) Elliptical gate for lymphocytes', nbin=128, nrpoints = 100, data=ff.rect, filter = n2f, smooth=FALSE, colramp = colorRampPalette(c("blue", "green", "yellow", "red", "red")), xlim=xlim, ylim=ylim)
    print(p, position = c(0.5, 0, 1, 1))
  }
  return(ff.lymph)
}