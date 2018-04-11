## fonction ###########################################################################
plot_MLMM.Ago<-function(x) {
	output1 <-x$out[order(x$out$Pos),]
	output_ok<-output1[order(output1$Chr),]

	maxpos<-c(0,cumsum(as.numeric(aggregate(output_ok$Pos,list(output_ok$Chr),max)$x+max(cumsum(as.numeric(aggregate(output_ok$Pos,list(output_ok$Chr),max)$x)))/200)))
	plot_col<-rep(c('gray10','gray60'),ceiling(max(unique(output_ok$Chr))/2))
#	plot_col<-c('blue','darkgreen','red','cyan','purple')
	size<-aggregate(output_ok$Pos,list(output_ok$Chr),length)$x
	a<-rep(maxpos[1],size[1])
	b<-rep(plot_col[1],size[1])
		for (i in 2:max(unique(output_ok$Chr))){
	a<-c(a,rep(maxpos[i],size[i]))
	b<-c(b,rep(plot_col[i],size[i]))}

	output_ok$xpos<-output_ok$Pos+a
	output_ok$col<-b
	output_ok$col[output_ok$mk=='qtl']<-'cyan'
	output_ok$col[output_ok$SNP %in% x$cof]<-'red'

	d<-(aggregate(output_ok$xpos,list(output_ok$Chr),min)$x+aggregate(output_ok$xpos,list(output_ok$Chr),max)$x)/2

	plot(output_ok$xpos,-log10(output_ok$pval),col=output_ok$col,pch=20,ylab='-log10(pval)',xaxt='n',xlab='chromosome')
	axis(1,tick=FALSE,at=d,labels=c(1:max(unique(output_ok$Chr))))
  xline(output_ok$xpos[output_ok$mk=='qtl'], col='cyan', lwd=0.1)
	abline(h=x$bonf_thresh,lty=3,col='black')


    if (length(output_ok$pval[-log10(output_ok$pval) > x$bonf_thresh]) > 0) {
      text(output_ok$xpos[-log10(output_ok$pval) > x$bonf_thresh], -log10(output_ok$pval[-log10(output_ok$pval) > x$bonf_thresh]), output_ok$SNP[-log10(output_ok$pval) > x$bonf_thresh], pos=3, cex=0.7)
      legend("topright", lty=3, paste("bonf thresh :", x$bonf_thresh ,sep=" "))
      } else {
      legend("topright", lty=3, paste("bonf thresh :", x$bonf_thresh ,sep=" "))
      }

}
