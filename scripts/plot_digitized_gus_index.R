out <- 'results/summary/pident95_cov90/gus_index_digitized'
coords <- read.delim(file.path(out,'digitized_activity.tsv'))
p <- read.delim(file.path(out,'subject_predictions.tsv'))
d <- read.delim('results/summary/pident95_cov90/assay_order_stacked_rpkm/input_abundance.tsv',row.names=1,check.names=FALSE)
X <- t(as.matrix(d[,coords$subject])); y <- p$estimated_activity
full <- list(pred=matrix(p$in_sample_prediction,ncol=1)); j <- 1
cv <- list(pred=p$nested_LOO_prediction)
fig<-function(){
 par(mfrow=c(1,3),oma=c(1.5,0,0,0),mar=c(4.6,4.3,3,1.0),mgp=c(2.8,.7,0),family='sans')
 panels<-list(list(x=rowSums(X),title='Unweighted total abundance',xl='Total GUS abundance (RPKM)',color='#64748B'),list(x=full$pred[,j],title='Weighted index: fitted data',xl='Fitted activity (nM)',color='#C17C36'),list(x=cv$pred,title='Weighted index: held-out subjects',xl='Nested LOO predicted activity (nM)',color='#247C82'))
 for(p in panels){
  plot(p$x,y,pch=21,bg=p$color,col='white',cex=1.35,xlab=p$xl,ylab='Estimated 15-min activity (nM)',main=p$title,cex.main=.95,ylim=c(0,25),bty='l')
  if(p$title!='Unweighted total abundance')abline(0,1,col='#C3C8CC',lty=2)
  selids<-which(coords$subject %in% c('S11','S10','S14','S21'));text(p$x[selids],y[selids],coords$subject[selids],pos=3,cex=.65)
  mtext(sprintf('Pearson r = %.2f; Spearman rho = %.2f',cor(p$x,y),cor(p$x,y,method='spearman')),side=3,line=.25,cex=.7)
 }
 mtext('Exploratory only: activity digitized from a figure; weights fitted to 17 subjects',outer=TRUE,side=1,line=0,cex=.8)
}
png(file.path(out,'gus_index_exploration.png'),width=3300,height=1200,res=250,type='cairo');fig();dev.off()
cairo_pdf(file.path(out,'gus_index_exploration.pdf'),width=13.2,height=4.8);fig();dev.off()
