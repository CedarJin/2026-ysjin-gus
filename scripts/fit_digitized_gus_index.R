#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(glmnet))
out <- 'results/summary/pident95_cov90/gus_index_digitized'
coords <- read.delim(file.path(out,'digitization_coordinates.tsv'))
coords$estimated_activity_nM <- round((237.5-coords$estimated_mean_y_px)*60/(237.5-58.5),2)
coords$assumed_reading_halfwidth_nM <- 0.5
write.table(coords,file.path(out,'digitized_activity.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
d <- read.delim('results/summary/pident95_cov90/assay_order_stacked_rpkm/input_abundance.tsv',row.names=1,check.names=FALSE)
X <- t(as.matrix(d[,coords$subject])); y <- coords$estimated_activity_nM
stopifnot(nrow(X)==17,all(is.finite(X)),all(is.finite(y)))
lambdas <- 10^seq(3,-6,length.out=19)
path <- function(x,y,test,positive=TRUE){
 mu<-colMeans(x); s<-apply(x,2,sd);keep<-s>1e-10
 ym<-mean(y);ys<-sd(y)
 xx<-sweep(sweep(x[,keep,drop=FALSE],2,mu[keep]),2,s[keep],'/')
 tt<-sweep(sweep(test[,keep,drop=FALSE],2,mu[keep]),2,s[keep],'/')
 f<-glmnet(xx,(y-ym)/ys,alpha=0,lambda=lambdas,standardize=FALSE,intercept=TRUE,lower.limits=if(positive)0 else -Inf,thresh=1e-10,maxit=100000)
 predictions<-as.matrix(predict(f,newx=tt,s=lambdas))*ys+ym
 cc<-as.matrix(coef(f,s=lambdas));w<-matrix(0,ncol(x),length(lambdas),dimnames=list(colnames(x),NULL));w[keep,]<-sweep(cc[-1,,drop=FALSE],1,s[keep],'/')*ys
 intercept<-ym+ys*cc[1,]-as.numeric(mu%*%w)
 list(pred=predictions,w=w,intercept=intercept)
}
select_lambda <- function(x,y){
 pred<-matrix(NA,length(y),length(lambdas))
 for(i in seq_along(y))pred[i,]<-path(x[-i,,drop=FALSE],y[-i],x[i,,drop=FALSE])$pred
 loss<-colMeans((pred-y)^2)
 # In ties, prefer stronger regularization (lambdas decreasing).
 list(j=which.min(loss),loss=loss)
}
nested <- function(x,y,details=FALSE){
 n<-length(y);pred<-numeric(n);lam<-numeric(n);weights<-matrix(0,n,ncol(x));base<-numeric(n);total<-numeric(n)
 for(i in seq_len(n)){
  tr<-setdiff(seq_len(n),i);sel<-select_lambda(x[tr,,drop=FALSE],y[tr]);fit<-path(x[tr,,drop=FALSE],y[tr],x[i,,drop=FALSE]);pred[i]<-fit$pred[1,sel$j];lam[i]<-lambdas[sel$j];weights[i,]<-fit$w[,sel$j];base[i]<-mean(y[tr])
  a<-rowSums(x[tr,,drop=FALSE]);b<-lm(y[tr]~a);total[i]<-as.numeric(coef(b)[1]+coef(b)[2]*sum(x[i,]))
 }
 list(pred=pred,lambda=lam,weights=weights,baseline=base,total=total,improvement=mean((y-base)^2)-mean((y-pred)^2))
}
metrics <- function(name,p,y){data.frame(model=name,pearson=cor(p,y),spearman=cor(p,y,method='spearman'),RMSE=sqrt(mean((p-y)^2)),MAE=mean(abs(p-y)),Q2_vs_global_mean=1-sum((p-y)^2)/sum((y-mean(y))^2))}
cat('Starting nested leave-one-subject-out cross-validation\n');flush.console()
cv<-nested(X,y)
sel<-select_lambda(X,y);full<-path(X,y,X);j<-sel$j
ols<-lm(y~X);bestnonneg<-full$pred[,length(lambdas)]
summary<-rbind(metrics('Unweighted total RPKM (in-sample linear fit)',fitted(lm(y~rowSums(X))),y),metrics('Unconstrained OLS (in-sample upper-bound diagnostic)',fitted(ols),y),metrics('Near-unpenalized nonnegative fit (in-sample)',bestnonneg,y),metrics('Selected nonnegative ridge (in-sample)',full$pred[,j],y),metrics('Nested LOO nonnegative ridge (held-out)',cv$pred,y),metrics('LOO total RPKM linear model (held-out)',cv$total,y),metrics('LOO training-mean baseline',cv$baseline,y))
write.table(summary,file.path(out,'model_performance.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
w<-full$w[,j];wt<-data.frame(GUS=colnames(X),weight_raw_RPKM=w,weight_sum1=if(sum(w)>0)w/sum(w) else 0,near_unpenalized_nonnegative=full$w[,length(lambdas)],OLS_signed=coef(ols)[-1],outer_fold_min=apply(cv$weights,2,min),outer_fold_max=apply(cv$weights,2,max))
write.table(wt,file.path(out,'gus_weights.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
preds<-data.frame(subject=coords$subject,estimated_activity=y,total_RPKM=rowSums(X),index=as.numeric(X%*%w),in_sample_prediction=full$pred[,j],nested_LOO_prediction=cv$pred,LOO_total_prediction=cv$total,LOO_mean_prediction=cv$baseline,outer_lambda=cv$lambda)
write.table(preds,file.path(out,'subject_predictions.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
write.table(data.frame(lambda=lambdas,inner_LOO_MSE=sel$loss),file.path(out,'full_data_lambda_selection.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
print(summary);print(wt);cat('Starting 99 full-pipeline permutations and 30 reading-error sensitivity runs\n');flush.console()
set.seed(20261015);permy<-replicate(99,sample(y),simplify=FALSE)
workers<-min(2,as.integer(Sys.getenv('SLURM_CPUS_PER_TASK','1')))
perm<-unlist(parallel::mclapply(permy,function(yp)nested(X,yp)$improvement,mc.cores=workers))
pval<-(1+sum(perm>=cv$improvement))/(1+length(perm))
write.table(data.frame(iteration=seq_along(perm),MSE_improvement_over_LOO_mean=perm),file.path(out,'permutation_results.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
set.seed(20261016);jittery<-replicate(30,pmax(0,y+runif(length(y),-.5,.5)),simplify=FALSE)
jit<-parallel::mclapply(jittery,function(yp){z<-nested(X,yp);c(pearson=cor(z$pred,yp),spearman=cor(z$pred,yp,method='spearman'),MSE_improvement=z$improvement)},mc.cores=workers)
write.table(as.data.frame(do.call(rbind,jit)),file.path(out,'reading_error_sensitivity.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
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
notes<-c('EXPLORATORY ANALYSIS USING APPROXIMATE FIGURE-DIGITIZED ACTIVITY',
'Input activity: visually estimated mean-bar y coordinates in original 1650x963 image; zero=237.5 px, 60 nM=58.5 px. Markers obscure exact bar tops. Coordinates and estimated means saved. Readout is PGE2 formed after 15 min, not initial velocity. No replicate reconstruction or control subtraction. +/-0.5 nM is an assumed sensitivity interval, not a measured confidence interval.',
'Predictors: the supplied 8-GUS RPKM matrix, 17 assay subjects. Primary model: nonnegative ridge on training-set standardized predictors and response; centering/scaling and zero-variance feature removal repeated inside every training fold. Lambda grid 10^seq(3,-6,length=19). Inner LOOCV selects minimum MSE; outer LOOCV estimates predictions. No label-based feature selection or model switching. Full-data weights are for exploration and are not independent validation.',
'OLS is an in-sample diagnostic with signed coefficients, not a biological activity index; it maximizes in-sample Pearson correlation among unconstrained linear combinations (with intercept). Neither fitted correlation nor individual coefficients are given inferential P values.',
'Permutation test: 99 permutations of subject activity labels, repeating the entire nested model-fitting procedure; statistic is MSE improvement over the training-mean LOO baseline. Pooled LOOCV correlations can be misleading, especially for nearly constant predictions; use errors and comparison to baseline. No independent external validation.',
'30 sensitivity runs jitter each readout independently by Uniform(-0.5,+0.5) nM; this does not capture systematic reading errors or biological assay uncertainty.',
paste('Selected lambda=',lambdas[j],'intercept=',full$intercept[j]),paste('MSE improvement over LOO mean=',cv$improvement,'Permutation P=',pval),capture.output(print(summary)),capture.output(print(wt)),capture.output(print(apply(do.call(rbind,jit),2,range))))
writeLines(notes,file.path(out,'analysis_notes.txt'));cat(paste(notes,collapse='\n'));writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
