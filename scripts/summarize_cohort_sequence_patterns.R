#!/usr/bin/env Rscript
out <- 'results/summary/pident95_cov90/cohort_sequence_audit'
g <- read.delim(file.path(out,'gene_reference_comparisons.tsv'),check.names=FALSE)
s <- read.delim(file.path(out,'subject_sequence_summary.tsv'),check.names=FALSE)
v <- read.delim(file.path(out,'variable_sites_long.tsv'),check.names=FALSE)
h <- read.delim(file.path(out,'protein_haplotypes.tsv'),check.names=FALSE)
stopifnot(nrow(s)==17,!anyDuplicated(s$subject),!anyDuplicated(paste(g$subject,g$reference)))
write_tsv<-function(x,name)write.table(x,file.path(out,name),sep='\t',quote=FALSE,row.names=FALSE,na='NA')
tests <- list()
add_test <- function(family,reference,feature,n,effect,statistic,p,details=''){
 tests[[length(tests)+1]] <<- data.frame(family=family,reference=reference,feature=feature,n=n,effect=effect,statistic=statistic,p=p,details=details)
}
for(trait in c('abundance_weighted_substitutions_per100aa','exact_native_RPKM_fraction','terminal_deficit_ge5_RPKM_fraction')){
 ct<-cor.test(s[[trait]],s$estimated_activity_nM,method='spearman',exact=FALSE)
 no11<-s$subject!='S11'
 add_test('Across-subject sequence burden','all',trait,17,unname(ct$estimate),'Spearman rho',ct$p.value,paste('Without S11 rho=',round(cor(s[[trait]][no11],s$estimated_activity_nM[no11],method='spearman'),3)))
}
# Within-reference mutation burden, adjusting only for that homolog's DNA abundance.
for(ref in unique(g$reference)){
 z<-g[g$reference==ref,]
 if(nrow(z)>=6 && length(unique(z$native_substitutions))>=2){
  m<-lm(estimated_activity_nM~log1p(gene_rpkm)+native_substitutions,z);co<-coef(summary(m))
  add_test('Within-reference substitution burden',ref,'native_substitutions',nrow(z),co['native_substitutions','Estimate'],'Adjusted linear coefficient',co['native_substitutions','Pr(>|t|)'],'Adjusts log1p(same-reference RPKM); not other enzyme expression/activity.')
 }
}
patterns<-list();seen<-character()
for(ref in unique(v$reference)){
 zz<-v[v$reference==ref,]
 for(pos in unique(zz$native_position)){
  z<-zz[zz$native_position==pos,];native<-unique(z$native_residue);stopifnot(length(native)==1)
  if(any(z$sample_residue %in% c('-','X')))next
  for(alt in setdiff(unique(z$sample_residue),native)){
   carriers<-sort(z$subject[z$sample_residue==alt]);others<-sort(z$subject[z$sample_residue!=alt])
   if(min(length(carriers),length(others))<3)next
   key<-paste(ref,paste(carriers,collapse=','),sep=':');feature<-paste0(native,pos,alt)
   if(key %in% seen){patterns[[key]]$features<-paste(patterns[[key]]$features,feature,sep=';');next}
   seen<-c(seen,key);z$group<-as.integer(z$subject %in% carriers)
   m<-lm(estimated_activity_nM~log1p(gene_rpkm)+group,z);co<-coef(summary(m))
   patterns[[key]]<-data.frame(reference=ref,features=feature,n=nrow(z),carriers=paste(carriers,collapse=','),comparison_subjects=paste(others,collapse=','),carrier_n=length(carriers),other_n=length(others),adjusted_effect_nM=co['group','Estimate'],p=co['group','Pr(>|t|)'])
  }
 }
}
if(length(patterns)){
 pp<-do.call(rbind,patterns);rownames(pp)<-NULL
 for(i in seq_len(nrow(pp)))add_test('Within-reference common variant pattern',pp$reference[i],pp$features[i],pp$n[i],pp$adjusted_effect_nM[i],'Adjusted carrier coefficient',pp$p[i],paste('Carriers:',pp$carriers[i],'; others:',pp$comparison_subjects[i]))
 write_tsv(pp,'common_variant_patterns.tsv')
}
alltests<-do.call(rbind,tests);alltests$BH_q<-p.adjust(alltests$p,method='BH');write_tsv(alltests,'exploratory_sequence_associations.tsv')
# Summarize identical amino-acid sequences across people without treating gene rows as independent people.
alleles<-lapply(split(h,h$haplotype),function(z){
 data.frame(reference=z$reference[1],haplotype=z$haplotype[1],n=length(unique(z$subject)),subjects=paste(z$subject,collapse=','),activity_min=min(z$estimated_activity_nM),activity_max=max(z$estimated_activity_nM),activity_fold=max(z$estimated_activity_nM)/min(z$estimated_activity_nM),RPKM_min=min(z$gene_rpkm),RPKM_max=max(z$gene_rpkm),subjects_activity_RPKM=paste(sprintf('%s:%.2fnM/%.3fRPKM',z$subject,z$estimated_activity_nM,z$gene_rpkm),collapse='; '))
})
a<-do.call(rbind,alleles);a<-a[order(-a$n,-a$activity_fold),];write_tsv(a,'identical_sequence_groups_summary.tsv')
# Complete 17 x 8 descriptive matrix; absent genes remain distinct from zero substitutions.
refs<-c('WP_012740861','6NCY_A','WP_022501306','6MVG_A','WP_270534885','WP_007841259','WP_217761416','3CMG_A')
labels<-c('E. eligens','F. saccharivorans','Eubacterium','R. gnavus','B. mass','B. dorei','P. copri','B. fragilis')
subjects<-s$subject
mat<-matrix(NA_real_,length(refs),length(subjects),dimnames=list(refs,subjects));cell<-matrix('',length(refs),length(subjects))
for(i in seq_len(nrow(g))){r<-match(g$reference[i],refs);c<-match(g$subject[i],subjects);mat[r,c]<-g$native_substitutions[i];cell[r,c]<-as.character(g$native_substitutions[i]);if(g$native_N_terminal_region_missing[i]+g$native_C_terminal_region_missing[i]>=5)cell[r,c]<-paste0(cell[r,c],'*')}
plot_summary <- function(){
 layout(matrix(c(1,2,3),ncol=1),heights=c(1,1,2.0));par(oma=c(2.2,0,0,0),family='sans',mgp=c(2.7,.6,0))
 par(mar=c(.6,8,2,1));plot(1:17,s$estimated_activity_nM,type='n',xlim=c(.5,17.5),xaxs='i',axes=FALSE,xlab='',ylab='Activity (nM)',ylim=c(0,25),main='17 subjects: activity, abundance and protein-sequence differences');rect((1:17)-.35,0,(1:17)+.35,s$estimated_activity_nM,col='#438C91',border=NA);axis(2,las=1);box(bty='l')
 par(mar=c(.8,8,.6,1));plot(1:17,s$total_RPKM,type='n',xlim=c(.5,17.5),xaxs='i',axes=FALSE,xlab='',ylab='Total RPKM',ylim=c(0,50));rect((1:17)-.35,0,(1:17)+.35,s$total_RPKM,col='#8793A6',border=NA);axis(2,las=1);box(bty='l')
 par(mar=c(3.2,8.0,.6,1));plot.new();plot.window(xlim=c(.5,17.5),ylim=c(.5,8.5),xaxs='i',yaxs='i')
 pal<-colorRampPalette(c('#FFF8ED','#FBCB8B','#DA7A35','#873817'))(max(mat,na.rm=TRUE)+1)
 for(i in seq_len(8))for(j in seq_len(17)){yy<-9-i;bg<-if(is.na(mat[i,j]))'#EAECEF' else pal[mat[i,j]+1];rect(j-.47,yy-.47,j+.47,yy+.47,col=bg,border=NA);if(!is.na(mat[i,j]))text(j,yy,cell[i,j],cex=.83,col=if(mat[i,j]>10)'white' else '#222222')}
 axis(1,at=seq_len(17),labels=subjects,tick=FALSE,cex.axis=.85);axis(2,at=8:1,labels=labels,tick=FALSE,las=1,cex.axis=.9)
 mtext('Cell = substitutions vs native/core reference; * = terminal-region deficit >=5 aa; grey = not detected',side=1,outer=TRUE,line=.4,cex=.8)
 mtext('Activity estimated from the 15-min figure; sequence features do not establish functional effects',side=1,outer=TRUE,line=1.4,cex=.8)
}
png(file.path(out,'cohort_sequence_overview.png'),width=3300,height=2100,res=250,type='cairo');plot_summary();dev.off()
cairo_pdf(file.path(out,'cohort_sequence_overview.pdf'),width=13.2,height=8.4);plot_summary();dev.off()
cat('Exploratory association tests (shared BH family):\n');print(alltests);cat('Identical sequence groups (n>=2):\n');print(a[a$n>=2,]);cat('Functional site checks:\n');f<-read.delim(file.path(out,'functional_site_checks.tsv'));print(f[tolower(as.character(f$unchanged)) != 'true',])
writeLines(c('All outcomes are approximate figure-digitized 15-min product formation, not raw assay measurements. These post-hoc comparisons are exploratory, not causal tests.','Primary across-subject tests: Spearman associations of activity with RPKM-weighted substitution burden, exact-native fraction, and terminal-deficit abundance fraction. Within-reference substitution models require >=6 subjects and adjust log1p(same-reference RPKM). Common non-gap alternative residues require >=3 carriers and >=3 comparison subjects within the same reference; perfectly co-occurring variants are grouped. Activity is modeled as log1p(gene RPKM) plus carrier status. All test p-values share one BH correction family. No weighting of genes optimized to activity. No activity/abundance ratio used as primary outcome. Native-reference distance is not a functional damage score.','Global rank associations are confounded by reference composition; within-reference models cannot adjust all community, expression, assay or assembly effects. Lack of significant association does not establish absence of effects; n=17 has low power and variant patterns may be linked.','Exact-protein comparisons refer to a single GUS homolog, not identical whole-community enzyme repertoires. Identical protein sequence does not imply identical expression, protein abundance, cellular localization or activity.'),file.path(out,'statistics_methods.txt'))
