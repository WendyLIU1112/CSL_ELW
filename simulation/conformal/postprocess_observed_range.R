# Alternative explicitly bounded search domain; preserve original results.
# Bounds use OBSERVED outcomes only, never Y0/Y1 potential-outcome extrema.
out<-Sys.getenv('EXACT_SIM_OUTPUT');stopifnot(nzchar(out))
dest<-file.path(out,'observed-range');dir.create(dest,showWarnings=FALSE)
files<-list.files(file.path(out,'jobs'),pattern='[.]rds$',full.names=TRUE)
parts<-parallel::mclapply(split(files,rep(1:4,length.out=length(files))),function(paths)do.call(rbind,lapply(paths,function(p){
 z<-readRDS(p);stopifnot(z$success);v<-z$individual;d<-z$job;n<-length(v$Y)
 lo<-min(v$Y);hi<-max(v$Y);end<-v$endpoints
 lower<-pmax(end[,,,1],lo);dim(lower)<-c(n,2L,2L)
 upper<-pmin(end[,,,2],hi);dim(upper)<-c(n,2L,2L)
 # Original sets [L,U); intersect with CLOSED [min observed Y,max observed Y].
 upper_closed<-end[,,,2]>hi
 empty<-v$empty|lower>upper|(lower==upper&!upper_closed)
 cov<-array(FALSE,c(n,2L,2L))
 for(t in 1:2)for(g in 1:2){y<-if(t==1)v$Y0 else v$Y1;cov[,t,g]<-!empty[,t,g]&y>=lower[,t,g]&(y<upper[,t,g]|(upper_closed[,t,g]&y==upper[,t,g]))}
 mis_lower<-ifelse(v$A==1,lower[,1,2],lower[,2,2]);mis_upper<-ifelse(v$A==1,upper[,1,2],upper[,2,2])
 mis_empty<-ifelse(v$A==1,empty[,1,2],empty[,2,2]);mis_uc<-ifelse(v$A==1,upper_closed[,1,2],upper_closed[,2,2])
 il<-ifelse(v$A==1,v$Y-mis_upper,mis_lower-v$Y);iu<-ifelse(v$A==1,v$Y-mis_lower,mis_upper-v$Y)
 lc<-ifelse(v$A==1,mis_uc,TRUE);uc<-ifelse(v$A==1,TRUE,mis_uc)
 truth<-v$Y1-v$Y0
 covered<-!mis_empty&(truth>il|(lc&truth==il))&(truth<iu|(uc&truth==iu))
 missing<-ifelse(v$A==1,cov[,1,2],cov[,2,2]);stopifnot(identical(covered,missing),all(is.finite(il)),all(is.finite(iu)),all(!covered|v$ite_covered))
 # Empty sets have width zero, with empty counts separately reported.
 w<-ifelse(mis_empty,0,iu-il)
 data.frame(panel=d$panel,scenario=d$scenario,n=d$n,iter=d$iter,range_min=lo,range_max=hi,coverage_y0=mean(cov[,1,1]),coverage_y1=mean(cov[,2,1]),coverage_ite_one_missing=mean(covered),coverage_ite_treated=mean(covered[v$A==1]),coverage_ite_control=mean(covered[v$A==0]),width_ite_one_missing=median(w),empty_ite=sum(mis_empty),original_coverage=z$result$coverage_ite_one_missing,lost_coverage=sum(v$ite_covered&!covered))
})),mc.cores=4,mc.set.seed=FALSE)
stopifnot(!any(vapply(parts,inherits,logical(1),'try-error')))
a<-do.call(rbind,parts);stopifnot(nrow(a)==8L*as.integer(Sys.getenv('EXACT_SIM_REPS','1000'))*length(as.integer(strsplit(Sys.getenv('EXACT_SIM_N','5000,10000'),',')[[1]]))/2L);write.csv(a,file.path(dest,'raw.csv'),row.names=FALSE)
metrics<-c('coverage_y0','coverage_y1','coverage_ite_one_missing','coverage_ite_treated','coverage_ite_control','width_ite_one_missing','empty_ite','original_coverage','lost_coverage')
s<-do.call(rbind,lapply(split(a,interaction(a$panel,a$scenario,a$n)),function(d)data.frame(panel=d$panel[1],scenario=d$scenario[1],n=d$n[1],reps=nrow(d),metric=metrics,mean=vapply(d[metrics],mean,numeric(1)),mcse=vapply(d[metrics],function(x)sd(x)/sqrt(nrow(d)),numeric(1)))))
write.csv(s,file.path(dest,'summary.csv'),row.names=FALSE)
writeLines(c('Alternative method: prediction set intersected with [min observed Y,max observed Y] within each replication.','Bounds never use unobserved potential outcomes. This is a data-dependent search restriction, not a known outcome support.','Original calibrated sets and results remain unchanged. Endpoint closure is retained in coverage calculations.','All requested records passed finite-endpoint, subset coverage, and missing-potential-outcome / ITE identity checks.','An empty set has width zero; empty counts are separately saved.','Reproduce exact individual bounds/sets from parent jobs plus range_min/range_max and this frozen script.'),file.path(dest,'README.txt'))
writeLines('All requested saved replications postprocessed and verified.',file.path(dest,'COMPLETE.txt'))
print(s[s$metric%in%c('coverage_ite_one_missing','width_ite_one_missing'),],row.names=FALSE)
