# Aggregate complete v2 records and attach exact binomial Monte Carlo intervals.
args<-commandArgs(TRUE);stopifnot(length(args)==1L)
f<-sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1]);root<-dirname(normalizePath(gsub('~+~',' ',f,fixed=TRUE)));source(file.path(root,'scalar_results.R'))
files<-list.files(args[1],pattern='^rep_[0-9]+[.]rds$',full.names=TRUE);stopifnot(length(files)>0)
z<-lapply(files,readRDS);stopifnot(all(vapply(z,function(x)isTRUE(x$ok)&&isTRUE(x$marginals_complete),logical(1))))
d<-do.call(rbind,lapply(z,`[[`,'rows'));s<-gbm_summarize_rows(d)
s$coverage_CI_lower<-s$coverage_CI_upper<-NA_real_
for(i in seq_len(nrow(s))){x<-subset(d,scenario==s$scenario[i]&estimand==s$estimand[i]&interval_method==s$interval_method[i]);ci<-binom.test(sum(x$covered),nrow(x))$conf.int;s$coverage_CI_lower[i]<-ci[1];s$coverage_CI_upper[i]<-ci[2]}
write.csv(s,file.path(args[1],'audited-summary.csv'),row.names=FALSE)
