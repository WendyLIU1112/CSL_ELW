#!/usr/bin/env Rscript
# GBM interval assessment: scalar effects, marginals, and bootstrap draws.
args<-commandArgs(FALSE);script<-normalizePath(gsub('~+~',' ',sub('^--file=','',grep('^--file=',args,value=TRUE)[1]),fixed=TRUE))
root<-dirname(script);source(file.path(root,'..','shared','application_functions.R'));source(file.path(root,'scalar_results.R'))
reference<-Sys.getenv('CSL_GBM_REFERENCE',file.path(root,'results','calibration'))
calfile<-Sys.getenv('CSL_GBM_CALIBRATION',file.path(reference,'calibration.rds'))
out<-Sys.getenv('CSL_GBM_OUTPUT',file.path(root,'results-with-marginals-v2'))
n<-as.integer(Sys.getenv('CSL_GBM_N','10000'));R<-as.integer(Sys.getenv('CSL_GBM_REPS','100'));B<-as.integer(Sys.getenv('CSL_GBM_BOOT','500'));cores<-as.integer(Sys.getenv('CSL_GBM_CORES','12'))
method<-Sys.getenv('CSL_GBM_INTERVAL','percentile')
stopifnot(n>=500,R>0,B>=2,cores>0,method%in%c('percentile','normal','basic'))
cal<-readRDS(calfile);truth<-gbm_scalar_truth(cal);targets<-gbm_scalar_targets()
config<-list(schema_version=2L,n=n,R=R,B=B,cores=cores,spec=gbm_spec,data_seed=20260910L,point_seed=20270621L,calibration_md5=unname(tools::md5sum(calfile)),function_md5=unname(tools::md5sum(file.path(root,'..','shared','application_functions.R'))),storage_md5=unname(tools::md5sum(file.path(root,'scalar_results.R'))),targets=targets,save_bootstrap_ps=TRUE,threshold='Y < 2500',interval_method=method)
cat('PLAN schema=2 n=',n,' R=',R,' B=',B,' workers=',cores,' output=',out,'\n',sep='')
if(!'--run'%in%args){cat('No simulation started. Explicit --run is required.\n');quit(status=0)}
dir.create(out,recursive=TRUE,showWarnings=FALSE)
cfg<-file.path(out,'config.rds');if(file.exists(cfg))stopifnot(identical(readRDS(cfg),config))else saveRDS(config,cfg)
atomic<-function(x,p){saveRDS(x,paste0(p,'.tmp'),compress=TRUE);stopifnot(file.rename(paste0(p,'.tmp'),p))}
# Completed records can be reused when only the total replication count or
# execution parallelism changes. All scientific settings must match exactly.
same_scientific_config<-function(a,b){
 execution<-c('R','cores')
 identical(a[!names(a)%in%execution],b[!names(b)%in%execution])
}
valid<-function(z,j){
 if(!isTRUE(z$ok)||!identical(z$schema_version,2L)||!same_scientific_config(z$config,config)||!identical(z$rep,j)||!isTRUE(z$marginals_complete))return(FALSE)
 stopifnot(j<=z$config$R,j<=R)
 stopifnot(identical(dim(z$bootstrap_ps),c(n,B)),all(is.finite(z$bootstrap_ps)),length(z$point_ps)==n)
 for(s in names(truth)){
  stopifnot(identical(dim(z$boot[[s]]),c(B,length(targets))),identical(colnames(z$boot[[s]]),targets))
  gbm_check_scalar_identities(z$point[[s]])
  for(b in seq_len(B))gbm_check_scalar_identities(z$boot[[s]][b,])
 }
 TRUE
}
lock<-file.path(out,'running.lock');if(!dir.create(lock,showWarnings=FALSE))stop('Output locked; verify process before restarting.')
writeLines(as.character(Sys.getpid()),file.path(lock,'pid'))
make_sample<-function(j){
 e<-new.env(parent=globalenv());e$N<-n;seed<-20260910L+n*100L+j;set.seed(seed)
 id<-sample.int(cal$N,n,replace=TRUE);x<-cal$Xraw[id,,drop=FALSE];ss<-apply(x,2,sd);keep<-is.finite(ss)&ss>0
 e$X_main<-scale(x[,keep,drop=FALSE]);e$A<-rbinom(n,1,cal$e[id]);stopifnot(sum(e$A)>1,sum(e$A==0)>1)
 y0<-cal$mu0[id]+sample(cal$residual,n,replace=TRUE);y1<-y0+cal$mu1[id]-cal$mu0[id]
 e$ys<-list(equality=y0,shift=ifelse(e$A==1,y1,y0));e$Y<-y0
 e$sample<-list(id=id,A=e$A,Y0=y0,Y1=list(equality=y0,shift=y1),observed=e$ys,true_ps=cal$e[id],design_columns=colnames(e$X_main),scale_center=attr(e$X_main,'scaled:center'),scale_scale=attr(e$X_main,'scaled:scale'))
 sys.source(file.path(root,'..','shared','application_functions.R'),envir=e)
 e$evaluate<-eval(quote(function(p,w){f1<-elw_side(p,1L,w);f0<-elw_side(p,0L,w);lapply(ys,function(y)gbm_scalar_from_sides(y[A==0],f0$weight,y[A==1],f1$weight,weighted_quantile))}),envir=e)
 e
}
one<-function(j){
 dest<-file.path(out,sprintf('rep_%03d.rds',j));if(file.exists(dest)){stopifnot(valid(readRDS(dest),j));return(NULL)}
 start<-proc.time()[3];e<-make_sample(j);p<-e$fit_propensity(fit_seed=20270621L);point<-e$evaluate(p,rep(1,n))
 boot<-lapply(names(truth),function(s)matrix(NA_real_,B,length(targets),dimnames=list(NULL,targets)));names(boot)<-names(truth)
 ps_boot<-matrix(NA_real_,n,B);done<-0L;cpfile<-file.path(out,sprintf('checkpoint_%03d.rds',j))
 if(file.exists(cpfile)){cp<-readRDS(cpfile);stopifnot(identical(cp$config,config),isTRUE(all.equal(cp$point,point,tolerance=1e-12)));boot<-cp$boot;ps_boot<-cp$bootstrap_ps;done<-cp$done}
 # Recorded seeds and propensity scores support reconstruction of the ELW weights.
 boot_seeds<-20260910L+n*10000L+j*1000L+seq_len(B)
 cat('REP_START rep=',j,' done=',done,'\n',sep='');flush.console()
 if(done<B)for(b in seq.int(done+1L,B)){
  set.seed(boot_seeds[b]);w<-rexp(n);pb<-e$fit_propensity(w);v<-e$evaluate(pb,w);ps_boot[,b]<-pb
  for(s in names(truth))boot[[s]][b,]<-v[[s]]
  if(b%%50L==0L||b==B){atomic(list(config=config,point=point,boot=boot,bootstrap_ps=ps_boot,done=b),cpfile);cat('BOOT rep=',j,' done=',b,'/',B,'\n',sep='');flush.console()}
 }
 rows<-do.call(rbind,lapply(names(truth),function(s)gbm_interval_rows(point[[s]],boot[[s]],truth[[s]],n,j,s,method)))
 z<-list(schema_version=2L,ok=TRUE,rep=j,rows=rows,point=point,boot=boot,truth=truth,sample=e$sample,point_ps=p,bootstrap_ps=ps_boot,seeds=list(data=20260910L+n*100L+j,point=20270621L,bootstrap=boot_seeds),spec=gbm_spec,config=config,marginals_complete=TRUE,elapsed=unname(proc.time()[3]-start))
 stopifnot(valid(z,j));atomic(z,dest);unlink(cpfile)
 cat('REP_DONE rep=',j,' elapsed=',round(z$elapsed,1),'s\n',sep='');flush.console();NULL
}
selected<-seq_len(R)
res<-parallel::mclapply(selected,function(j)tryCatch(one(j),error=function(e){atomic(list(rep=j,error=conditionMessage(e)),file.path(out,sprintf('error_%03d.rds',j)));message('FAILED rep=',j,': ',conditionMessage(e));NULL}),mc.cores=cores,mc.preschedule=FALSE,mc.set.seed=FALSE)
newfiles<-file.path(out,sprintf('rep_%03d.rds',seq_len(R)));has<-file.exists(newfiles)
rows<-lapply(which(has),function(j){z<-readRDS(newfiles[j]);stopifnot(valid(z,j));z$rows})
if(length(rows)){d<-do.call(rbind,rows);write.csv(d,file.path(out,'scalar-raw.csv'),row.names=FALSE);write.csv(gbm_summarize_rows(d),file.path(out,'scalar-summary.csv'),row.names=FALSE)}
# Each estimand carries its own R. Missing single-side results are not counted as completed.
write.csv(data.frame(rep=seq_len(R),complete=has),file.path(out,'availability.csv'),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
missing_selected<-selected[!has[selected]]
if(length(missing_selected)){
 writeLines(paste('Incomplete selected replications:',paste(missing_selected,collapse=',')),file.path(out,'INCOMPLETE.txt'))
 unlink(lock,recursive=TRUE);stop('Some selected replications failed; retained results/checkpoints and error records.')
}
if(all(has))writeLines(paste('Completed',format(Sys.time(),tz='Asia/Singapore',usetz=TRUE),R,'/',R,'replications; all single-side and effect draws saved.'),file.path(out,'COMPLETE.txt'))
unlink(lock,recursive=TRUE);cat('DONE selected=',length(selected),' schema2_complete=',sum(has),'/',R,'\n',sep='')
