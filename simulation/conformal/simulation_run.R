args<-commandArgs(FALSE);root<-dirname(normalizePath(gsub('~+~',' ',sub('^--file=','',grep('^--file=',args,value=TRUE)[1]),fixed=TRUE)))
source(file.path(root,'engine.R'));source(file.path(root,'simulation_rep.R'))
legacy<-file.path(root,'support')
source(file.path(legacy,'csl_smoking_all_functions.R'));source(file.path(legacy,'pca_kernel_utils.R'))
output<-Sys.getenv('EXACT_SIM_OUTPUT');stopifnot(nzchar(output));dir.create(output,recursive=TRUE,showWarnings=FALSE);dir.create(file.path(output,'jobs'),showWarnings=FALSE)
R<-as.integer(Sys.getenv('EXACT_SIM_REPS','1000'));cores<-as.integer(Sys.getenv('EXACT_SIM_CORES','12'));K<-5L
spec<-list(schema='exact_conformal_sim_v1',input_md5=unname(tools::md5sum(Sys.getenv('CSL_SIM_DATA'))),R=R,K=K,n=as.integer(strsplit(Sys.getenv('EXACT_SIM_N','5000,10000'),',')[[1]]),panels=c('regular','limited'),ps_CE_L='PLE',ps_CE_NL=list(trees=500L,depth=2L,node=100L,shrinkage=.02,bag=.5),kernel='PCA3 Gaussian normal reference, direct CDF, training-fold transforms',seed_data='20260802+10000*nonlinear+iter',seed_fold_gbm='91000+iter',code_md5=tools::md5sum(file.path(root,c('engine.R','simulation_rep.R','simulation_run.R','postprocess_observed_range.R'))),legacy_md5=tools::md5sum(file.path(legacy,c('csl_smoking_all_functions.R','pca_kernel_utils.R'))))
path<-file.path(output,'settings.rds');if(file.exists(path))stopifnot(identical(spec,readRDS(path)))else saveRDS(spec,path)
writeLines(capture.output(sessionInfo()),file.path(output,'sessionInfo.txt'))
grid<-expand.grid(n=as.integer(strsplit(Sys.getenv('EXACT_SIM_N','5000,10000'),',')[[1]]),nonlinear=c(FALSE,TRUE),panel=c('regular','limited'),iter=seq_len(R),stringsAsFactors=FALSE)
grid$scenario<-ifelse(grid$nonlinear,'CE-NL','CE-L');grid$file<-sprintf('%s_%s_n%d_r%04d.rds',grid$panel,grid$scenario,grid$n,grid$iter)
write.csv(grid,file.path(output,'job_manifest.csv'),row.names=FALSE)
base<-make_csl_simulation_base(read.csv(Sys.getenv('CSL_SIM_DATA'),check.names=FALSE))
run_job<-function(j){
 d<-grid[j,];p<-file.path(output,'jobs',d$file);if(file.exists(p))return(NULL)
 target_prevalence<<-if(d$panel=='regular').30 else NA_real_;ps_lower<<-if(d$panel=='regular').05 else 0;ps_upper<<-if(d$panel=='regular').95 else 1
 start<-proc.time()[3];warnings<-character()
 ans<-tryCatch(withCallingHandlers({
  z<-one_rep(d$iter,d$n,d$nonlinear,d$scenario)
  stopifnot(!anyNA(z$result),all(is.finite(z$result$coverage_ite_one_missing)))
  c(list(success=TRUE,error=NULL),z)
 },warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')}),error=function(e)list(success=FALSE,error=conditionMessage(e)))
 ans$seconds<-unname(proc.time()[3]-start);ans$warnings<-unique(warnings);ans$job<-d
 saveRDS(ans,paste0(p,'.tmp'));stopifnot(file.rename(paste0(p,'.tmp'),p))
 cat(format(Sys.time()),d$file,'success=',ans$success,'seconds=',round(ans$seconds,1),'\n')
 NULL
}
missing<-which(!file.exists(file.path(output,'jobs',grid$file)))
cat('Exact Conformal jobs:',nrow(grid),'pending:',length(missing),'workers:',cores,'\n')
invisible(parallel::mclapply(missing,run_job,mc.cores=cores,mc.preschedule=FALSE,mc.set.seed=FALSE))
res<-lapply(file.path(output,'jobs',grid$file),readRDS)
status<-do.call(rbind,lapply(res,function(z)data.frame(z$job,success=z$success,seconds=z$seconds,error=if(is.null(z$error))''else z$error,warnings=paste(z$warnings,collapse=' | '))))
write.csv(status,file.path(output,'job_status.csv'),row.names=FALSE)
raw<-do.call(rbind,lapply(res,function(z)if(z$success)cbind(panel=z$job$panel,z$result)))
write.csv(raw,file.path(output,'raw.csv'),row.names=FALSE)
metrics<-setdiff(names(raw),c('panel','scenario','n','iter'))
summary<-do.call(rbind,lapply(split(raw,interaction(raw$panel,raw$scenario,raw$n)),function(d)data.frame(panel=d$panel[1],scenario=d$scenario[1],n=d$n[1],completed=nrow(d),metric=metrics,mean=vapply(d[metrics],mean,numeric(1)),mcse=vapply(d[metrics],function(x)sd(x)/sqrt(length(x)),numeric(1)))))
write.csv(summary,file.path(output,'summary.csv'),row.names=FALSE)
if(any(!status$success))stop('At least one job failed; inspect job_status.csv')
writeLines('All requested exact-formula Conformal simulation jobs completed successfully.',file.path(output,'COMPLETE.txt'))

# Final reporting version selected by the user: observed-outcome search range.
source(file.path(root,'postprocess_observed_range.R'),local=FALSE)
