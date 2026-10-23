# Pure storage/evaluation helpers. These functions neither fit models nor draw RNG.
gbm_effect_targets <- function() c('DTE',paste0('QTE_',c(.1,.25,.5,.75)),paste0('PTE_',1:4))
gbm_side_targets <- function() c('F0_2500minus','F1_2500minus',paste0('Q0_',c(.1,.25,.5,.75)),paste0('Q1_',c(.1,.25,.5,.75)),paste0('P0_',1:4),paste0('P1_',1:4))
gbm_scalar_targets <- function() c(gbm_effect_targets(),gbm_side_targets())
gbm_check_scalar_identities <- function(v,tol=1e-10) {
 stopifnot(identical(names(v),gbm_scalar_targets()),all(is.finite(v)))
 stopifnot(abs(v['DTE']-(v['F1_2500minus']-v['F0_2500minus']))<tol)
 for(t in c(.1,.25,.5,.75))stopifnot(abs(v[paste0('QTE_',t)]-(v[paste0('Q1_',t)]-v[paste0('Q0_',t)]))<tol)
 for(k in 1:4)stopifnot(abs(v[paste0('PTE_',k)]-(v[paste0('P1_',k)]-v[paste0('P0_',k)]))<tol)
 for(a in 0:1){
  pp<-v[paste0('P',a,'_',1:4)];stopifnot(all(pp>=-tol&pp<=1+tol),abs(sum(pp)-1)<tol,abs(v[paste0('F',a,'_2500minus')]-sum(pp[1:2]))<tol)
  stopifnot(all(diff(v[paste0('Q',a,'_',c(.1,.25,.5,.75))])>=0))
 }
 invisible(TRUE)
}
gbm_scalar_from_sides <- function(y0,w0,y1,w1,quantile_fn) {
 tau<-c(.1,.25,.5,.75);bounds<-c(-Inf,1500,2500,4000,Inf)
 stopifnot(length(y0)==length(w0),length(y1)==length(w1),abs(sum(w0)-1)<1e-9,abs(sum(w1)-1)<1e-9)
 f0<-sum(w0[y0<2500]);f1<-sum(w1[y1<2500])
 q0<-quantile_fn(y0,w0,tau);q1<-quantile_fn(y1,w1,tau)
 p0<-vapply(1:4,function(k)sum(w0[y0>=bounds[k]&y0<bounds[k+1]]),numeric(1))
 p1<-vapply(1:4,function(k)sum(w1[y1>=bounds[k]&y1<bounds[k+1]]),numeric(1))
 v<-setNames(c(f1-f0,q1-q0,p1-p0,f0,f1,q0,q1,p0,p1),gbm_scalar_targets())
 gbm_check_scalar_identities(v);v
}
gbm_scalar_truth <- function(cal) {
 stopifnot(identical(cal$design,'empirical_residual_v1'))
 lapply(cal$truth,function(tr){
  lookup<-function(F,y){i<-findInterval(y,tr$grid);c(0,F)[i+1L]}
  tau<-c(.1,.25,.5,.75)
  q0<-vapply(tau,function(t)tr$grid[which(tr$F0>=t)[1]],numeric(1));q1<-vapply(tau,function(t)tr$grid[which(tr$F1>=t)[1]],numeric(1))
  f0<-lookup(tr$F0,2499);f1<-lookup(tr$F1,2499)
  p0<-diff(c(0,lookup(tr$F0,c(1499,2499,3999)),1));p1<-diff(c(0,lookup(tr$F1,c(1499,2499,3999)),1))
  v<-setNames(c(f1-f0,q1-q0,p1-p0,f0,f1,q0,q1,p0,p1),gbm_scalar_targets())
  gbm_check_scalar_identities(v);stopifnot(max(abs(v[2:5]-tr$qte))<1e-10,max(abs(v[6:9]-tr$pte))<1e-10);v
 })
}
gbm_interval_rows <- function(point,boot,truth,n,rep,scenario,method='percentile',level=.95) {
 # Raw draws are authoritative. This choice is explicitly recorded and can be changed later.
 stopifnot(identical(names(point),colnames(boot)),identical(names(truth),names(point)),all(is.finite(boot)),nrow(boot)>=2L)
 alpha<-1-level;se<-apply(boot,2,sd);q<-apply(boot,2,quantile,c(alpha/2,1-alpha/2),names=FALSE)
 if(method=='percentile'){lo<-q[1,];hi<-q[2,]}else if(method=='normal'){lo<-point-qnorm(1-alpha/2)*se;hi<-point+qnorm(1-alpha/2)*se}else if(method=='basic'){lo<-2*point-q[2,];hi<-2*point-q[1,]}else stop('Unsupported interval method: ',method)
 data.frame(n=n,rep=rep,B=nrow(boot),scenario=scenario,estimand=names(point),truth=as.numeric(truth),estimate=as.numeric(point),se=se,lower=lo,upper=hi,covered=lo<=truth&truth<=hi,interval_method=method,level=level,row.names=NULL)
}
gbm_summarize_rows <- function(d) {
 do.call(rbind,lapply(split(d,interaction(d$scenario,d$estimand,d$interval_method,drop=TRUE)),function(x){
  ok<-is.finite(x$estimate)&is.finite(x$se)&!is.na(x$covered);a<-x[ok,];err<-a$estimate-a$truth;r<-nrow(a)
  if(!r)return(NULL)
  data.frame(scenario=a$scenario[1],estimand=a$estimand[1],interval_method=a$interval_method[1],n=a$n[1],R=r,B=a$B[1],bias=mean(err),SD=sqrt(mean((err-mean(err))^2)),mean_boot_SE=mean(a$se),coverage=mean(a$covered),MCSE=sqrt(mean(a$covered)*(1-mean(a$covered))/r),mean_width=mean(a$upper-a$lower))
 }))
}
