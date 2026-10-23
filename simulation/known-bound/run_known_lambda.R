a<-commandArgs(FALSE);root<-dirname(normalizePath(gsub('~+~',' ',sub('^--file=','',grep('^--file=',a,value=TRUE)[1]),fixed=TRUE)))
source(file.path(root,'sensitivity_functions.R'))
out<-file.path(Sys.getenv('KNOWN_OUTPUT',file.path(root,'results')),'point');dir.create(out,recursive=TRUE,showWarnings=FALSE)
R<-as.integer(Sys.getenv('KNOWN_REPS','100'));n<-as.integer(Sys.getenv('KNOWN_N','2000'));Ltrue<-1.5;grid<-seq(1,3,by=.05);taus<-c(.1,.25,.5,.75)
cells<-expand.grid(x1=0:1,x2=0:1)
b0<-uniroot(function(b)mean(plogis(b+.25*(cells$x1-.5)+.25*(cells$x2-.5)))-.058,c(-10,0))$root
cells$e<-plogis(b0+.25*(cells$x1-.5)+.25*(cells$x2-.5))
cells$lo<-plogis(qlogis(cells$e)-log(Ltrue));cells$hi<-plogis(qlogis(cells$e)+log(Ltrue));cells$q<-(cells$e-cells$lo)/(cells$hi-cells$lo)
stopifnot(max(abs(cells$q*cells$hi+(1-cells$q)*cells$lo-cells$e))<1e-12)
F0<-function(y)mean((1-cells$q)*punif(y,2950+100*cells$x1,3450+100*cells$x1)+cells$q*punif(y,2250+100*cells$x1,2750+100*cells$x1))
truth<-lapply(c(0,-150),function(delta)c(rep(delta,4),F0(2500-delta)-F0(2500)));names(truth)<-c('null','shift')
cross<-function(x,lo,hi){hit<-which(lo<=0&hi>=0);if(!length(hit))return(NA_real_);j<-hit[1];if(j==1)return(x[1]);prev<-j-1;v<-if(lo[prev]>0)lo else hi;x[prev]-v[prev]*(x[j]-x[prev])/(v[j]-v[prev])}
one<-function(r){dest<-file.path(out,sprintf('rep_%03d.rds',r));if(file.exists(dest))return(NULL);set.seed(20261009L+r)
 X<-matrix(rbinom(n*12,1,.5),n,12);e<-plogis(b0+.25*(X[,1]-.5)+.25*(X[,2]-.5));lo<-plogis(qlogis(e)-log(Ltrue));hi<-plogis(qlogis(e)+log(Ltrue));q<-(e-lo)/(hi-lo)
 U<-rbinom(n,1,q);rho<-ifelse(U==1,hi,lo);A<-rbinom(n,1,rho);Y0<-3200+100*X[,1]-700*U+runif(n,-250,250)
 stopifnot(max(abs(pmax((rho/(1-rho))/(e/(1-e)),(e/(1-e))/(rho/(1-rho)))-Ltrue))<1e-10)
 rows<-list();curves<-list()
 for(sc in c('null','shift')){
 delta<-if(sc=='null')0 else -150;Y<-Y0+delta*A
 z<-elw_msm_sensitivity(A,Y,e,grid,taus,threshold=2500)
 for(k in 1:5){v<-z[z$tau==taus[if(k<5)k else 1],];lower<-if(k<5)v$effect_lower else v$lbw_rd_lower;upper<-if(k<5)v$effect_upper else v$lbw_rd_upper;t<-truth[[sc]][k];j<-which(abs(v$Lambda-Ltrue)<1e-9);nm<-c('QTE10','QTE25','QTE50','QTE75','DTE')[k]
 rows[[length(rows)+1]]<-data.frame(rep=r,scenario=sc,estimand=nm,truth=t,baseline=lower[1],tip=cross(v$Lambda,lower,upper),truth_in_interval_at_true_Lambda=lower[j]<=t&t<=upper[j],treated=sum(A))
 curves[[length(curves)+1]]<-data.frame(rep=r,scenario=sc,estimand=nm,Lambda=v$Lambda,lower,upper)
 }
 }
 saveRDS(list(rows=do.call(rbind,rows),curves=do.call(rbind,curves)),dest);message('MNAR ',r);NULL
}
invisible(parallel::mclapply(1:R,one,mc.cores=as.integer(Sys.getenv('ELW_CORES','2')),mc.preschedule=FALSE,mc.set.seed=FALSE))
z<-lapply(file.path(out,sprintf('rep_%03d.rds',1:R)),readRDS);d<-do.call(rbind,lapply(z,`[[`,'rows'));write.csv(d,file.path(out,'raw.csv'),row.names=FALSE)
s<-do.call(rbind,lapply(split(d,list(d$scenario,d$estimand),drop=TRUE),function(x){f<-ifelse(is.finite(x$tip),x$tip,Inf);data.frame(scenario=x$scenario[1],estimand=x$estimand[1],truth=x$truth[1],mean_baseline=mean(x$baseline),median_tip=quantile(f,.5,type=1),q25=quantile(f,.25,type=1),q75=quantile(f,.75,type=1),no_cross_by_3=sum(!is.finite(x$tip)),tip_le_true_Lambda=sum(x$tip<=Ltrue,na.rm=TRUE)/R,truth_in_interval_at_true_Lambda=mean(x$truth_in_interval_at_true_Lambda))}))
write.csv(s,file.path(out,'summary.csv'),row.names=FALSE);saveRDS(list(n=n,R=R,Lambda=Ltrue,grid=grid,seed=20261009L,b0=b0,cells=cells,truth=truth,code_md5=tools::md5sum(file.path(root,c('run_known_lambda.R','sensitivity_functions.R')))),file.path(out,'config.rds'));print(s);capture.output(sessionInfo(),file=file.path(out,'sessionInfo.txt'))
