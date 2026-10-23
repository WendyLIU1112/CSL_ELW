# Reconstruct the fixed empirical-residual calibration from authorized cleaned data.
args<-commandArgs(TRUE);stopifnot(length(args)==2L)
script<-normalizePath(gsub('~+~',' ',sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1]),fixed=TRUE));root<-dirname(script)
source(file.path(root,'..','shared','application_functions.R'));dat<-read.csv(args[1],check.names=FALSE)
required<-c('Momage','Height','PreBMI','Momrace','Marital','Parity','Hxpreterm','Diabetes','Hxasthma','HIV','Alcohol','Sex')
stopifnot(all(c(required,'treatment_smoking','outcome_birthweight')%in%names(dat)))
covariates<-setdiff(names(dat),c('row_id','Smoke','Birthweight','treatment_smoking','outcome_birthweight'));stopifnot(setequal(covariates,required));x<-dat[covariates];x[setdiff(required,c('Momage','Height','PreBMI'))]<-lapply(x[setdiff(required,c('Momage','Height','PreBMI'))],factor)
Xraw<-model.matrix(~.-1,data=x);ss<-apply(Xraw,2,sd);X_main<-scale(Xraw[,is.finite(ss)&ss>0,drop=FALSE]);N<-nrow(dat);A<-as.integer(dat$treatment_smoking);Y<-dat$outcome_birthweight
stopifnot(all(A%in%0:1),all(is.finite(Y)),all(is.finite(X_main)))
efit<-fit_propensity(fit_seed=20270621L);offset<-uniroot(function(a)mean(plogis(a+qlogis(efit)))-mean(A),c(-20,20),tol=1e-12)$root;e<-plogis(offset+qlogis(efit))
beta<-lm.fit(cbind(1,X_main[A==0,,drop=FALSE]),Y[A==0])$coefficients;beta[is.na(beta)]<-0;mu<-as.vector(cbind(1,X_main)%*%beta)
delta<--120-50*rank(-mu,ties.method='average')/N
res<-Y[A==0]-mu[A==0]
cal<-list(Xraw=Xraw,mu=mu,sigma=sd(res),delta=delta,e=e,prevalence=mean(A),N=N,covariates=covariates,offset=offset,residual=round(res-mean(res)),mu0=round(mu),mu1=round(mu+delta),design='empirical_residual_v1')
pmf<-function(z){lo<-min(z);hi<-max(z);list(lo=lo,p=tabulate(z-lo+1,nbins=hi-lo+1)/length(z))}
pE<-pmf(cal$residual)
convolve_cdf<-function(m){pM<-pmf(m);size<-length(pE$p)+length(pM$p)-1L;L<-2^ceiling(log2(size));mass<-Re(fft(fft(c(pE$p,rep(0,L-length(pE$p))))*fft(c(pM$p,rep(0,L-length(pM$p)))),inverse=TRUE))/L;mass<-mass[seq_len(size)];stopifnot(min(mass)>-1e-12,abs(sum(mass)-1)<1e-10);mass<-pmax(mass,0);mass<-mass/sum(mass);list(y=pE$lo+pM$lo+seq_len(size)-1L,F=pmin(1,cumsum(mass)))}
f0<-convolve_cdf(cal$mu0);lookup<-function(f,y)c(0,f$F)[findInterval(y,f$y)+1L]
cal$truth<-lapply(list(equality=cal$mu0,shift=cal$mu1),function(m1){f1<-convolve_cdf(m1);grid<-seq(min(f0$y,f1$y)-1,max(f0$y,f1$y)+1);q<-function(f)vapply(c(.1,.25,.5,.75),function(t)f$y[which(f$F>=t)[1]],numeric(1));p<-function(f)diff(c(0,lookup(f,c(1499,2499,3999)),1));ee<-ecdf(cal$residual);set.seed(20260911L);checks<-unique(c(1499,2499,3999,sample(grid,100),q(f0),q(f0)-1,q(f1),q(f1)-1));err<-max(abs(lookup(f0,checks)-vapply(checks,function(y)mean(ee(y-cal$mu0)),numeric(1))),abs(lookup(f1,checks)-vapply(checks,function(y)mean(ee(y-m1)),numeric(1))));stopifnot(err<1e-10);list(grid=grid,F0=lookup(f0,grid),F1=lookup(f1,grid),qte=q(f1)-q(f0),pte=p(f1)-p(f0),numerical_error=err)})
dir.create(args[2],recursive=TRUE,showWarnings=FALSE);saveRDS(cal,file.path(args[2],'calibration.rds'));writeLines(capture.output(sessionInfo()),file.path(args[2],'sessionInfo.txt'));write.csv(data.frame(N=N,prevalence=mean(A),offset=offset),file.path(args[2],'calibration_summary.csv'),row.names=FALSE)
