a<-commandArgs(FALSE);root<-dirname(normalizePath(gsub('~+~',' ',sub('^--file=','',grep('^--file=',a,value=TRUE)[1]),fixed=TRUE)))
source(file.path(root,'..','shared','application_functions.R'))
base<-file.path(Sys.getenv('KNOWN_OUTPUT',file.path(root,'results')),'point');cfg<-readRDS(file.path(base,'config.rds'))
out<-file.path(Sys.getenv('KNOWN_OUTPUT',file.path(root,'results')),'widened');dir.create(out,recursive=TRUE,showWarnings=FALSE)
B<-as.integer(Sys.getenv('KNOWN_BOOT','500'));R<-cfg$R;n<-cfg$n;Ltrue<-cfg$Lambda;b0<-cfg$b0;taus<-c(.1,.25,.5,.75)
cross<-function(x,lo,hi){hit<-which(lo<=0&hi>=0);if(!length(hit))return(NA_real_);j<-hit[1];if(j==1)return(x[1]);v<-if(lo[j-1]>0)lo else hi;x[j-1]-v[j-1]*(x[j]-x[j-1])/(v[j]-v[j-1])}
one<-function(r){
 dest<-file.path(out,sprintf('rep_%03d.rds',r));if(file.exists(dest))return(NULL)
 N<-n;set.seed(20261009L+r)
 X<-matrix(rbinom(n*12,1,.5),n,12);e<-plogis(b0+.25*(X[,1]-.5)+.25*(X[,2]-.5));lo<-plogis(qlogis(e)-log(Ltrue));hi<-plogis(qlogis(e)+log(Ltrue));q<-(e-lo)/(hi-lo)
 U<-rbinom(n,1,q);rho<-ifelse(U==1,hi,lo);A<-rbinom(n,1,rho);Y<-3200+100*X[,1]-700*U+runif(n,-250,250)
 sys.source(file.path(root,'..','shared','application_functions.R'),envir=environment())
 summarize<-function(w){f1<-elw_side(e,1L,w);f0<-elw_side(e,0L,w);qte<-weighted_quantile(f1$y,f1$weight,taus)-weighted_quantile(f0$y,f0$weight,taus);control<-sum(f0$weight*(f0$y<2500));c(qte,sum(f1$weight*(f1$y<2500))-control,qte-150,sum(f1$weight*(f1$y-150<2500))-control)}
 point<-summarize(rep(1,n));boot<-matrix(NA_real_,B,10)
 for(b in 1:B){set.seed(20261010L+r*1000L+b);boot[b,]<-summarize(rexp(n))}
 margins<-apply(sweep(boot,2,point,'-'),2,quantile,c(.025,.975));margins[1,]<-pmin(0,margins[1,]);margins[2,]<-pmax(0,margins[2,])
 old<-readRDS(file.path(base,sprintf('rep_%03d.rds',r)));rows<-old$rows
 stopifnot(max(abs(rows$baseline-point))<1e-8)
 rows$widened_tip<-NA_real_;rows$widened_truth_included<-NA
 for(k in 1:10){z<-subset(old$curves,scenario==rows$scenario[k]&estimand==rows$estimand[k]);lower<-z$lower+margins[1,k];upper<-z$upper+margins[2,k];j<-which(abs(z$Lambda-Ltrue)<1e-8);rows$widened_tip[k]<-cross(z$Lambda,lower,upper);rows$widened_truth_included[k]<-lower[j]<=rows$truth[k]&rows$truth[k]<=upper[j]}
 # Widening cannot delay first zero inclusion.
 stopifnot(all(ifelse(is.na(rows$widened_tip),Inf,rows$widened_tip)<=ifelse(is.na(rows$tip),Inf,rows$tip)+1e-8))
 saveRDS(list(rows=rows,boot=boot,margins=margins),dest);message('Widening ',r);NULL
}
invisible(parallel::mclapply(1:R,one,mc.cores=as.integer(Sys.getenv('ELW_CORES','2')),mc.preschedule=FALSE,mc.set.seed=FALSE))
z<-lapply(file.path(out,sprintf('rep_%03d.rds',1:R)),readRDS);d<-do.call(rbind,lapply(z,`[[`,'rows'));write.csv(d,file.path(out,'raw.csv'),row.names=FALSE)
s<-do.call(rbind,lapply(split(d,list(d$scenario,d$estimand),drop=TRUE),function(x){p<-ifelse(is.finite(x$tip),x$tip,Inf);w<-ifelse(is.finite(x$widened_tip),x$widened_tip,Inf);data.frame(scenario=x$scenario[1],estimand=x$estimand[1],truth=x$truth[1],point_median=quantile(p,.5,type=1),point_q25=quantile(p,.25,type=1),point_q75=quantile(p,.75,type=1),widened_median=quantile(w,.5,type=1),widened_q25=quantile(w,.25,type=1),widened_q75=quantile(w,.75,type=1),point_over3=sum(!is.finite(p)),widened_over3=sum(!is.finite(w)),point_truth_included=mean(x$truth_in_interval_at_true_Lambda),widened_truth_included=mean(x$widened_truth_included))}))
write.csv(s,file.path(out,'summary.csv'),row.names=FALSE);print(s);capture.output(sessionInfo(),file=file.path(out,'sessionInfo.txt'));saveRDS(list(R=R,B=B,n=n,Lambda=Ltrue,seed=20261010L,code_md5=tools::md5sum(file.path(root,c('run_mnar_widening.R','application_functions.R')))),file.path(out,'config.rds'))
