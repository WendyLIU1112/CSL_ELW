# Direct implementation of Fhat(y|x) and {y: a/2 <= G(Fhat(y|x)) <= 1-a/2}.
# No outcome bins, quantile grid, spatial grid, clipping, or nearest-grid inversion.
exact_kernel_model <- function(z, y, mass, bandwidth) {
 z<-as.matrix(z);stopifnot(nrow(z)==length(y),length(y)==length(mass),all(is.finite(z)),all(is.finite(y)),all(is.finite(mass)),all(mass>=0),sum(mass)>0,is.finite(bandwidth),bandwidth>0)
 keep<-mass>0;o<-order(y[keep]);zz<-z[keep,,drop=FALSE][o,,drop=FALSE];yy<-y[keep][o]
 list(z=zz,y=yy,logmass=log(mass[keep][o]),h=bandwidth)
}
exact_local_weights <- function(model, x) {
 stopifnot(length(x)==ncol(model$z),all(is.finite(x)))
 logw<-model$logmass-rowSums(sweep(model$z,2,x,'-')^2)/(2*model$h^2)
 w<-exp(logw-max(logw));w/sum(w)
}
exact_cdf_score <- function(model, x, y) {
 # Held-out outcome is used only as the evaluation threshold, never in fitting.
 w<-exact_local_weights(model,x);j<-findInterval(y,model$y);if(j==0L)return(0);if(j==length(w))return(1);cdf<-cumsum(w);cdf[j]/tail(cdf,1)
}
exact_calibration <- function(scores, mass, alpha=.05) {
 stopifnot(length(scores)==length(mass),all(is.finite(scores)),all(scores>=0&scores<=1),all(is.finite(mass)),all(mass>=0),sum(mass)>0,alpha>0,alpha<1)
 keep<-mass>0;o<-order(scores[keep]);u<-scores[keep][o];w<-mass[keep][o]
 ends<-c(which(diff(u)!=0),length(u));cdf<-cumsum(w)/sum(w);cdf[length(cdf)]<-1
 u<-u[ends];cdf<-cdf[ends]
 # Exactly invert the two inequalities for the RIGHT-continuous empirical G.
 # G(v)>=alpha/2 iff v>=lower_score; G(v)<=1-alpha/2 iff v<upper_score.
 list(lower_score=u[which(cdf>=alpha/2)[1]],upper_score=u[which(cdf>1-alpha/2)[1]],scores=u,cdf=cdf,alpha=alpha)
}
exact_prediction_set <- function(model,x,calibration,weights=NULL) {
 w<-if(is.null(weights))exact_local_weights(model,x)else weights
 ends<-c(which(diff(model$y)!=0),length(model$y));yy<-model$y[ends];cw<-cumsum(w);ff<-(cw/tail(cw,1))[ends];ff[length(ff)]<-1
 l<-calibration$lower_score;u<-calibration$upper_score
 lower<-if(l==0)-Inf else yy[which(ff>=l)[1]]
 upper<-yy[which(ff>=u)[1]]
 empty<-u<=l||is.na(lower)||is.na(upper)||lower>=upper
 # [lower,upper): endpoint inclusion follows the formula, including ties.
 list(lower=lower,upper=upper,lower_closed=TRUE,upper_closed=FALSE,empty=empty)
}
exact_set_contains <- function(set,y) !set$empty&y>=set$lower&y<set$upper
exact_ite_set <- function(missing, observed, treatment) {
 if(treatment==1L)list(lower=observed-missing$upper,upper=observed-missing$lower,lower_closed=FALSE,upper_closed=TRUE,empty=missing$empty)
 else list(lower=missing$lower-observed,upper=missing$upper-observed,lower_closed=TRUE,upper_closed=FALSE,empty=missing$empty)
}
exact_ite_direction <- function(set) {
 if(set$empty)return('empty')
 if(set$upper<0||(set$upper==0&&!set$upper_closed))return('negative')
 if(set$lower>0||(set$lower==0&&!set$lower_closed))return('positive')
 'contains_zero'
}
