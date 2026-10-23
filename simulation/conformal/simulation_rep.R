# Kernel replacement of the existing CSL conformal simulation.
  one_rep <- function(iter, n, nonlinear, scenario) {
    sim <- generate_csl_simulation_data(base, n = n, b = 1, heteroskedastic = TRUE,
                                        seed = 20260802 + 10000L * nonlinear + iter)
    if (is.finite(target_prevalence)) {
      Xm <- attr(sim, "X.main")
      eta <- sim$ell
      if (nonlinear) {
        age <- if ("Momage" %in% colnames(Xm)) Xm[, "Momage"] else 0
        bmi <- if ("PreBMI" %in% colnames(Xm)) Xm[, "PreBMI"] else 0
        height <- if ("Height" %in% colnames(Xm)) Xm[, "Height"] else 0
        eta <- eta + 0.60 * (age^2 - 1) + 0.60 * sin(bmi) + 0.50 * age * bmi -
          0.40 * as.numeric(height > 0.75)
      }
      bounded_probability <- function(intercept)
        pmin(ps_upper, pmax(ps_lower, stats::plogis(intercept + eta)))
      intercept <- stats::uniroot(
        function(a) mean(bounded_probability(a)) - target_prevalence,
        c(-30, 30), tol = 1e-12)$root
      sim$pi.true <- bounded_probability(intercept)
      sim$T <- stats::rbinom(n, 1L, sim$pi.true)
      sim$Y <- ifelse(sim$T == 1L, sim$Y1, sim$Y0)
    } else if (nonlinear) {
      Xm <- attr(sim, "X.main")
      age <- if ("Momage" %in% colnames(Xm)) Xm[, "Momage"] else 0
      bmi <- if ("PreBMI" %in% colnames(Xm)) Xm[, "PreBMI"] else 0
      height <- if ("Height" %in% colnames(Xm)) Xm[, "Height"] else 0
      g.nl <- 0.60 * (age^2 - 1) + 0.60 * sin(bmi) + 0.50 * age * bmi -
        0.40 * as.numeric(height > 0.75)
      a.nl <- calibrate_logit_intercept(sim$ell + g.nl, base$smoking.prevalence)
      sim$pi.true <- stats::plogis(a.nl + sim$ell + g.nl)
      sim$T <- stats::rbinom(n, 1L, sim$pi.true)
      sim$Y <- ifelse(sim$T == 1L, sim$Y1, sim$Y0)
    }

    X<-as.data.frame(attr(sim,'X.main'),check.names=TRUE);A<-sim$T;Y<-sim$Y
    set.seed(91000+iter);fold<-integer(n)
    for(t in 0:1)fold[A==t]<-sample(rep(seq_len(K),length.out=sum(A==t)))
    p_oof<-score<-rep(NA_real_,n);models<-vector('list',K)
    for(k in seq_len(K)){
      tr<-which(fold!=k);va<-which(fold==k)
      if(!nonlinear){
        fit<-suppressWarnings(glm(T~.,data=data.frame(T=A[tr],X[tr,,drop=FALSE]),family=binomial()))
        pt<-as.numeric(predict(fit,newdata=X[tr,,drop=FALSE],type='response'))
        p_oof[va]<-as.numeric(predict(fit,newdata=X[va,,drop=FALSE],type='response'))
      }else{
        fit<-gbm::gbm(T~.,data=data.frame(T=A[tr],X[tr,,drop=FALSE]),distribution='bernoulli',n.trees=500L,interaction.depth=2L,shrinkage=.02,n.minobsinnode=100L,bag.fraction=.5,train.fraction=1,verbose=FALSE)
        pt<-as.numeric(predict(fit,newdata=X[tr,,drop=FALSE],n.trees=500L,type='response'))
        p_oof[va]<-as.numeric(predict(fit,newdata=X[va,,drop=FALSE],n.trees=500L,type='response'))
      }
      pca<-weighted_pca_transform(as.matrix(X[tr,,drop=FALSE]),as.matrix(X[va,,drop=FALSE]),n_components=3L)
      arms<-vector('list',2L)
      for(t in 0:1){
        ii<-which(A[tr]==t);w<-elw_group(A[tr],pt,if(t==1)'treated'else'control')$weight
        arms[[t+1L]]<-exact_kernel_model(pca$train[ii,,drop=FALSE],Y[tr][ii],w,normal_reference_bandwidth(length(ii),3L))
        for(j in which(A[va]==t))score[va[j]]<-exact_cdf_score(arms[[t+1L]],pca$test[j,],Y[va[j]])
      }
      models[[k]]<-list(test=va,z=pca$test,arms=arms)
    }
    stopifnot(all(is.finite(p_oof)),all(is.finite(score)))
    cal<-vector('list',2L)
    for(t in 0:1){
      ii<-which(A==t);w<-elw_group(A,p_oof,if(t==1)'treated'else'control')$weight
      fac<-if(t==0)p_oof[ii]else 1-p_oof[ii]
      cal[[t+1L]]<-list(population=exact_calibration(score[ii],w),opposite=exact_calibration(score[ii],w*fac))
    }
    # Every endpoint and empty-set flag is retained; sets are [L,U).
    endpoints<-array(NA_real_,c(n,2L,2L,2L),dimnames=list(NULL,c('Y0','Y1'),c('population','opposite'),c('lower','upper')))
    empty<-array(FALSE,c(n,2L,2L));covered<-array(FALSE,c(n,2L,2L))
    for(k in seq_len(K)){
      m<-models[[k]]
      for(j in seq_along(m$test)){
        i<-m$test[j]
        for(t in 0:1){
          w<-exact_local_weights(m$arms[[t+1L]],m$z[j,])
          truth<-if(t==0)sim$Y0[i]else sim$Y1[i]
          for(g in 1:2){
            s<-exact_prediction_set(m$arms[[t+1L]],m$z[j,],cal[[t+1L]][[g]],weights=w)
            endpoints[i,t+1L,g,]<-c(s$lower,s$upper);empty[i,t+1L,g]<-s$empty
            covered[i,t+1L,g]<-exact_set_contains(s,truth)
          }
        }
      }
    }
    ite_lower<-ifelse(A==1,Y-endpoints[,1,2,2],endpoints[,2,2,1]-Y)
    ite_upper<-ifelse(A==1,Y-endpoints[,1,2,1],endpoints[,2,2,2]-Y)
    ite_empty<-ifelse(A==1,empty[,1,2],empty[,2,2])
    ite_true<-sim$Y1-sim$Y0
    ite_cover<-!ite_empty&ifelse(A==1,ite_true>ite_lower&ite_true<=ite_upper,ite_true>=ite_lower&ite_true<ite_upper)
    # Independent identity: ITE coverage equals coverage of the missing potential outcome.
    missing_cover<-ifelse(A==1,covered[,1,2],covered[,2,2]);stopifnot(identical(ite_cover,missing_cover))
    raw<-data.frame(scenario=scenario,n=n,iter=iter,coverage_y0=mean(covered[,1,1]),coverage_y1=mean(covered[,2,1]),coverage_y0_target1=mean(covered[A==1,1,2]),coverage_y1_target0=mean(covered[A==0,2,2]),coverage_ite_treated=mean(ite_cover[A==1]),coverage_ite_control=mean(ite_cover[A==0]),coverage_ite_one_missing=mean(ite_cover),width_y0=median(endpoints[,1,1,2]-endpoints[,1,1,1]),width_y1=median(endpoints[,2,1,2]-endpoints[,2,1,1]),width_ite_one_missing=median(ite_upper-ite_lower),treated_n=sum(A==1),empty_ite=sum(ite_empty),empty_y0=sum(empty[,1,1]),empty_y1=sum(empty[,2,1]))
    list(result=raw,individual=list(id=sim$id,A=A,Y=Y,Y0=sim$Y0,Y1=sim$Y1,pi_true=sim$pi.true,fold=fold,pi_oof=p_oof,score=score,endpoints=endpoints,empty=empty,ite_lower=ite_lower,ite_upper=ite_upper,ite_empty=ite_empty,ite_covered=ite_cover),calibration=cal)
  }
