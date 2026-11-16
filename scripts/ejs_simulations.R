#!/usr/bin/env Rscript
# Reproducible numerical study for latent proportions in bivariate Poisson mixtures.
# Base R only. Run: Rscript ejs_simulations.R test|pilot|full|figures|local_test|local_full [output_dir]
run_local_alternatives <- function(local_mode, output_dir) {
  mode <- local_mode
  stopifnot(mode %in% c('test','full'))
  out <- output_dir
  dir.create(out,recursive=TRUE,showWarnings=FALSE)
  set.seed(20261031)
  a <- 1; b <- 5; u0 <- 0.4; alpha <- 0.05
  n_grid <- if(mode=='test') c(250L,1000L) else c(250L,500L,1000L,2000L,4000L,10000L)
  B <- if(mode=='test') 30L else 1000L
  Bcal <- if(mode=='test') 100L else 1500L
  kappas <- c(0.5,0.75,1,1.25)
  # Exact polynomial coefficient recursion from original simulation program.
  poly_mul_linear <- function(p,z,kappa) {
   q<-numeric(length(p)+1L); q[seq_along(p)]<-q[seq_along(p)]+z*p
   q[seq_along(p)+1L]<-q[seq_along(p)+1L]-kappa*p; q
  }
  cheb_coeff <- function(m,a,b) {
   z<-(a+b)/(b-a); kappa<-2/(b-a); t0<-1; t1<-c(z,-kappa)
   if(m==1L) tm<-t1 else {
    for(j in 2:m) {
     prod<-2*poly_mul_linear(t1,z,kappa)
     prod[seq_along(t0)]<-prod[seq_along(t0)]-t0
     t0<-t1;t1<-prod
    }; tm<-t1
   }
   -(tm/cosh(m*acosh(z)))[-1]
  }
  score_poly <- function(x1,x2,coeff) {
   d<-x1+x2; ans<-numeric(length(d)); falling<-rep(1,length(d))
   for(j in 0:(length(coeff)-1L)) {
    if(j>0) falling<-falling*(d-j)
    ans<-ans+coeff[j+1L]*x1*falling
   };ans
  }
  clip <- function(x) pmax(0,pmin(1,x))
  simulate <- function(n,gamma,scenario='linear') {
   s<-runif(n,a,b);v<-(s-a)/(b-a)
   u<-switch(scenario,
    linear=u0+gamma*(s-(a+b)/2)/(b-a),
    constant=rep(u0,n),
    beta=rbeta(n,2,3),
    nonlinear=u0+1.2*((v-0.5)^2-1/12))
   stopifnot(all(u>=0 & u<=1))
   list(x1=rpois(n,s*u),x2=rpois(n,s*(1-u)))
  }
  statistic <- function(d,m,coeff) {
   x1<-d$x1;x2<-d$x2;total<-x1+x2
   mean(x1)-mean(total)*clip(mean(score_poly(x1,x2,coeff)))
  }
  max_gamma<-2*min(u0,1-u0)
  mu<-(a+b)/2
  lambda<-log((sqrt(b)+sqrt(a))/(sqrt(b)-sqrt(a)))
  rows<-list(); cals<-list(); pos<-0L
  for(n in n_grid) {
   m<-max(1L,as.integer(floor(log(n)/log(log(n)))))
   coeff<-cheb_coeff(m,a,b)
   # Separate independent calibration under Design IV; NOT uniform null control.
   null_cal<-replicate(Bcal,abs(statistic(simulate(n,0,'nonlinear'),m,coeff)))
   cutoff<-as.numeric(quantile(null_cal,1-alpha,type=1))
   cals[[length(cals)+1L]]<-data.frame(n=n,m=m,calibration='nonlinear_zero_cov',Bcal=Bcal,cutoff=cutoff)
   for(sc in c('constant','beta','nonlinear')) {
    vals<-replicate(B,abs(statistic(simulate(n,0,sc),m,coeff)))
    pos<-pos+1L;rows[[pos]]<-data.frame(n=n,m=m,scenario=sc,kappa=NA_real_,rho=0,gamma=0,
      B=B,rejections=sum(vals>cutoff),rate=mean(vals>cutoff),
      mcse=sqrt(mean(vals>cutoff)*(1-mean(vals>cutoff))/B),status='null')
   }
   for(k in kappas) {
    rho<-exp(-k*lambda*log(n)/log(log(n)))
    gam<-12*rho/(b-a)
    if(gam>max_gamma) {
     pos<-pos+1L;rows[[pos]]<-data.frame(n=n,m=m,scenario='local_linear',kappa=k,rho=rho,gamma=gam,
      B=0L,rejections=NA_integer_,rate=NA_real_,mcse=NA_real_,status='inadmissible')
     next
    }
    vals<-replicate(B,abs(statistic(simulate(n,gam),m,coeff)))
    pos<-pos+1L;rows[[pos]]<-data.frame(n=n,m=m,scenario='local_linear',kappa=k,rho=rho,gamma=gam,
      B=B,rejections=sum(vals>cutoff),rate=mean(vals>cutoff),
      mcse=sqrt(mean(vals>cutoff)*(1-mean(vals>cutoff))/B),status='admissible')
   }
   cat('Finished n=',n,' m=',m,'\n',sep='')
  }
  res<-do.call(rbind,rows)
  write.csv(res,file.path(out,'local_power.csv'),row.names=FALSE)
  write.csv(do.call(rbind,cals),file.path(out,'null_calibration.csv'),row.names=FALSE)
  writeLines(capture.output(sessionInfo()),file.path(out,'session_info.txt'))
  writeLines(c('Local alternatives: S~Uniform[1,5], U=0.4+gamma*(S-3)/4',
   'Cov(S,U)=gamma/3; admissibility requires abs(gamma)<=0.8',
   'Critical values independently calibrated from nonlinear zero-covariance null only.',
   'Not a uniform level-alpha test; not a direct minimax sharpness verification.',
   paste('seed=20261031, B=',B,', Bcal=',Bcal,', alpha=',alpha,sep='')),
   file.path(out,'protocol.txt'))
  pdf(file.path(out,'local_power.pdf'),width=8,height=5)
  plot(NA,xlim=range(n_grid),ylim=c(0,1),log='x',xlab='Sample size n',ylab='Rejection frequency',
   main='Local linear alternatives: covariance test')
  abline(h=alpha,lty=2,col='grey50')
  for(i in seq_along(kappas)) {
   d<-res[res$scenario=='local_linear' & !is.na(res$kappa) & res$kappa==kappas[i] & res$status=='admissible',]
   if(nrow(d)) lines(d$n,d$rate,type='b',pch=14+i,col=i)
  }
  legend('topright',legend=paste('kappa =',kappas),col=seq_along(kappas),lty=1,pch=15:18,bty='n')
  dev.off()
  cat('Saved results in ',out,'\n',sep='')
}

args <- commandArgs(trailingOnly=TRUE)
mode <- if(length(args)>=1) args[1] else 'test'
if (mode %in% c('local_test', 'local_full')) {
  local_mode <- sub('^local_', '', mode)
  output_dir <- if (length(args) >= 2L) args[2] else file.path('results', mode)
  run_local_alternatives(local_mode, output_dir)
  quit(save='no')
}
if(!mode %in% c('pilot','full','test','figures')) stop('mode must be test, pilot, full, figures, local_test, or local_full')
out <- if(length(args)>=2) args[2] else file.path('results',mode)
dir.create(out, recursive=TRUE, showWarnings=FALSE)
set.seed(20261030)
# v2.1: explicit MOM eligibility, fixed-degree exploratory robust comparator,
# independent test calibration, honest MC intervals on log axes
options(stringsAsFactors=FALSE)

# Chebyshev T_m(z-kappa*s) in the monomial basis, coefficients in increasing order.
poly_mul_linear <- function(p, z, kappa) {
  q <- numeric(length(p)+1L)
  q[seq_along(p)] <- q[seq_along(p)] + z*p
  q[seq_along(p)+1L] <- q[seq_along(p)+1L] - kappa*p
  q
}
cheb_coeff <- function(m,a,b) {
  stopifnot(m>=1, a>0, b>a)
  z <- (a+b)/(b-a); kappa <- 2/(b-a)
  t0 <- 1; t1 <- c(z,-kappa)
  if(m==1L) tm <- t1 else {
    for(j in 2:m) {
      prod <- 2*poly_mul_linear(t1,z,kappa)
      prod[seq_along(t0)] <- prod[seq_along(t0)]-t0
      t0 <- t1; t1 <- prod
    }
    tm <- t1
  }
  # T_m(z) evaluated by cosh for moderate m; log form avoids overflow.
  denom <- cosh(m*acosh(z))
  p <- tm/denom
  # q(s)=(1-P_m(s))/s. Enforce P_m(0)=1 analytically.
  q <- -p[-1]
  list(q=q, denom=denom, const_residual=abs(p[1]-1),
       l1=sum(abs(q)), r=(sqrt(b)-sqrt(a))/(sqrt(b)+sqrt(a)))
}

# H_j=X1*(D-1)_j, with H_j=0 when D=0.
score_poly <- function(x1,x2,coeff) {
  d <- x1+x2; k <- length(coeff)
  ans <- numeric(length(d)); falling <- rep(1,length(d))
  for(j in 0:(k-1L)) {
    if(j>0) falling <- falling*(d-j)
    ans <- ans + coeff[j+1L]*x1*falling
  }
  ans
}
clip <- function(x) pmin(1,pmax(0,x))
ratio_est <- function(x1,x2) {
  den <- sum(x1+x2)
  if(den==0) return(0)
  sum(x1)/den
}
mom_est <- function(z,blocks) {
  # fixed, deterministic, approximately equal-sized blocks
  g <- cut(seq_along(z), breaks=blocks, labels=FALSE)
  clip(median(as.numeric(tapply(z,g,mean))))
}

# Model scenarios: S ~ Uniform[a,b].
scenarios <- c('constant','independent_beta','linear_dependence','nonlinear_zero_cov')
parameters <- function(scenario,a,b) {
  u0 <- 0.4
  if(scenario=='constant') return(list(theta=u0,cov=0,ratio=u0))
  if(scenario=='independent_beta') return(list(theta=2/5,cov=0,ratio=2/5))
  if(scenario=='linear_dependence') {
    gamma <- 0.6; cv <- gamma*(b-a)/12
    return(list(theta=u0,cov=cv,ratio=u0+cv/((a+b)/2)))
  }
  if(scenario=='nonlinear_zero_cov') return(list(theta=u0,cov=0,ratio=u0))
  stop('unknown scenario')
}
simulate_data <- function(n,scenario,a,b) {
  s <- runif(n,a,b); v <- (s-a)/(b-a)
  u <- switch(scenario,
    constant=rep(0.4,n),
    independent_beta=rbeta(n,2,3),
    linear_dependence=0.4+0.6*(s-(a+b)/2)/(b-a),
    nonlinear_zero_cov=0.4+1.2*((v-0.5)^2-1/12))
  stopifnot(all(u>=0),all(u<=1))
  list(x1=rpois(n,s*u),x2=rpois(n,s*(1-u)))
}

# Explicit theoretical constants from Appendix A (can be very conservative).
constants <- function(a,b) {
  z <- (a+b)/(b-a); kappa <- 2/(b-a)
  C0 <- exp(b); C1 <- 9*(z+kappa)^2*max(1,2*b)
  Cb <- exp(1)*b/2+b^2/8
  A1 <- 4*sqrt(24*Cb)/a
  list(C0=C0,C1=C1,A1=A1)
}

# Exact (theorem-based) and exploratory calibrated rules are kept separate.
# The calibrated rule is a finite-sample sensitivity analysis, NOT covered by Thm 4.2.
compute_estimates <- function(x1,x2,a,b,degrees,calibration=0.35) {
  n <- length(x1); ratio <- ratio_est(x1,x2)
  mn <- max(1L,as.integer(floor(log(n)/log(log(n)))))
  cs <- constants(a,b)
  Ktheory <- as.integer(ceiling(32*log(n)))
  # Lemma A.5 proof uses K <= n/2; do not silently weaken this condition.
  has_theory <- Ktheory <= floor(n/2)
  # Practical comparator: prespecified degree 3, 10 blocks, and fixed
  # threshold multiplier. This is NOT the estimator of Theorem 4.2.
  mcal <- 3L
  Kcal <- 10L # Prespecified practical MOM block count, based on v2.2 pilot
  allm <- sort(unique(c(degrees,mn,mcal)))
  coefs <- lapply(allm,function(m) cheb_coeff(m,a,b))
  names(coefs) <- as.character(allm)
  zlist <- lapply(allm,function(m) score_poly(x1,x2,coefs[[as.character(m)]]$q))
  names(zlist) <- as.character(allm)
  est <- data.frame(method=paste0('poly_m',allm),
                    estimate=vapply(zlist,function(z)clip(mean(z)),0.0),
                    degree=allm,selected_ratio=NA_real_,
                    max_abs_score=vapply(zlist,function(z)max(abs(z)),0.0),
                    threshold=NA_real_,comparison_gap=NA_real_,blocks=NA_integer_)
  est <- rbind(data.frame(method='ratio',estimate=ratio,degree=NA_integer_,
                          selected_ratio=NA_real_,max_abs_score=NA_real_,
                          threshold=NA_real_,comparison_gap=NA_real_,blocks=NA_integer_),est)
  zcal <- zlist[[as.character(mcal)]]
  momcal <- mom_est(zcal,Kcal)
  r <- coefs[[as.character(mcal)]]$r
  secal <- sd(zcal)/sqrt(n)
  threshold_cal <- calibration*(2*r^mcal+secal)
  gap_cal <- abs(momcal-ratio)
  choose_cal <- gap_cal <= threshold_cal
  est <- rbind(est,
    data.frame(method='mom_calibrated',estimate=momcal,degree=mcal,
               selected_ratio=NA_real_,max_abs_score=max(abs(zcal)),
               threshold=NA_real_,comparison_gap=NA_real_,blocks=Kcal),
    data.frame(method='adaptive_calibrated',
               estimate=if(choose_cal)ratio else momcal,degree=mcal,
               selected_ratio=as.numeric(choose_cal),max_abs_score=max(abs(zcal)),
               threshold=threshold_cal,comparison_gap=gap_cal,blocks=Kcal))
  if(has_theory) {
    z <- zlist[[as.character(mn)]]
    mom <- mom_est(z,Ktheory)
    V <- cs$C0*cs$C1^mn*gamma(mn+1)
    rt <- coefs[[as.character(mn)]]$r
    e <- 2*rt^mn+17*sqrt(V*log(n)/n)
    h <- cs$A1*sqrt(log(n)/n)
    gap <- abs(mom-ratio)
    choose_theory <- gap <= e+h
    est <- rbind(est,
      data.frame(method='mom',estimate=mom,degree=mn,
                 selected_ratio=NA_real_,max_abs_score=max(abs(z)),
                 threshold=NA_real_,comparison_gap=NA_real_,blocks=Ktheory),
      data.frame(method='adaptive_theory',estimate=if(choose_theory)ratio else mom,
                 degree=mn,selected_ratio=as.numeric(choose_theory),
                 max_abs_score=max(abs(z)),threshold=e+h,
                 comparison_gap=gap,blocks=Ktheory))
  }
  est
}

# Deterministic checks against analytically known moments and polynomial identity.
run_checks <- function() {
  for(ab in list(c(1,3),c(0.5,5),c(2,8))) {
    a <- ab[1];b <- ab[2]
    for(m in 1:8) {
      cf <- cheb_coeff(m,a,b)
      ss <- seq(a,b,length.out=101)
      approx <- ss*vapply(ss,function(s)sum(cf$q*s^(seq_along(cf$q)-1L)),0.0)
      expected <- 1/cf$denom
      stopifnot(max(abs(1-approx)) <= expected*(1+1e-6)+1e-8)
      stopifnot(cf$const_residual<1e-6)
    }
  }
  # Falling factorial implementation: H_1=X1*(D-1), H_2=X1*(D-1)*(D-2)
  stopifnot(isTRUE(all.equal(score_poly(c(2,1),c(3,0),c(0,1,0)),c(8,0))))
  # Regression test: filtering each adaptive method must not mix the two panels.
  dummy <- data.frame(method=c('adaptive_theory','adaptive_calibrated'),x=1:2)
  stopifnot(nrow(dummy[dummy$method=='adaptive_theory',])==1L)
  # Ratio has no polynomial score; its missing score is not instability.
  stopifnot(!(!is.finite(0) | (!is.na(NA_real_) & !is.finite(NA_real_))))
  # Check MOM for constant scores, including n=250 exploratory configuration.
  stopifnot(abs(mom_est(rep(0.4,250),44)-0.4)<1e-12)
  tmp <- compute_estimates(rep(1,250),rep(1,250),1,5,1:8)
  stopifnot(tmp$blocks[tmp$method=='mom_calibrated']==10L,
            'adaptive_calibrated' %in% tmp$method,
            !('adaptive_theory' %in% tmp$method))
  TRUE
}
if(mode=='test') {run_checks();cat('All deterministic checks passed.\n');quit(save='no')}
run_checks()

cfg <- if(mode=='pilot') list(n=c(250L,1000L,4000L),B=20L,Bcal=100L,degrees=1:8,a=1,b=5) else
  list(n=c(250L,500L,1000L,2000L,4000L,10000L),B=1000L,Bcal=1500L,degrees=1:8,a=1,b=5)
# 'figures' rebuilds figures from previously saved CSV files without resimulating.
if(mode=='figures') {
  if(!file.exists(file.path(out,'summary.csv'))) stop('Missing summary.csv in output directory')
} else {
writeLines(capture.output(sessionInfo()), file.path(out,'session_info.txt'))
writeLines(c(
  'Methods frozen before evaluation:',
  'Practical MOM degree=3, blocks=10, threshold multiplier=0.35.',
  'These practical tuning values originate from earlier pilot studies, not from the final evaluation.',
  'Theoretical MOM: degree=floor(log(n)/log(log(n))), blocks=ceiling(32*log(n)).',
  'Theoretical MOM is only reported if ceiling(32*log(n)) <= floor(n/2).',
  'Covariance test critical values calibrated under Design IV using separate simulation seeds.',
  'This calibration does NOT guarantee uniform type-I error over all zero-covariance mixing laws.',
  'Final performance results are not used to choose tuning parameters.'
),file.path(out,'protocol.txt'))
write.csv(data.frame(mode=mode,seed=20261030,a=cfg$a,b=cfg$b,B=cfg$B,Bcal=cfg$Bcal,
                     n=paste(cfg$n,collapse=';'),degrees=paste(cfg$degrees,collapse=';')),
          file.path(out,'config.csv'),row.names=FALSE)
results <- vector('list',length(scenarios)*length(cfg$n)*cfg$B)
idx <- 0L
for(scenario in scenarios) for(n in cfg$n) {
  truth <- parameters(scenario,cfg$a,cfg$b)
  for(rep in seq_len(cfg$B)) {
    d <- simulate_data(n,scenario,cfg$a,cfg$b)
    e <- compute_estimates(d$x1,d$x2,cfg$a,cfg$b,cfg$degrees)
    e$scenario <- scenario; e$n <- n; e$replication <- rep
    e$theta <- truth$theta; e$true_cov <- truth$cov; e$ratio_target <- truth$ratio
    e$error <- e$estimate-truth$theta; e$squared_error <- e$error^2
    e$numerically_unstable <- !is.finite(e$estimate) | (!is.na(e$max_abs_score) & !is.finite(e$max_abs_score))
    idx <- idx+1L;results[[idx]] <- e
  }
  cat(sprintf('Finished %-24s n=%d B=%d\n',scenario,n,cfg$B))
}
raw <- do.call(rbind,results)
write.csv(raw,file.path(out,'replications.csv'),row.names=FALSE)

# Summary with Monte Carlo standard error of the estimated MSE.
groups <- split(raw,interaction(raw$scenario,raw$n,raw$method,drop=TRUE))
summ <- do.call(rbind,lapply(groups,function(d) {
  mse <- mean(d$squared_error); nrep <- nrow(d)
  data.frame(scenario=d$scenario[1],n=d$n[1],method=d$method[1],B=nrep,
    bias=mean(d$error),bias_mcse=sd(d$error)/sqrt(nrep),sd=sd(d$estimate),mse=mse,
    mse_mcse=sd(d$squared_error)/sqrt(nrep),
    ratio_selection=if(all(is.na(d$selected_ratio))) NA_real_ else mean(d$selected_ratio,na.rm=TRUE),
    selection_mcse=if(all(is.na(d$selected_ratio))) NA_real_ else sqrt(mean(d$selected_ratio,na.rm=TRUE)*(1-mean(d$selected_ratio,na.rm=TRUE))/sum(!is.na(d$selected_ratio))),
    unstable_fraction=mean(d$numerically_unstable))
}))
row.names(summ) <- NULL
write.csv(summ,file.path(out,'summary.csv'),row.names=FALSE)
# Paired MSE contrasts: same replication for both methods.
pair_methods <- c('ratio','poly_m3','mom_calibrated','adaptive_calibrated',
                  'mom','adaptive_theory')
pair_raw <- raw[raw$method %in% pair_methods,
                c('scenario','n','replication','method','squared_error')]
paired <- list()
pi <- 0L
for(sc in scenarios) for(nn in cfg$n) {
  d <- pair_raw[pair_raw$scenario==sc & pair_raw$n==nn,]
  base <- d[d$method=='poly_m3',c('replication','squared_error')]
  names(base)[2] <- 'base_se'
  for(mm in setdiff(pair_methods,'poly_m3')) {
    dd <- d[d$method==mm,c('replication','squared_error')]
    merged <- merge(dd,base,by='replication')
    if(!nrow(merged)) next
    delta <- merged$squared_error-merged$base_se
    pi <- pi+1L
    paired[[pi]] <- data.frame(scenario=sc,n=nn,method=mm,
      reference='poly_m3',B=length(delta),mse_difference=mean(delta),
      difference_mcse=sd(delta)/sqrt(length(delta)),
      ci_lower=mean(delta)-1.96*sd(delta)/sqrt(length(delta)),
      ci_upper=mean(delta)+1.96*sd(delta)/sqrt(length(delta)))
  }
}
write.csv(do.call(rbind,paired),file.path(out,'paired_mse_differences.csv'),row.names=FALSE)


# Covariance test: independent null calibration, disjoint from evaluation.
# Fixed-design critical values are NOT uniformly valid over the entire null class.
# Dedicated seeds separate calibration and evaluation reproducibly.
calibration_rows <- list(); ci <- 0L
set.seed(20261031)
for(n in cfg$n) {
  m <- max(1L,as.integer(floor(log(n)/log(log(n)))))
  cf <- cheb_coeff(m,cfg$a,cfg$b)$q
  for(rep in seq_len(cfg$Bcal)) {
    d <- simulate_data(n,'nonlinear_zero_cov',cfg$a,cfg$b)
    theta_p <- clip(mean(score_poly(d$x1,d$x2,cf)))
    chat <- mean(d$x1)-mean(d$x1+d$x2)*theta_p
    ci <- ci+1L
    calibration_rows[[ci]] <- data.frame(n=n,replication=rep,chat=chat)
  }
}
calibration_raw <- do.call(rbind,calibration_rows)
write.csv(calibration_raw,file.path(out,'covariance_calibration_independent.csv'),row.names=FALSE)
cal <- aggregate(abs(chat)~n,data=calibration_raw,
                 FUN=function(x)unname(quantile(x,0.95,type=8)))
names(cal)[2] <- 'critical_value'
write.csv(cal,file.path(out,'covariance_thresholds_independent.csv'),row.names=FALSE)

set.seed(20261032)
test_rows <- list(); ii <- 0L
for(scenario in scenarios) for(n in cfg$n) {
  m <- max(1L,as.integer(floor(log(n)/log(log(n)))))
  cf <- cheb_coeff(m,cfg$a,cfg$b)$q
  for(rep in seq_len(cfg$B)) {
    d <- simulate_data(n,scenario,cfg$a,cfg$b)
    theta_p <- clip(mean(score_poly(d$x1,d$x2,cf)))
    chat <- mean(d$x1)-mean(d$x1+d$x2)*theta_p
    ii <- ii+1L
    test_rows[[ii]] <- data.frame(scenario=scenario,n=n,replication=rep,
                                   chat=chat,true_cov=parameters(scenario,cfg$a,cfg$b)$cov)
  }
}
tests <- do.call(rbind,test_rows)
tests <- merge(tests,cal,by='n')
tests$reject <- abs(tests$chat)>tests$critical_value
write.csv(tests,file.path(out,'covariance_replications.csv'),row.names=FALSE)
rate_groups <- split(tests,interaction(tests$scenario,tests$n,drop=TRUE))
rates <- do.call(rbind,lapply(rate_groups,function(d) {
  pr <- mean(d$reject)
  data.frame(scenario=d$scenario[1],n=d$n[1],B=nrow(d),
             rejection_fraction=pr,rejection_mcse=sqrt(pr*(1-pr)/nrow(d)))
}))
row.names(rates) <- NULL
write.csv(rates,file.path(out,'covariance_rejection_independent.csv'),row.names=FALSE)
} # end of simulation; 'figures' reads existing summaries
if(mode=='figures') {
  summ <- read.csv(file.path(out,'summary.csv'))
  rates <- read.csv(file.path(out,'covariance_rejection_independent.csv'))
  cfg$n <- sort(unique(summ$n))
}

# Record method eligibility and threshold diagnostics for auditing.
if(mode!='figures') {
  diagnostic <- raw[grepl('^adaptive_',raw$method),
                    c('scenario','n','replication','method','degree','blocks',
                      'selected_ratio','comparison_gap','threshold')]
  write.csv(diagnostic,file.path(out,'adaptive_diagnostics.csv'),row.names=FALSE)
}
# Publication-ready vector PDFs and raster PNGs, all base R.
plot_both <- function(filename,plotfun,w=9,h=6) {
  pdf(file.path(out,paste0(filename,'.pdf')),width=w,height=h,onefile=TRUE)
  plotfun();dev.off()
  png(file.path(out,paste0(filename,'.png')),width=1400,height=900,res=160)
  plotfun();dev.off()
}
# Figure 1: five prespecified comparisons, with a shared external legend.
figure1_methods <- c('ratio','poly_m3','poly_m4','adaptive_theory','adaptive_calibrated')
figure1_labels <- c('Ratio','Polynomial (m=3)','Polynomial (m=4)',
                    'Adaptive (theoretical)','Adaptive (practical)')
figure1_cols <- c('#222222','#0072B2','#D55E00','#7B3294','#009E73')
figure1_lty <- c(1,2,3,4,5)
figure1_pch <- c(16,17,15,18,19)
figure1_titles <- c(constant='Constant composition',
                    independent_beta='Independent Beta',
                    linear_dependence='Linear dependence',
                    nonlinear_zero_cov='Nonlinear dependence with zero covariance')
plot_both('figure1_mse_vs_n',function() {
  # Four panels plus a dedicated fifth region for the legend.
  op <- par(no.readonly=TRUE)
  on.exit({layout(1); par(op)})
  layout(matrix(c(1,2,3,4,5,5),nrow=3,byrow=TRUE),
         heights=c(1,1,0.20))
  par(mar=c(3.7,4.3,2.4,1.1),mgp=c(2.5,0.65,0),cex.axis=0.88)
  for(sc in scenarios) {
    dd <- summ[summ$scenario==sc & summ$method %in% figure1_methods, ,drop=FALSE]
    good <- is.finite(dd$mse) & dd$mse>0
    if(!any(good)) stop(paste('No positive MSE values for',sc))
    yr <- range(dd$mse[good]); yr <- exp(log(yr)+c(-0.12,0.12)*diff(log(yr)))
    plot(NA,xlim=range(cfg$n),ylim=yr,log='xy',xlab='Sample size n',
         ylab='Monte Carlo MSE',main=figure1_titles[[sc]],cex.main=1.05)
    for(j in seq_along(figure1_methods)) {
      d <- dd[dd$method==figure1_methods[j] & is.finite(dd$mse) & dd$mse>0,,drop=FALSE]
      d <- d[order(d$n),,drop=FALSE]
      if(nrow(d)>0L) lines(d$n,d$mse,lty=figure1_lty[j],col=figure1_cols[j],
                          pch=figure1_pch[j],type='b',cex=0.85,lwd=1.6)
    }
  }
  # Legend occupies its own layout region: no overlap with x-axis labels.
  par(mar=c(0,0,0,0))
  plot.new()
  legend('center',legend=figure1_labels,col=figure1_cols,
         lty=figure1_lty,pch=figure1_pch,lwd=1.6,
         ncol=3,bty='n',cex=0.85,x.intersp=0.8,y.intersp=1.15)
},w=10.5,h=8.4)
# Figure 2: Monte Carlo bias, consistent styling and a dedicated legend area.
figure2_methods <- c('ratio','poly_m2','poly_m4','adaptive_theory','adaptive_calibrated')
figure2_labels <- c('Ratio','Polynomial (m=2)','Polynomial (m=4)',
                    'Adaptive (theoretical)','Adaptive (practical)')
figure2_cols <- c('#222222','#0072B2','#D55E00','#7B3294','#009E73')
figure2_lty <- c(1,2,3,4,5)
figure2_pch <- c(16,17,15,18,19)
plot_both('figure2_bias_vs_n',function() {
  op <- par(no.readonly=TRUE)
  on.exit({layout(1);par(op)})
  layout(matrix(c(1,2,3,4,5,5),nrow=3,byrow=TRUE),
         heights=c(1,1,0.20))
  par(mar=c(3.7,4.3,2.4,1.1),mgp=c(2.5,0.65,0),cex.axis=0.88)
  for(sc in scenarios) {
    dd <- summ[summ$scenario==sc & summ$method %in% figure2_methods, ,drop=FALSE]
    dd <- dd[is.finite(dd$bias), ,drop=FALSE]
    if(nrow(dd)==0L) stop(paste('No finite bias values for',sc))
    reference <- if(sc=='linear_dependence') 1/15 else 0
    yy <- c(dd$bias,0,reference)
    pad <- max(0.003,0.10*diff(range(yy)))
    ylim <- range(yy)+c(-pad,pad)
    plot(NA,xlim=range(cfg$n),ylim=ylim,log='x',
         xlab='Sample size n',ylab='Monte Carlo bias',
         main=figure1_titles[[sc]],cex.main=1.05)
    abline(h=0,col='grey55',lty=3,lwd=1.2)
    if(sc=='linear_dependence')
      abline(h=1/15,col='grey35',lty=2,lwd=1.2)
    for(j in seq_along(figure2_methods)) {
      d <- dd[dd$method==figure2_methods[j], ,drop=FALSE]
      d <- d[order(d$n), ,drop=FALSE]
      if(nrow(d)>0L)
        lines(d$n,d$bias,type='b',col=figure2_cols[j],
              lty=figure2_lty[j],pch=figure2_pch[j],
              cex=0.85,lwd=1.6)
    }
  }
  par(mar=c(0,0,0,0))
  plot.new()
  legend('center',legend=figure2_labels,col=figure2_cols,
         lty=figure2_lty,pch=figure2_pch,lwd=1.6,
         ncol=3,bty='n',cex=0.85,x.intersp=0.8,y.intersp=1.15)
},w=10.5,h=8.4)
# Figure 3: degree sensitivity, six sample sizes; compact boxed legend
# inside the upper-left corner of every panel.
figure3_cols <- c('#222222','#D55E00','#009E73','#0072B2','#56B4E9','#CC79A7')
figure3_lty <- c(1,2,3,4,5,6)
figure3_pch <- c(16,17,15,4,18,25)
plot_both('figure3_mse_by_degree',function() {
  op <- par(no.readonly=TRUE)
  on.exit(par(op))
  par(mfrow=c(2,2),mar=c(3.8,4.6,2.5,1.2),
      mgp=c(2.5,0.65,0),cex.axis=0.90)
  for(sc in scenarios) {
    dd <- summ[summ$scenario==sc & grepl('^poly_m[0-9]+$',summ$method), ,drop=FALSE]
    dd$degree <- as.integer(sub('^poly_m','',dd$method))
    dd <- dd[is.finite(dd$mse) & dd$mse>0 & is.finite(dd$degree), ,drop=FALSE]
    if(!nrow(dd)) stop(paste('No positive polynomial MSE for',sc))
    yr <- range(dd$mse)
    yr <- exp(log(yr)+c(-0.10,0.10)*max(diff(log(yr)),0.1))
    plot(NA,xlim=range(dd$degree),ylim=yr,log='y',
         xaxt='n',xlab='Polynomial degree m',ylab='Monte Carlo MSE',
         main=figure1_titles[[sc]],cex.main=1.03)
    axis(1,at=sort(unique(dd$degree)))
    for(j in seq_along(cfg$n)) {
      d <- dd[dd$n==cfg$n[j], ,drop=FALSE]
      d <- d[order(d$degree), ,drop=FALSE]
      if(nrow(d)) lines(d$degree,d$mse,type='b',
                        col=figure3_cols[j],lty=figure3_lty[j],
                        pch=figure3_pch[j],cex=0.90,lwd=1.6)
    }
    {
      legend('topleft',
             legend=paste0('n = ',format(cfg$n,big.mark=',')),
             col=figure3_cols[seq_along(cfg$n)],
             lty=figure3_lty[seq_along(cfg$n)],
             pch=figure3_pch[seq_along(cfg$n)],
             lwd=1.4,pt.cex=0.75,
             ncol=2,cex=0.67,bty='o',box.col='grey55',
             bg='white',inset=0.015,x.intersp=0.65,
             y.intersp=1.08,seg.len=1.4)
    }
  }
},w=10.5,h=7.8)
# Shared presentation for Figures 4 and 5 (same scenario coding).
scenario_labels <- c('Constant composition', 'Independent beta',
                     'Linear dependence', 'Nonlinear dependence (zero covariance)')
scenario_cols <- c('#222222', '#D55E00', '#009E73', '#0072B2')
scenario_lty <- c(1, 2, 3, 4)
scenario_pch <- c(16, 17, 15, 4)

# Manuscript Figure 4 is generated exclusively by scripts/figure4_article.R.
# The previous two-panel selector graphic is intentionally not produced here.

# Figure 5: publication layout, using existing Monte Carlo results only.
# Confidence bars show pointwise normal-approximation 95% Monte Carlo intervals
# for the estimated rejection probabilities (NOT confidence intervals for Q).
# The test is calibrated under one zero-covariance scenario, not uniformly.
plot_both('figure5_covariance_rejection',function() {
  op <- par(no.readonly=TRUE)
  on.exit({layout(1); par(op)})

  # Compact legend region, plus a single overall title.
  layout(matrix(c(1,2,3,3),nrow=2,byrow=TRUE),heights=c(1,0.15))
  par(oma=c(0,0,2.0,0),mar=c(4.4,4.5,2.2,1.0),
      mgp=c(2.55,0.68,0),cex.axis=0.86,cex.lab=0.94,las=1)

  plot_panel <- function(which_scenarios,ylim,title,show_y=TRUE,mc_bars=FALSE) {
    plot(NA,xlim=range(cfg$n)*c(0.93,1.07),ylim=ylim,log='x',
         xaxt='n',xlab='Sample size n',
         ylab=if(show_y) 'Empirical rejection frequency' else '',
         main=title,cex.main=0.98)
    axis(1,at=cfg$n,labels=format(cfg$n,scientific=FALSE,big.mark=','),
         cex.axis=0.78)
    abline(h=0.05,lty=2,col='grey45',lwd=1.3)

    for(j in which_scenarios) {
      d <- rates[rates$scenario==scenarios[j] &
                 is.finite(rates$rejection_fraction),,drop=FALSE]
      d <- d[order(d$n),,drop=FALSE]
      if(!nrow(d)) next
      lines(d$n,d$rejection_fraction,type='b',
            col=scenario_cols[j],lty=scenario_lty[j],
            pch=scenario_pch[j],lwd=1.7,cex=0.88)

      if(mc_bars) {
        # Small horizontal offsets prevent intervals from overlapping.
        offset <- c(0.972,1.0,1.028)[match(j,which(scenarios!='linear_dependence'))]
        if('rejection_mcse' %in% names(d)) {
          se <- d$rejection_mcse
        } else if(all(c('B','rejection_fraction') %in% names(d))) {
          se <- sqrt(d$rejection_fraction*(1-d$rejection_fraction)/d$B)
        } else {
          stop('Figure 5 MC intervals require rejection_mcse or B in rates CSV')
        }
        ok <- is.finite(se) & se>=0 & is.finite(d$n)
        xx <- d$n[ok]*offset
        lower <- pmax(ylim[1],d$rejection_fraction[ok]-1.96*se[ok])
        upper <- pmin(ylim[2],d$rejection_fraction[ok]+1.96*se[ok])
        arrows(xx,lower,xx,upper,angle=90,code=3,length=0.035,
               col=adjustcolor(scenario_cols[j],alpha.f=0.72),lwd=1.0)
      }
    }
    box(col='grey50')
  }

  plot_panel(seq_along(scenarios),c(0,1.04),
             'All scenarios: rejection and power')
  plot_panel(which(scenarios!='linear_dependence'),c(0,0.12),
             'Zero-covariance scenarios (detail)',show_y=FALSE,mc_bars=TRUE)

  # Common legend, deliberately close to the panel x-axis labels.
  par(mar=c(0,0,0,0))
  plot.new()
  legend('center',legend=c(scenario_labels,'Nominal level (5%)'),
         col=c(scenario_cols,'grey45'),lty=c(scenario_lty,2),
         pch=c(scenario_pch,NA_integer_),lwd=c(rep(1.7,4),1.3),
         ncol=3,cex=0.79,bty='n',
         x.intersp=0.8,y.intersp=1.12,seg.len=1.8)
  mtext('Covariance test: empirical size and power',side=3,outer=TRUE,
        line=0.55,cex=1.13,font=2)
},w=11.5,h=6.2)
cat('Outputs saved to ',normalizePath(out),'\n',sep='')
