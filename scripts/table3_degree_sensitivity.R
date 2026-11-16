#!/usr/bin/env Rscript
# Complementary paired experiment for Table 3: n=10000, m=3 versus m=4.
# Base R only; does not overwrite manuscript reference results.
args <- commandArgs(trailingOnly=TRUE)
out <- if (length(args)>=1) args[1] else 'results/table3_degree_sensitivity'
dir.create(out, recursive=TRUE, showWarnings=FALSE)
n <- 10000L; B <- 1000L; Bcal <- 1500L; alpha <- 0.05
set.seed(20261108)
a <- 1; b <- 5; u0 <- 0.4
kappas <- c(0.5,0.75,1,1.25); degrees <- c(3L,4L)
poly_mul_linear <- function(p,z,k) {
  q <- numeric(length(p)+1L)
  q[seq_along(p)] <- q[seq_along(p)] + z*p
  q[seq_along(p)+1L] <- q[seq_along(p)+1L] - k*p
  q
}
cheb_coeff <- function(m,a,b) {
  z <- (a+b)/(b-a); k <- 2/(b-a)
  t0 <- 1; t1 <- c(z,-k)
  if(m==1L) tm <- t1 else {
    for(j in 2:m) {
      prod <- 2*poly_mul_linear(t1,z,k)
      prod[seq_along(t0)] <- prod[seq_along(t0)] - t0
      t0 <- t1; t1 <- prod
    }
    tm <- t1
  }
  -(tm/cosh(m*acosh(z)))[-1]
}
coeffs <- lapply(degrees,cheb_coeff,a=a,b=b)
simulate <- function(n,gamma,scenario='linear') {
  s <- runif(n,a,b); v <- (s-a)/(b-a)
  u <- switch(scenario,
    linear=u0+gamma*(s-(a+b)/2)/(b-a),
    nonlinear=u0+1.2*((v-0.5)^2-1/12))
  stopifnot(all(u>=0 & u<=1))
  list(x1=rpois(n,s*u),x2=rpois(n,s*(1-u)))
}
clip <- function(x) pmax(0,pmin(1,x))
score_poly <- function(x1,x2,coeff) {
  d <- x1+x2; ans <- numeric(length(d)); falling <- rep(1,length(d))
  for(j in 0:(length(coeff)-1L)) {
    if(j>0) falling <- falling*(d-j)
    ans <- ans+coeff[j+1L]*x1*falling
  }
  ans
}
statistic <- function(d,coeff) {
  total <- d$x1+d$x2
  mean(d$x1)-mean(total)*clip(mean(score_poly(d$x1,d$x2,coeff)))
}
# Paired calibration: same null datasets for m=3 and m=4, separate cutoffs.
cal <- matrix(NA_real_,nrow=Bcal,ncol=length(degrees))
for(i in seq_len(Bcal)) {
  dat <- simulate(n,0,'nonlinear')
  for(j in seq_along(degrees)) cal[i,j] <- abs(statistic(dat,coeffs[[j]]))
}
cutoffs <- apply(cal,2,quantile,probs=1-alpha,type=1)
write.csv(data.frame(degree=degrees,cutoff=cutoffs,Bcal=Bcal,alpha=alpha),
          file.path(out,'cutoffs.csv'),row.names=FALSE)
# Independent holdout evaluation of rejection frequencies under the calibration null.
# This evaluates pointwise empirical size only, not uniform size control.
Bnull <- 2000L
null_reject <- matrix(FALSE,nrow=Bnull,ncol=length(degrees))
for(i in seq_len(Bnull)) {
  dat <- simulate(n,0,'nonlinear')
  for(j in seq_along(degrees))
    null_reject[i,j] <- abs(statistic(dat,coeffs[[j]])) > cutoffs[j]
}
null_summary <- data.frame(
  degree=degrees, n=n, Bnull=Bnull,
  rejections=colSums(null_reject),
  rejection_rate=colMeans(null_reject),
  mcse=sqrt(colMeans(null_reject)*(1-colMeans(null_reject))/Bnull),
  cutoff=as.numeric(cutoffs))
write.csv(null_summary,file.path(out,'null_size_check.csv'),row.names=FALSE)
write.csv(data.frame(replication=seq_len(Bnull),
  reject_m3=as.integer(null_reject[,1]),
  reject_m4=as.integer(null_reject[,2])),
  file.path(out,'null_paired_rejections.csv'),row.names=FALSE)
cat('Independent null rejection checks:\n')
print(null_summary)

# Paired alternatives: one sample per replicate and kappa, evaluated at both degrees.
lambda <- log((sqrt(b)+sqrt(a))/(sqrt(b)-sqrt(a)))
raw <- list(); summary <- list(); idx <- 0L
for(kappa in kappas) {
  rho <- exp(-kappa*lambda*log(n)/log(log(n)))
  gamma <- 12*rho/(b-a)
  stopifnot(gamma <= 2*min(u0,1-u0))
  reject <- matrix(FALSE,nrow=B,ncol=length(degrees))
  for(i in seq_len(B)) {
    dat <- simulate(n,gamma,'linear')
    for(j in seq_along(degrees))
      reject[i,j] <- abs(statistic(dat,coeffs[[j]])) > cutoffs[j]
  }
  for(j in seq_along(degrees)) {
    idx <- idx+1L
    summary[[idx]] <- data.frame(n=n,kappa=kappa,rho=rho,gamma=gamma,
      degree=degrees[j],B=B,rejections=sum(reject[,j]),power=mean(reject[,j]),
      mcse=sqrt(mean(reject[,j])*(1-mean(reject[,j]))/B))
  }
  raw[[length(raw)+1L]] <- data.frame(kappa=kappa,replication=seq_len(B),
    reject_m3=as.integer(reject[,1]),reject_m4=as.integer(reject[,2]))
  cat('Completed kappa =',kappa,'\n')
}
write.csv(do.call(rbind,summary),file.path(out,'degree_power.csv'),row.names=FALSE)
write.csv(do.call(rbind,raw),file.path(out,'paired_rejections.csv'),row.names=FALSE)
paired <- do.call(rbind,raw)
paired_summary <- do.call(rbind,lapply(split(paired,paired$kappa),function(z) {
  delta <- z$reject_m3-z$reject_m4
  data.frame(kappa=z$kappa[1],B=nrow(z),
    power_m3=mean(z$reject_m3),power_m4=mean(z$reject_m4),
    difference=mean(delta),paired_mcse=sd(delta)/sqrt(nrow(z)),
    ci_lower=mean(delta)-1.96*sd(delta)/sqrt(nrow(z)),
    ci_upper=mean(delta)+1.96*sd(delta)/sqrt(nrow(z)))
}))
write.csv(paired_summary,file.path(out,'paired_power_differences.csv'),row.names=FALSE)
cat('Paired power differences:\n')
print(paired_summary)
writeLines(c('Paired degree sensitivity, n=10000, degrees=3,4',
  'B=1000, Bcal=1500, Bnull=2000, alpha=0.05, seed=20261108',
  'Same simulated dataset evaluated at both degrees per replicate.',
  'Separate critical values by degree; calibrated under nonlinear covariance-null design.',
  'Independent holdout samples under the same null evaluate pointwise rejection rates.',
  'Local covariance rho_n=exp(-kappa*lambda*log(n)/log(log(n))).',
  'Not a uniform-size calibration; does not prove minimax boundary.'),
  file.path(out,'protocol.txt'))
writeLines(capture.output(sessionInfo()),file.path(out,'session_info.txt'))
print(do.call(rbind,summary))
