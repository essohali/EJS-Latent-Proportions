#!/usr/bin/env Rscript
# Reproduce and audit Table 4: support-bound sensitivity of degree-3 estimator.
# Usage: Rscript scripts/table4_support_sensitivity.R [full|test] [output_directory]
# Run from repository root. No external R packages required.
args <- commandArgs(trailingOnly=TRUE)
mode <- if(length(args)) args[1] else 'test'
stopifnot(mode %in% c('test','full'))
out <- if(length(args)>=2) args[2] else file.path('results','generated','table4_support_sensitivity',mode)
dir.create(out,recursive=TRUE,showWarnings=FALSE)

# Import EXACT function definitions from the main simulation script, without
# sourcing it (sourcing would run its simulation pipeline and possibly quit R).
main <- file.path('scripts','ejs_simulations.R')
if(!file.exists(main)) stop('Run this script from the repository root: ',main)
needed <- c('poly_mul_linear','cheb_coeff','score_poly','clip','simulate_data')
exprs <- parse(main,keep.source=FALSE)
for(e in exprs) {
  if(is.call(e) && identical(e[[1]],as.name('<-')) && is.symbol(e[[2]]) &&
     as.character(e[[2]]) %in% needed && is.call(e[[3]]) && identical(e[[3]][[1]],as.name('function'))) {
    eval(e,envir=environment())
  }
}
if(!all(vapply(needed, function(nm) exists(nm, envir=.GlobalEnv, mode='function', inherits=FALSE), logical(1))))
  stop('Cannot load all required functions from ',main)

n <- if(mode=='full') 10000L else 500L
B <- if(mode=='full') 1000L else 20L
m <- 3L
seed <- 20261108L
set.seed(seed)
# True support and polynomial working support are separate by design.
cases <- data.frame(
  configuration=c('correct_2_5','correct_1_5','correct_0.2_5',
                  'working_1.5_4.5','working_0.2_6'),
  true_a=c(2,1,0.2,1,1), true_b=c(5,5,5,5,5),
  working_a=c(2,1,0.2,1.5,0.2), working_b=c(5,5,5,4.5,6),
  stringsAsFactors=FALSE)
scenarios <- c('constant','linear_dependence')
# Historical manuscript MSEs (rounded to 6 decimal places), not reference truth.
historical <- matrix(c(0.000017,0.000027,
                       0.000049,0.000122,
                       0.000330,0.001723,
                       0.000032,0.000105,
                       0.000223,0.004669),ncol=2,byrow=TRUE,
                     dimnames=list(cases$configuration,scenarios))

rows <- list(); errors <- list(); k <- 0L
for(i in seq_len(nrow(cases))) {
  cs <- cases[i,]
  coeff <- cheb_coeff(m,cs$working_a,cs$working_b)$q
  for(sc in scenarios) {
    squared_errors <- numeric(B)
    estimates <- numeric(B)
    for(j in seq_len(B)) {
      d <- simulate_data(n,sc,cs$true_a,cs$true_b)
      estimates[j] <- clip(mean(score_poly(d$x1,d$x2,coeff)))
      squared_errors[j] <- (estimates[j]-0.4)^2
    }
    k <- k+1L
    mse <- mean(squared_errors)
    mcse <- sd(squared_errors)/sqrt(B)
    ref <- historical[cs$configuration,sc]
    # A six-decimal historical number is rounded, with rounding error <= 5e-7.
    # 95% CI is an approximate Monte Carlo CI for a fresh, independent run.
    rows[[k]] <- data.frame(configuration=cs$configuration,scenario=sc,
      true_a=cs$true_a,true_b=cs$true_b,
      working_a=cs$working_a,working_b=cs$working_b,
      n=n,B=B,m=m,seed=seed,mse=mse,mcse=mcse,
      ci_lower=max(0,mse-1.96*mcse),ci_upper=mse+1.96*mcse,
      historical_mse=ref,difference=mse-ref,
      historical_in_mc_ci=(ref>=mse-1.96*mcse-5e-7 && ref<=mse+1.96*mcse+5e-7))
    errors[[k]] <- data.frame(configuration=cs$configuration,scenario=sc,
      replication=seq_len(B),estimate=estimates,squared_error=squared_errors)
    cat(sprintf('%-18s %-18s mse=%.9g  MCSE=%.3g  historical=%.9g\n',
      cs$configuration,sc,mse,mcse,ref))
  }
}
result <- do.call(rbind,rows)
write.csv(result,file.path(out,'table4_audit.csv'),row.names=FALSE)
write.csv(do.call(rbind,errors),file.path(out,'table4_replications.csv'),row.names=FALSE)
writeLines(c(sprintf('Mode: %s; n=%d; B=%d; degree=%d; seed=%d',mode,n,B,m,seed),
  'Estimator: clip(mean(score_poly(X1,X2,cheb_coeff(3,working_a,working_b)$q))).',
  'True target theta=0.4 in both scenarios.',
  'Functions imported verbatim at runtime from scripts/ejs_simulations.R.',
  'Historical manuscript MSEs are rounded to six decimal places.',
  '95% Monte Carlo intervals use normal approximation based on squared-error SD.',
  'Matching is statistical, not bitwise, because original random stream is unknown.',
  'Changing the true support also changes the data-generating law.'),
  file.path(out,'protocol.txt'))
writeLines(capture.output(sessionInfo()),file.path(out,'session_info.txt'))
cat('Saved results in ',out,'\n',sep='')
