#!/usr/bin/env Rscript
# Rebuild manuscript Figure 4 from the original replication-level Monte Carlo results.
# Usage: Rscript scripts/figure4_article.R path/to/replications.csv output/directory
args <- commandArgs(trailingOnly=TRUE)
input <- if(length(args)>=1L) args[1] else 'results/reference/replications.csv'
out <- if(length(args)>=2L) args[2] else 'results/generated/figure4'
stopifnot(file.exists(input))
dir.create(out,recursive=TRUE,showWarnings=FALSE)
set.seed(20261108)
Bboot <- 2000L
raw <- read.csv(input,stringsAsFactors=FALSE)
required <- c('scenario','n','replication','method','selected_ratio','squared_error')
stopifnot(all(required %in% names(raw)))
scenarios <- c('constant','independent_beta','linear_dependence','nonlinear_zero_cov')
ns <- sort(unique(raw$n))
rows <- list(); k <- 0L
for(sc in scenarios) for(nn in ns) {
  a <- raw[raw$scenario==sc & raw$n==nn & raw$method=='adaptive_calibrated',c('replication','selected_ratio','squared_error')]
  p <- raw[raw$scenario==sc & raw$n==nn & raw$method=='poly_m3',c('replication','squared_error')]
  names(a)[3] <- 'se_ad'; names(p)[2] <- 'se_p3'
  d <- merge(a,p,by='replication',sort=TRUE)
  if(nrow(d)==0L || anyNA(d[,c('selected_ratio','se_ad','se_p3')])) stop(paste('Missing/invalid data',sc,nn))
  if(anyDuplicated(d$replication)) stop(paste('Duplicate replications',sc,nn))
  B <- nrow(d)
  sel <- mean(d$selected_ratio)
  ratio <- mean(d$se_ad)/mean(d$se_p3)
  # Paired nonparametric bootstrap: same resampled replication indices for both estimators.
  boot_sel <- numeric(Bboot); boot_ratio <- numeric(Bboot)
  for(j in seq_len(Bboot)) {
    ix <- sample.int(B,B,replace=TRUE)
    boot_sel[j] <- mean(d$selected_ratio[ix])
    boot_ratio[j] <- mean(d$se_ad[ix])/mean(d$se_p3[ix])
  }
  qsel <- unname(quantile(boot_sel,c(.025,.975),type=7))
  qrat <- unname(quantile(boot_ratio,c(.025,.975),type=7))
  k <- k+1L
  rows[[k]] <- data.frame(scenario=sc,n=nn,B=B,ratio_selection=sel,selection_lo=qsel[1],selection_hi=qsel[2],relative_mse=ratio,relative_mse_lo=qrat[1],relative_mse_hi=qrat[2])
}
res <- do.call(rbind,rows)
write.csv(res,file.path(out,'figure4_audit.csv'),row.names=FALSE)
labels <- c('Constant composition','Independent beta','Linear dependence','Nonlinear dependence (zero covariance)')
cols <- c('#222222','#D55E00','#009E73','#0072B2')
ltys <- c(1,2,2,2); pchs <- c(16,17,15,3)
# Equal spacing of the six sample sizes, as in the manuscript figure.
plot_fig <- function() {
  op <- par(no.readonly=TRUE)
  on.exit({layout(1); par(op)})
  # Third plotting region is dedicated entirely to the legend.
  layout(matrix(c(1,2,3,3),nrow=2,byrow=TRUE), heights=c(1,0.23))
  x <- seq_along(ns)
  par(mar=c(4.0,4.2,2.1,0.8), mgp=c(2.35,0.65,0),
      cex.axis=0.82,cex.lab=0.95,cex.main=0.95,tcl=-0.25)
  plot(NA,xlim=c(0.75,length(ns)+0.25),ylim=c(-0.05,1.05),
       xlab='Sample size n',ylab='Frequency of selecting the ratio estimator',
       main='(A) Practical selection',xaxt='n',bty='l')
  axis(1,at=x,labels=ns)
  for(j in seq_along(scenarios)) {
    d <- res[res$scenario==scenarios[j],]; d <- d[order(d$n),]
    lines(x,d$ratio_selection,type='b',col=cols[j],lty=ltys[j],
          pch=pchs[j],lwd=1.15,cex=0.65)
    ok <- is.finite(d$selection_lo) & is.finite(d$selection_hi) &
          (d$selection_hi > d$selection_lo)
    if(any(ok)) arrows(x[ok],d$selection_lo[ok],x[ok],d$selection_hi[ok],
                       angle=90,code=3,length=0.018,lwd=0.7,col=cols[j])
  }
  par(mar=c(4.0,4.2,2.1,0.8), mgp=c(2.35,0.65,0),
      cex.axis=0.82,cex.lab=0.95,cex.main=0.95,tcl=-0.25)
  plot(NA,xlim=c(0.75,length(ns)+0.25),ylim=c(0.05,2.12),
       xlab='Sample size n',
       
       ylab = expression(
         MSE(hat(theta)[n]^{"ad,pr"}) /
           MSE(hat(theta)[n,3])
       )
       ,
       main='(B) Relative mean squared error',xaxt='n',bty='l')
  axis(1,at=x,labels=ns)
  abline(h=1,lty=2,col='grey60',lwd=0.8)
  for(j in seq_along(scenarios)) {
    d <- res[res$scenario==scenarios[j],]; d <- d[order(d$n),]
    lines(x,d$relative_mse,type='b',col=cols[j],lty=ltys[j],
          pch=pchs[j],lwd=1.15,cex=0.65)
    ok <- is.finite(d$relative_mse_lo) & is.finite(d$relative_mse_hi) &
          (d$relative_mse_hi > d$relative_mse_lo)
    if(any(ok)) arrows(x[ok],d$relative_mse_lo[ok],x[ok],d$relative_mse_hi[ok],
                       angle=90,code=3,length=0.018,lwd=0.7,col=cols[j])
  }
  # A full-width legend, safely below both x-axis titles.
  par(mar=c(0,0,0,0));plot.new()
  legend('center',legend=labels,col=cols,lty=ltys,pch=pchs,
         lwd=1.15,pt.cex=0.7,ncol=2,bty='n',cex=0.78,
         x.intersp=0.7,y.intersp=1.15)
}
pdf(file.path(out,'figure4.pdf'),width=11,height=6.5)
plot_fig();dev.off()
png(file.path(out,'figure4.png'),width=1650,height=975,res=150)
plot_fig();dev.off()
cat('Saved figure and audit CSV in:',normalizePath(out),'\n')
