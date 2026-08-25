suppressMessages(devtools::load_all(quiet = TRUE))
source(file.path("dev", "analytic_certification", "certify_analytic.R"))

# 1) re-validate v2 on the full grids used before
check <- function(l,u,n,d, step=NULL) {
  W <- u-l
  means <- round(seq(l, u, by = 10^-d), d)
  sds <- round(seq(0, (W/2)*sqrt(n/(n-1))*1.05, by = 10^-d), d)
  g <- expand.grid(mean=means, sd=sds)
  if (!is.null(step)) g <- g[seq(1,nrow(g),by=step),]
  t0 <- Sys.time(); ref <- brimmest(l,u,n,mean=g$mean,sd=g$sd,digits=d)
  tb <- as.numeric(Sys.time()-t0, units="secs")
  t0 <- Sys.time(); v2 <- certify_analytic2(l,u,n,g$mean,g$sd,d)
  ta <- as.numeric(Sys.time()-t0, units="secs")
  cat(sprintf("grid %d-%d n=%d d=%d: %d tuples | brimmest %.2fs, analytic-v2 %.2fs | disagreements %d\n",
              l,u,n,d,nrow(g),tb,ta,sum(ref$possible!=v2$possible)))
  sum(ref$possible!=v2$possible)
}
tot <- 0
tot <- tot + check(1,5,9,1)
tot <- tot + check(1,7,9,1)
tot <- tot + check(0,4,15,2, step=7)

# 2) edge-tuple stress incl. the previous worst cases
set.seed(7)
stress <- function(l,u,n,d,N=200) {
  W <- u-l
  sd_max_of <- function(m){S<-round((m-l)*n);a<-S%/%W;b<-S-a*W;sqrt((a*W^2+b^2-S^2/n)/(n-1))}
  means <- round(runif(N, l+0.02*(u-l), u-0.02*(u-l)), d)
  sds <- vapply(means, function(m) round(max(0, sd_max_of(m)-runif(1,0,2)*10^-d), d), 0)
  t_b <- t_a <- numeric(N); ag <- logical(N)
  for (i in seq_len(N)) {
    t0<-Sys.time(); rb<-brimmest(l,u,n,mean=means[i],sd=sds[i],digits=d); t_b[i]<-as.numeric(Sys.time()-t0,units="secs")
    t0<-Sys.time(); ra<-certify_analytic2(l,u,n,means[i],sds[i],d);      t_a[i]<-as.numeric(Sys.time()-t0,units="secs")
    ag[i] <- rb$possible==ra$possible
  }
  cat(sprintf("edge %d-%d n=%d d=%d: brimmest med %.1f max %.0f ms | v2 med %.2f max %.1f ms | disagreements %d\n",
              l,u,n,d, 1000*median(t_b),1000*max(t_b), 1000*median(t_a),1000*max(t_a), sum(!ag)))
  sum(!ag)
}
tot <- tot + stress(0,63,50,2)
tot <- tot + stress(0,63,200,2)
tot <- tot + stress(1,7,2000,2)
# the previous 7.9s outlier:
t0<-Sys.time(); r<-certify_analytic2(0,63,200,27.10,31.25,2); t1<-as.numeric(Sys.time()-t0,units="secs")
cat(sprintf("previous 7.9s outlier (0-63 n=200 mean 27.10 sd 31.25): v2 %.1f ms, possible=%s\n", 1000*t1, r$possible))
cat(sprintf("TOTAL disagreements: %d\n", tot))
