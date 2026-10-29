cols <- function(d,cells) { x <- match(cells,d$cells);stopifnot(!anyNA(x));x }
select_features <- function(d,donor,comparison,visibility_other,pool) {
  nd <- rowSums(d$observed[,donor,drop=FALSE])[d$pos]
  no <- rowSums(d$observed[,visibility_other,drop=FALSE])[d$pos]
  visible <- 1-(length(donor)+length(visibility_other)-nd-no)/(length(donor)+length(visibility_other))>=.1 &
    1-(length(donor)-nd)/length(donor)>=.1 & 1-(length(visibility_other)-no)/length(visibility_other)>=.1
  present <- Matrix::rowSums(d$raw_present[,pool,drop=FALSE])>0
  nmut <- Matrix::rowSums(d$mutated[,pool,drop=FALSE])
  indices <- unname(which(visible & present & nmut>=5))
  a <- Matrix::rowSums(d$mutated[indices,donor,drop=FALSE])
  c <- Matrix::rowSums(d$mutated[indices,comparison,drop=FALSE])
  n1 <- nd[indices]
  n2 <- rowSums(d$observed[,comparison,drop=FALSE])[d$pos[indices]]
  total <- n1+n2;pr <- (a+c)/total
  observed <- cbind(a,n1-a,c,n2-c)
  expected <- cbind(n1*pr,n1*(1-pr),n2*pr,n2*(1-pr))
  delta <- abs(observed-expected)
  correction <- pmin(.5,delta[,1],delta[,2],delta[,3],delta[,4])
  p <- pchisq(rowSums((delta-correction)^2/expected),df=1,lower.tail=FALSE)
  q <- p.adjust(p,'BH')
  keep <- !is.na(p) & p<=.05 & q<=.1 & a*n2>c*n1
  order <- order(q)
  features <- indices[order[keep[order]]]
  list(features=features,observable_variants=sum(visible & present),
       stats=data.frame(ID=d$ids[indices],index=indices,a=a,b=n1-a,c=c,d=n2-c,p=p,q=q))
}

background <- function(d,features,reference) {
  m <- Matrix::rowSums(d$alt[features,reference,drop=FALSE])
  n <- rowSums(d$coverage[d$pos[features],reference,drop=FALSE])
  p <- (1-m/n)*(.001/3)+(m/n)*(1-.001)
  sd <- (1-4*.001/3)*sqrt(m*(n-m)/n^3)
  list(m=m,n=n,lower=pmax(p-1.96*sd,0),upper=pmin(p+1.96*sd,1))
}

neff <- function(d,features,candidates,bg) {
  s <- summary(d$af[features,candidates,drop=FALSE])
  dp <- d$coverage[cbind(d$pos[features[s$i]],candidates[s$j])]
  keep <- dp>=5 & s$x>0
  s <- s[keep,];dp <- dp[keep];mi <- s$x*dp
  w <- (1-pbeta(bg$lower[s$i],mi+1,dp-mi+1))/pbeta(bg$upper[s$i],mi+1,dp-mi+1)
  counts <- tabulate(s$j[which(w>1)],nbins=length(candidates))
  observable <- colSums(d$observed[d$pos[features],candidates,drop=FALSE])
  list(count=counts,observed=observable,nan_weights=sum(is.na(w)),infinite_weights=sum(is.infinite(w)))
}

unique_rank <- function(x,decreasing=FALSE) {
  rank <- x;rank[order(x,decreasing=decreasing)] <- seq_along(x);rank
}
write_gz <- function(d,path) {
  con <- gzfile(path,'wt');on.exit(close(con));write.table(d,con,sep='\t',quote=FALSE,row.names=FALSE)
}
confusion <- function(z,truth) {
  tp <- sum(z & truth);fp <- sum(z & !truth);tn <- sum(!z & !truth);fn <- sum(!z & truth)
  data.frame(TP=tp,FP=fp,TN=tn,FN=fn,
    sensitivity=if(sum(truth))tp/sum(truth) else NA_real_,
    specificity=if(sum(!truth))tn/sum(!truth) else NA_real_,
    precision=if(tp+fp)tp/(tp+fp) else NA_real_,
    F1=if(2*tp+fp+fn)2*tp/(2*tp+fp+fn) else NA_real_,
    FPR=if(sum(!truth))fp/sum(!truth) else NA_real_)
}

check_enrichment <- function(s) {
  p <- or <- rep(NA_real_,nrow(s$stats))
  for(i in seq_len(nrow(s$stats))) {
    x <- as.numeric(s$stats[i,c('a','b','c','d')])
    m <- matrix(x,2,2)
    p[i] <- suppressWarnings(chisq.test(m)$p.value)
    or[i] <- fisher.test(m)$estimate
  }
  q <- p.adjust(p,'BH');ord <- order(q)
  chosen <- s$stats$index[ord[which(p[ord]<=.05 & q[ord]<=.1 & or[ord]>1)]]
  stopifnot(identical(chosen,s$features),identical(is.na(p),is.na(s$stats$p)))
  max_error <- max(abs(p-s$stats$p),na.rm=TRUE)
  stopifnot(max_error<1e-10)
  list(all_actual_tables=nrow(s$stats),features_identical=TRUE,max_chisq_p_error=max_error,
       reference='R chisq.test/fisher.test/p.adjust on every first-split table; same upstream expressions')
}

check_neff <- function(d,features,donor,candidates,reference,fast) {
  ids <- unique(round(seq(1,length(features),length.out=min(12,length(features)))))
  fi <- features[ids]
  ti <- unique(round(seq(1,length(candidates),length.out=12)))
  target <- candidates[ti];all <- unique(c(donor,reference,target))
  m <- as.matrix(d$af[fi,all,drop=FALSE]);depth <- d$coverage[d$pos[fi],all,drop=FALSE]
  m[depth<5] <- NA_real_
  cov <- d$coverage[unique(d$pos[fi]),all,drop=FALSE]
  actual <- Cell_Neff_cal(m,d$variants,cov,d$cells[donor],d$cells[target],
                         Ref_nonReceivers=d$cells[reference],mutFeatures=d$ids[fi],
                         adjust=TRUE,adjust.method='LR',e=.001)
  partial <- neff(d,fi,target,background(d,fi,reference))
  stopifnot(identical(as.integer(actual$DMTvar_count),partial$count))
  bg <- background(d,features,reference)
  minfo <- data.frame(vaf_po_receiver=bg$m/bg$n,s_reads_receiver=bg$m,
                      t_reads_receiver=bg$n,row.names=d$ids[features])
  small <- as.matrix(d$af[features,target,drop=FALSE])
  dp <- d$coverage[d$pos[features],target,drop=FALSE];small[dp<5] <- NA_real_
  original <- vapply(d$cells[target],function(cell)
      sum(w_cal2(d$ids[features],cell,minfo,small,d$coverage,e=.001)>1,na.rm=TRUE),0L)
  stopifnot(identical(unname(original),fast$count[ti]))
  list(cells_checked=length(target),full_feature_count=length(features),
       original_Cell_Neff_cal_feature_fixture=length(fi),
       original_w_cal2_all_features_equal=TRUE,full_background_held_fixed=TRUE,
       original_full_Neff_chain_rerun=FALSE)
}
