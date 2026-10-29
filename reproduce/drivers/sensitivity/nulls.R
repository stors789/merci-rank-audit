# Original BIO-04/BIO-06 finite null mechanism; explicit output, no repository paths.
args <- commandArgs(trailingOnly=TRUE)
if (!(length(args) %in% c(3L,4L))) stop('Expected N, frozen R seed, new output directory, and optional draw count')
n <- as.integer(args[1]); seed <- as.integer(args[2]); out <- args[3]
draws <- if(length(args)==4L) as.integer(args[4]) else 10000L
if (is.na(n) || n<1L || is.na(seed) || is.na(draws) || draws<1L || dir.exists(out)) stop('Invalid input or existing output')
dir.create(out, recursive=TRUE)
RNGkind('Mersenne-Twister','Inversion','Rejection'); set.seed(seed)
cutoffs <- seq(100L,10L,-10L)
fast_random <- function(n,cutoff,b) {
  k <- ceiling(quantile(seq_len(n),cutoff/100))
  vapply(seq_len(b),function(i)sum(sample.int(n)>n-k & sample.int(n)>n-k),0L)
}
values <- vapply(cutoffs,function(cut)fast_random(n,cut,draws),integer(draws))
colnames(values) <- paste0('cutoff_',cutoffs)
con <- gzfile(file.path(out,'draws.tsv.gz'),'wt')
write.table(values,con,sep='\t',row.names=FALSE,col.names=TRUE,quote=FALSE); close(con)
summary <- data.frame(N=n,cutoff=cutoffs,K=vapply(cutoffs,function(cut)as.integer(ceiling(quantile(seq_len(n),cut/100))),0L),
  draws=draws,R_seed=seed,maximum=apply(values,2,max),minimum=apply(values,2,min))
write.csv(summary,file.path(out,'maxima.csv'),row.names=FALSE,quote=FALSE)
writeLines(c(paste('R_version',R.version.string),paste('RNGkind',paste(RNGkind(),collapse=','))),file.path(out,'metadata.txt'))
