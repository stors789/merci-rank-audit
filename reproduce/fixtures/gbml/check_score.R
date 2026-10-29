#!/usr/bin/env Rscript
# Check the single-feature name repair against the original on synthetic inputs.
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=2L)stop('Expected upstream R-library')
userlib <- normalizePath(args[2],mustWork=TRUE);.Library.site <- character();.libPaths(c(userlib,.Library))
up <- file.path(normalizePath(args[1],mustWork=TRUE),'MERCI/R')
original <- new.env(parent=globalenv());adapted <- new.env(parent=globalenv())
for (f in c('Cell_Neff_cal.r','w_cal2.r','mut_inf_extr.r')) {
  sys.source(file.path(up,f),envir=original);sys.source(file.path(up,f),envir=adapted)
}
lines <- readLines(file.path(up,'mut_inf_extr.r'))
before <- 't.results <- data.frame(t(t.results))'
after <- 't.results <- data.frame(t(t.results)); rownames(t.results) <- muts'
stopifnot(sum(grepl(before,lines,fixed=TRUE))==1L)
lines <- gsub(before,after,lines,fixed=TRUE)
con <- textConnection(lines);source(con,local=adapted);close(con)
d <- paste0('d',1:10);r <- paste0('r',1:10);cells <- c(d,r)
features <- c('MT_1_A-C','MT_2_C-G')
m <- matrix(0,2,20,dimnames=list(features,cells))
m[1,d] <- .5;m[1,r[1:5]] <- .2;m[2,d[1:5]] <- .3;m[2,r[1]] <- .1
coverage <- matrix(10L,2,20,dimnames=list(c('1','2'),cells))
at <- which(m>0,arr.ind=TRUE)
variants <- data.frame(Cell=cells[at[,2]],ID=features[at[,1]],pos=at[,1],
                       s_reads=round(10*m[at]),t_reads=10,AF=m[at],Cell_reads=20)
pdf(tempfile(fileext='.pdf'))
for (subset in list(features[1],features)) {
  args <- list(varMatrix=m,MT_variants=variants,MT_coverage=coverage,donor_cells=d,mixed_cells=r,
               mutFeatures=subset,adjust=TRUE,adjust.method='LR',frac_tran=FALSE,e=.001)
  a <- suppressWarnings(do.call(original$Cell_Neff_cal,args))
  b <- suppressWarnings(do.call(adapted$Cell_Neff_cal,args))
  stopifnot(identical(a,b))
}
dev.off()
cat('PASS: single- and multi-feature Cell_Neff_cal match original exactly\n')
