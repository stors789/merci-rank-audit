args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=2L)stop('Expected coverage helper and synthetic TSV')
suppressPackageStartupMessages(library(data.table));setDTthreads(1L);source(args[1L])
x <- read_coverage(args[2L],c('cell_a','cell_b','cell_c'))
expected <- matrix(c(5L,4L,8L,10L,1L,5L,0L,6L,2L),nrow=3L,
                   dimnames=list(c('p1','p2','p3'),c('cell_a','cell_b','cell_c')))
stopifnot(identical(x$coverage,expected),identical(x$observed,expected>=5L),
          isTRUE(all.equal(unname(x$mean_MT_coverage),unname(colMeans(expected)))),
          identical(unname(x$observed_MT_positions),unname(colSums(expected>=5L))))
cat('PASS: explicit axes, integer coverage and original technical denominators\n')
