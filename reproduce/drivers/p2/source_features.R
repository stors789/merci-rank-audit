# Retain the pinned statistic, correcting only three matrix drop dimensions.
args <- commandArgs(trailingOnly=TRUE)
if(!length(args) %in% c(4L,5L)) stop('Expected input output upstream R-library [--self-test]')
input <- args[1];output <- args[2];upstream_path <- normalizePath(args[3],mustWork=TRUE)
userlib <- normalizePath(args[4],mustWork=TRUE);.Library.site <- character();.libPaths(c(userlib,.Library))
self_test <- length(args)==5L && args[5]=='--self-test'
stopifnot(identical(system2('git',c('-C',upstream_path,'rev-parse','HEAD'),stdout=TRUE),'79e105aaf13159bee325c3ceceb6ee679e430341'))
original <- readLines(file.path(upstream_path,'MERCI/R/Donor_enriched_vars.r'))
before <- c('t_varMatrix[t_index, ]','t_varMatrix[, Cluster_cells]','t_varMatrix[, Other_cells]')
after <- c('t_varMatrix[t_index, , drop=FALSE]','t_varMatrix[, Cluster_cells, drop=FALSE]','t_varMatrix[, Other_cells, drop=FALSE]')
patched <- original
for(i in seq_along(before)) {
  stopifnot(sum(grepl(before[i],patched,fixed=TRUE))==1L)
  patched <- gsub(before[i],after[i],patched,fixed=TRUE)
}
load_stat <- function(lines) {
  e <- new.env(parent=globalenv())
  connection <- textConnection(lines)
  source(connection,local=e);close(connection)
  e
}
capture <- function(e,m,d,r) {
  holder <- new.env(parent=emptyenv());e$.captured_statistic <- holder
  trace('Donor_enriched_vars',exit=quote(assign('statistics',t_stat,envir=.captured_statistic)),where=e,print=FALSE)
  on.exit(untrace('Donor_enriched_vars',where=e))
  warnings_seen <- list()
  ids <- withCallingHandlers(e$Donor_enriched_vars(m,d,r,Nmut_min=5,pvalue=.05,qvalue=.1,OR=1),warning=function(w) {
    msg <- conditionMessage(w)
    warnings_seen[[msg]] <<- if(is.null(warnings_seen[[msg]])) 1L else warnings_seen[[msg]]+1L
    invokeRestart('muffleWarning')
  })
  list(ids=ids,statistics=holder$statistics,warnings=warnings_seen)
}
if(self_test) {
  d <- paste0('D',1:20);r <- paste0('R',1:40)
  m <- matrix(0,3,60,dimnames=list(c('MT_100_A-C','MT_200_C-G','MT_300_G-T'),c(d,r)))
  m[1,d] <- .5;m[2,r] <- .5;m[3,c(d[1:10],r[1:20])] <- .5
  a <- capture(load_stat(original),m,d,r);b <- capture(load_stat(patched),m,d,r)
  stopifnot(identical(a,b))
  one <- capture(load_stat(patched),m[1,,drop=FALSE],d,r)
  stopifnot(identical(one$ids,rownames(m)[1]),nrow(one$statistics)==1L)
  single_cell <- capture(load_stat(patched),m[,c(d[1],r),drop=FALSE],d[1],r)
  stopifnot(is.data.frame(single_cell$statistics))
  cat('PASS: original multi-variant statistic, single-variant and single-donor shape\n')
} else {
  m <- as.matrix(read.delim(file.path(input,'VAF.tsv'),row.names=1,check.names=FALSE,na.strings='NA'))
  d <- readLines(file.path(input,'donors.txt'));r <- readLines(file.path(input,'targets.txt'))
  result <- capture(load_stat(patched),m,d,r);st <- result$statistics;st$ID <- rownames(st)
  stopifnot(setequal(result$ids,st$ID[st$catogary=='donor_enrich']))
  write.csv(st,output,row.names=FALSE,na='NA')
  jsonlite::write_json(list(actual_warnings=result$warnings,dimension_only_subsets_preserved=TRUE,equations_and_thresholds_unchanged=TRUE),paste0(output,'.warnings.json'),pretty=TRUE,auto_unbox=TRUE)
  writeLines(capture.output(sessionInfo()),file.path(input,'R_session.txt'))
}
