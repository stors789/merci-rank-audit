# Read the published wide integer coverage with at most one matrix copy.
read_coverage <- function(path, cells) {
  header <- strsplit(readLines(path, n=1L), '\t', fixed=TRUE)[[1]]
  table <- data.table::fread(path, skip=1L, header=FALSE, data.table=TRUE, nThread=1L)
  positions <- as.character(table[[1L]])
  if (length(header)==ncol(table)) header <- header[-1L]
  stopifnot(length(header)==ncol(table)-1L, !anyDuplicated(positions),
            !anyDuplicated(header), length(cells)==length(header),
            !anyDuplicated(cells), setequal(header,cells))
  data.table::setnames(table, c('position',header))
  data.table::setcolorder(table,c('position',cells))
  data.table::set(table,j='position',value=NULL)
  stopifnot(all(vapply(table,is.integer,logical(1))),
            all(vapply(table,function(x)all(x>=0L),logical(1))))
  coverage <- as.matrix(table)
  rownames(coverage) <- positions
  rm(table); invisible(gc())
  observed <- coverage>=5L
  list(coverage=coverage,observed=observed,
       mean_MT_coverage=colMeans(coverage),observed_MT_positions=colSums(observed))
}
