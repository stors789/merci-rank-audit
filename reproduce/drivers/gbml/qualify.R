#!/usr/bin/env Rscript
# The original MERCI transformation and source statistic, with recorded shape and exact-5% repairs.
args <- commandArgs(trailingOnly=TRUE)
if(!length(args) %in% c(4L,5L)) stop('Expected upstream input output R-library [--self-test]')
self_test <- length(args)==5L && args[5]=='--self-test'
up <- normalizePath(args[1],mustWork=TRUE);input <- args[2];out <- args[3]
userlib <- normalizePath(args[4],mustWork=TRUE)
.Library.site <- character();.libPaths(c(userlib,.Library))
suppressPackageStartupMessages(library(jsonlite))
if(!self_test) { if(dir.exists(out))stop('Protect output');dir.create(out,recursive=TRUE) }
stopifnot(identical(system2('git', c('-C', up, 'rev-parse', 'HEAD'), stdout=TRUE),
                    '79e105aaf13159bee325c3ceceb6ee679e430341'))
e <- new.env(parent=globalenv())
for (f in c('readMTvar_10x.r', 'MTmutMatrix_refined.r', 'MTmutMatrix_transform.r', 'Denrich_mtMut_extr.r'))
  sys.source(file.path(up, 'MERCI/R', f), envir=e)

load_filter <- function() {
  lines <- readLines(file.path(up, 'MERCI/R/filter_MTmuts.r'))
  first <- grep('t_index <- which\\(apply\\(observeRate_stat', lines)
  stopifnot(length(first)==1L, grepl('all( x >= min_observeRate )', lines[first+1L], fixed=TRUE),
            grepl('}) )', lines[first+2L], fixed=TRUE))
  replacement <- c(
    '\tt_index <- which(apply(filter_matrix2, 1, function(x) {',
    '\t\tall(20 * sum(!is.na(x)) >= cellNum,',
    '\t\t    20 * sum(!is.na(x[seq_len(doN)])) >= doN,',
    '\t\t    20 * sum(!is.na(x[doN + seq_len(reN)])) >= reN)',
    '\t}) )')
  lines <- append(lines[-(first:(first+2L))], replacement, after=first-1L)
  shape_before <- c('filter_matrix2 <- filter_matrix[, c(donor_cells, reciever_cells)]',
                    'New_mutMatrix <- filter_matrix[selected_variants, ]')
  shape_after <- c('filter_matrix2 <- filter_matrix[, c(donor_cells, reciever_cells), drop=FALSE]',
                   'New_mutMatrix <- filter_matrix[selected_variants, , drop=FALSE]')
  for (i in seq_along(shape_before)) {
    stopifnot(sum(grepl(shape_before[i], lines, fixed=TRUE))==1L)
    lines <- gsub(shape_before[i], shape_after[i], lines, fixed=TRUE)
  }
  con <- textConnection(lines); source(con, local=e); close(con)
}
load_filter()
lines <- readLines(file.path(up, 'MERCI/R/Donor_enriched_vars.r'))
before <- c('t_varMatrix[t_index, ]', 't_varMatrix[, Cluster_cells]', 't_varMatrix[, Other_cells]')
after <- c('t_varMatrix[t_index, , drop=FALSE]', 't_varMatrix[, Cluster_cells, drop=FALSE]',
           't_varMatrix[, Other_cells, drop=FALSE]')
for (i in seq_along(before)) {
  stopifnot(sum(grepl(before[i], lines, fixed=TRUE))==1L)
  lines <- gsub(before[i], after[i], lines, fixed=TRUE)
}
con <- textConnection(lines); source(con, local=e); close(con)

if (self_test) {
  cells <- paste0('c', seq_len(40)); donor <- cells[1:20]; target <- cells[21:40]
  m <- matrix(NA_real_, 3, 40, dimnames=list(c('MT_10_A-C','MT_20_A-G','MT_30_C-T'),cells))
  cov <- matrix(0L, 3, 40, dimnames=list(c('10','20','30'),cells))
  cov[1,c(1,21)] <- 2L; m[1,c(1,21)] <- c(.5,0)
  cov[2,1] <- 2L; m[2,1] <- .5
  cov[3,] <- 2L; m[3,] <- 0; m[3,1:5] <- .5
  result <- e$filter_MTmuts(m,cov,donor,target,min_d=2,min_observeRate=.05)
  stopifnot(identical(rownames(result),c('MT_10_A-C','MT_30_C-T')),
            is.matrix(result), is.na(result[1,2]), result[1,21]==0)
  one <- e$filter_MTmuts(m[3,,drop=FALSE],cov,donor,target,min_d=2,min_observeRate=.05)
  stopifnot(is.matrix(one),nrow(one)==1L)
  d <- paste0('d',1:20); r <- paste0('r',1:40)
  fixture <- matrix(0,3,60,dimnames=list(c('MT_1_A-C','MT_2_C-G','MT_3_G-T'),c(d,r)))
  fixture[1,d[1:15]] <- .5; fixture[1,r[1:5]] <- .5
  fixture[2,d[1:5]] <- .5; fixture[2,r[1:30]] <- .5
  fixture[3,c(d[1:10],r[1:20])] <- .5
  original <- new.env(parent=globalenv())
  sys.source(file.path(up,'MERCI/R/Donor_enriched_vars.r'),envir=original)
  sys.source(file.path(up,'MERCI/R/Denrich_mtMut_extr.r'),envir=original)
  expected <- suppressWarnings(original$Denrich_mtMut_extr(fixture,d,r))
  actual <- suppressWarnings(e$Denrich_mtMut_extr(fixture,d,r))
  stopifnot(identical(expected,actual),identical(actual,'MT_1_A-C'))
  full <- rbind(fixture, 'MT_4_T-A'=rep(.5,60))
  coverage <- matrix(2L,4,60,dimnames=list(as.character(1:4),c(d,r)))
  at <- which(full>0,arr.ind=TRUE)
  records <- data.frame(Cell=colnames(full)[at[,2]],ID=rownames(full)[at[,1]],
                        pos=at[,1],s_reads=1,t_reads=2,AF=full[at],Cell_reads=11)
  records <- rbind(records,data.frame(Cell='low',ID='MT_4_T-A',pos=4,s_reads=1,t_reads=2,AF=.5,Cell_reads=10))
  path <- tempfile();write.table(records,path,sep='\t',quote=FALSE,row.names=FALSE)
  loaded <- e$readMTvar_10x(path,minReads=10);unlink(path)
  transformed <- e$MTmutMatrix_transform(loaded,coverage,d,r,min_d=2,min_observeRate=.05)
  stopifnot(!('low'%in%loaded$Cell),setequal(colnames(transformed),colnames(full)),
            isTRUE(all.equal(unname(transformed[rownames(full),colnames(full)]),unname(full))))
  cat('PASS: strict >10, original transform, exact 5% per group, covered-zero, single-row shape, source defaults\n')
  quit(status=0L)
}

if (file.exists(file.path(out, 'source_qualification.json')) || file.exists(file.path(input, 'qualification.rds')))
  stop('Protect first source qualification')
started <- proc.time()
warnings_seen <- list()
result <- withCallingHandlers(tryCatch({
  donors <- readLines(file.path(input,'donors.txt')); targets <- readLines(file.path(input,'targets.txt'))
  variants <- e$readMTvar_10x(file.path(input,'MT_variants.tsv'), minReads=10)
  present <- unique(variants$Cell)
  donors_used <- intersect(donors,present); targets_used <- intersect(targets,present)
  common <- list(frozen_donors=length(donors), frozen_targets=length(targets),
                 variant_records_after_minReads=nrow(variants), donors_with_variant_record=length(donors_used),
                 targets_with_variant_record=length(targets_used))
  if (!nrow(variants) || !length(donors_used) || !length(targets_used)) {
    c(common, list(status='not_evaluable_no_variant_records_in_both_roles', source_features=0L,
                   downstream='not_run'))
  } else {
    coverage <- as.matrix(read.delim(gzfile(file.path(input,'Coverage.tsv.gz')), row.names=1,
                                    check.names=FALSE))
    stopifnot(is.numeric(coverage),all(coverage>=0),all(coverage==floor(coverage)),
              !anyDuplicated(variants[c('Cell','ID')]),all(present %in% colnames(coverage)))
    m <- e$MTmutMatrix_transform(variants,coverage,donors,targets,min_d=2,min_observeRate=.05)
    stopifnot(is.matrix(m),setequal(colnames(m),present))
    mutated <- if (nrow(m)) rowSums(m>0,na.rm=TRUE) else integer(0)
    stats <- NULL; features <- character(0)
    if (any(mutated>=5)) {
      holder <- new.env(parent=emptyenv()); e$.captured_statistic <- holder
      trace('Donor_enriched_vars', exit=quote(assign('statistics',t_stat,envir=.captured_statistic)),
            where=e, print=FALSE)
      features <- e$Denrich_mtMut_extr(m,donors_used,targets_used)
      untrace('Donor_enriched_vars',where=e)
      stats <- holder$statistics; stats$ID <- rownames(stats)
      stopifnot(setequal(features,stats$ID[stats$catogary=='donor_enrich']))
      write.csv(stats,file.path(out,'source_feature_statistics.csv'),row.names=FALSE,na='NA')
    }
    saveRDS(list(matrix=m,variants=variants,coverage=coverage,donors=donors_used,
                 targets=targets_used,features=features),file.path(input,'qualification.rds'),compress=FALSE)
    c(common, list(status=if(length(features)) 'source_features_identified' else 'no_qualified_source_features',
                   variants_before_observability=length(unique(variants$ID)),
                   variants_after_observability=nrow(m), variants_with_mutated_cells_ge5=sum(mutated>=5),
                   source_features=length(features), source_feature_ids=features,
                   downstream=if(length(features)) 'eligible_for_score_preparation' else 'not_run'))
  }
},error=function(err) list(status='run_failed', error=conditionMessage(err), downstream='not_run')),
warning=function(w) {
  msg <- conditionMessage(w)
  warnings_seen[[msg]] <<- if (is.null(warnings_seen[[msg]])) 1L else warnings_seen[[msg]]+1L
  invokeRestart('muffleWarning')
})
report <- list(result=result, warnings=warnings_seen, minReads_strict_gt=10, min_d=2, min_observeRate=.05,
               group_observability='20 * observed >= group size for donor, target, total',
               source_defaults=list(OR=2,Nmut_min=5,pvalue=1,qvalue=1),
               shape_repairs=before, seconds=unname((proc.time()-started)['elapsed']))
write_json(report,file.path(out,'source_qualification.json'),pretty=TRUE,auto_unbox=TRUE,na='null')
cat(toJSON(report,auto_unbox=TRUE,na='null'), '\n')
if (identical(result$status,'run_failed')) quit(status=1L)
