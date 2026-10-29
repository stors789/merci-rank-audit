#!/usr/bin/env Rscript
# B1 frozen strict chain from four explicit author-processed files.
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=11L) stop('Expected info expression variants coverage upstream engine coverage-helper output split original-maxima R-library')
names(args) <- c('info','expression','variants','coverage','upstream','engine','coverage_helper','output','split','maxima','r_library')
split <- as.integer(args['split'])
if(is.na(split) || split<21L || split>100L) stop('B1 frozen split must be 21..100')
userlib <- normalizePath(args['r_library'],mustWork=TRUE)
.Library.site <- character(); .libPaths(c(userlib,.Library))
required <- c('Matrix','data.table','jsonlite','e1071','proxy')
required_paths <- vapply(required,function(x)normalizePath(find.package(x,lib.loc=userlib),mustWork=TRUE),'')
if(any(dirname(required_paths)!=userlib)) stop('R dependencies must come from the explicit user library')
suppressPackageStartupMessages(library(Matrix));suppressPackageStartupMessages(library(data.table));suppressPackageStartupMessages(library(jsonlite))
setDTthreads(1L)
source(args['engine']);source(args['coverage_helper'])
for(name in c('MERCI_MT_est.r','MERCI_ReceiverPre.r'))source(file.path(args['upstream'],'MERCI','R',name))
started <- proc.time()
info <- fread(args['info'],skip=1L,header=FALSE,col.names=c('cell','sample','label','MT_reads'),data.table=FALSE,nThread=1L)
stopifnot(nrow(info)==12145L,!anyDuplicated(info$cell),all(info$MT_reads>1000))
donor <- info$cell[info$label=='MC_T_cell'];positive <- info$cell[info$label=='CC_cancer_cell'];negative <- info$cell[info$label=='MC_cancer_cell']
stopifnot(length(donor)==1669L,length(positive)==2904L,length(negative)==7572L)
RNGkind('Mersenne-Twister','Inversion','Rejection');split_seed <- 202610260L+split;set.seed(split_seed)
reference <- sort(sample(negative,3786L));test <- sort(setdiff(negative,reference));mixed_names <- sort(c(positive,test))
stopifnot(length(reference)==3786L,length(test)==3786L,length(mixed_names)==6690L,!length(intersect(reference,test)))
saveRDS(list(donor=donor,positive=positive,negative=negative,reference=reference,test=test,split=split,seed=split_seed,RNG=RNGkind()),
        file.path(args['output'],'roles.rds'),compress=FALSE)
cat('ROLES_FIXED split',split,'reference',length(reference),'test',length(test),'CC',length(positive),'\n');flush.console()
cells <- sort(info$cell)
# Preserve the original ^mt- feature order from the public normalized expression file.
con <- file(args['expression'],'rt');header <- strsplit(readLines(con,n=1L),'\t',fixed=TRUE)[[1]]
stopifnot(length(header)==length(cells),!anyDuplicated(header),setequal(header,cells))
selected <- list();gene_rows <- 0L
repeat {
 chunk <- readLines(con,n=128L);if(!length(chunk))break
 gene_rows <- gene_rows+length(chunk);genes <- sub('\t.*','',chunk)
 for(j in which(startsWith(genes,'mt-'))) {
  values <- strsplit(chunk[j],'\t',fixed=TRUE)[[1]]
  stopifnot(length(values)==length(header)+1L,!values[1L] %in% names(selected))
  selected[[values[1L]]] <- as.numeric(values[-1L])
 }
}
close(con);stopifnot(length(selected)==13L)
expr <- do.call(rbind,selected);colnames(expr) <- header;expr <- expr[,cells,drop=FALSE]
stopifnot(all(is.finite(expr)))
cat('EXPRESSION_READY genes',nrow(expr),'rows_scanned',gene_rows,'\n');flush.console()
v <- fread(args['variants'],data.table=FALSE,nThread=1L)
names(v) <- c('Cell','ID','pos','s_reads','t_reads','AF');v <- v[v$Cell %in% cells,]
stopifnot(!anyDuplicated(v[c('ID','Cell')]),all(is.finite(v$AF)),all(v$AF>=0 & v$AF<=1))
variant_ids <- sort(unique(v$ID));ri <- match(v$ID,variant_ids);ci <- match(v$Cell,cells)
cat('VARIANTS_READY records',nrow(v),'variants',length(variant_ids),'\n');flush.console()
coverage_data <- read_coverage(args['coverage'],cells);coverage <- coverage_data$coverage;observed <- coverage_data$observed
positions <- vapply(strsplit(variant_ids,'_',fixed=TRUE),`[`,'',2L)
pos <- match(positions,rownames(coverage));stopifnot(!anyNA(pos))
dims <- c(length(variant_ids),length(cells));dn <- list(variant_ids,cells)
af <- sparseMatrix(i=ri,j=ci,x=v$AF,dims=dims,dimnames=dn)
alt <- sparseMatrix(i=ri,j=ci,x=v$s_reads,dims=dims,dimnames=dn)
present <- sparseMatrix(i=ri,j=ci,x=1,dims=dims,dimnames=dn)
keep <- coverage[cbind(pos[ri],ci)]>=5L & v$AF>0
mutated <- sparseMatrix(i=ri[keep],j=ci[keep],x=1,dims=dims,dimnames=dn)
technical <- data.frame(cell=cells,MT_reads=info$MT_reads[match(cells,info$cell)],
  mean_MT_coverage=coverage_data$mean_MT_coverage,observed_MT_positions=coverage_data$observed_MT_positions,nuclear_UMI=NA_real_)
d <- list(cells=cells,ids=variant_ids,pos=pos,coverage=coverage,observed=observed,af=af,alt=alt,
          raw_present=present,mutated=mutated,expression=expr,variants=v,technical=technical,info=info)
rm(coverage_data,coverage,observed,af,alt,present,mutated,ri,ci,keep);invisible(gc())
cat('INPUTS_READY coverage',nrow(d$coverage),'x',ncol(d$coverage),'\n');flush.console()
donor_index <- cols(d,donor);reference_index <- cols(d,reference)
features <- select_features(d,donor_index,reference_index,reference_index,unique(c(donor_index,reference_index)))
if(!length(features$features))stop('No informative strict-reference source features')
saveRDS(features,file.path(args['output'],'feature_selection.rds'))
bg <- background(d,features$features,reference_index)
pdf(NULL)
regular_rna <- MERCI_MT_est(d$expression,mixed_names,donor,reference,organism='mouse')
stopifnot(identical(rownames(regular_rna),mixed_names),all(is.finite(regular_rna$Donor_MT_ind)))
cutoffs <- seq(100L,10L,-10L);old_maxima <- read.csv(args['maxima'],stringsAsFactors=FALSE)
case_receipts <- list()
for(kind in c('mixed','negative')) {
 candidate_names <- if(kind=='mixed')mixed_names else test;candidates <- cols(d,candidate_names)
 stopifnot(!length(intersect(unique(c(donor_index,reference_index)),candidates)))
 dna <- neff(d,features$features,candidates,bg)
 values <- regular_rna$Donor_MT_ind[match(candidate_names,rownames(regular_rna))]
 rna <- data.frame(Donor_MT_ind=values,W_Donor_rank=unique_rank(values),row.names=candidate_names)
 rank <- unique_rank(dna$count);dn <- data.frame(MTvar_rank=rank,row.names=candidate_names)
 predictions <- sapply(cutoffs,function(cut)MERCI_ReceiverPre(dn,rna,top_rank=cut)$prediction=='Receiver')
 counts <- colSums(predictions)
 maxima <- old_maxima[old_maxima$case==paste0('strict_regular_',kind),]
 maxima <- maxima[match(cutoffs,maxima$cutoff),]
 stopifnot(nrow(maxima)==10L,all(maxima$N==length(candidate_names)),!anyNA(maxima$maximum))
 gate <- any(counts>maxima$maximum);truth <- candidate_names %in% positive;tech <- d$technical[candidates,]
 score <- data.frame(cell=candidate_names,truth=truth,Neff=dna$count,DNA_rank=rank,RNA_score=rna$Donor_MT_ind,
   RNA_rank=rna$W_Donor_rank,observed_source_features=dna$observed,source_marker_observable=dna$observed>0,
   MT_reads=tech$MT_reads,mean_MT_coverage=tech$mean_MT_coverage,observed_MT_positions=tech$observed_MT_positions,
   nuclear_UMI=NA_real_,raw_top50=predictions[,cutoffs==50L])
 write_gz(score,file.path(args['output'],paste0('strict_regular_',kind,'_cells.tsv.gz')))
 case_receipts[[kind]] <- list(N=length(candidate_names),positive=sum(truth),gate=gate,counts=unname(as.integer(counts)),
    maxima=unname(as.integer(maxima$maximum)),raw_top50=sum(score$raw_top50),observed_features=sum(dna$observed>0L),
    NaN_weights=dna$nan_weights,infinite_weights=dna$infinite_weights)
 cat('SCORE_DONE',kind,'N',length(candidate_names),'raw',sum(score$raw_top50),'gate',gate,'\n');flush.console()
}
dev.off()
receipt <- list(status='SCORES_GENERATED_FROM_FOUR_PROCESSED_INPUTS',split=split,split_seed=split_seed,RNG=RNGkind(),
   counts=list(donor=length(donor),reference=length(reference),test=length(test),positive=length(positive),
               variants=length(variant_ids),features=length(features$features),coverage_positions=nrow(d$coverage)),
   RNA_genes=rownames(expr),source_selection='strict donor+reference only; no target visibility',
   cases=case_receipts,R_version=R.version.string,R_packages=as.list(vapply(required,function(x)as.character(packageVersion(x)),'')),
   R_library=userlib,elapsed_seconds=unname((proc.time()-started)['elapsed']))
write_json(receipt,file.path(args['output'],'generation_receipt.json'),pretty=TRUE,auto_unbox=TRUE,na='null')
cat(toJSON(receipt,auto_unbox=TRUE),'\n')
