#!/usr/bin/env Rscript
# Adapt the accepted bench2 preparation and scoring to one explicit public-input case.
args <- commandArgs(trailingOnly=TRUE)
if (length(args)!=7L) stop('Expected info, expression, variants, coverage, MERCI root, engine and output')
names(args) <- c('info','expression','variants','coverage','upstream','engine','output')
userlib <- normalizePath(Sys.getenv('R_LIBS_USER'),mustWork=TRUE)
.Library.site <- character()
.libPaths(c(userlib,.Library))
required <- c('Matrix','data.table','jsonlite','e1071','proxy')
required_paths <- vapply(required,function(x) normalizePath(find.package(x,lib.loc=userlib),mustWork=TRUE),'')
if (any(dirname(required_paths)!=userlib)) stop('Required R packages must come from the explicit isolated library')
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(jsonlite))
setDTthreads(1)
write_json(list(R=R.version.string, user_library=userlib, library_paths=.libPaths(),
                required_packages=as.list(required_paths), data_table_threads=getDTthreads()),
           file.path(args['output'],'R_library_receipt.json'),pretty=TRUE,auto_unbox=TRUE)
source(args['engine'])
for (name in c('MERCI_MT_est.r','MERCI_ReceiverPre.r')) source(file.path(args['upstream'],'MERCI/R',name))
markers <- c('mt-Nd1','mt-Nd2','mt-Co1','mt-Co2','mt-Atp8','mt-Atp6','mt-Co3',
             'mt-Nd3','mt-Nd4l','mt-Nd4','mt-Nd5','mt-Nd6','mt-Cytb')
started <- proc.time()

info <- fread(args['info'],skip=1,header=FALSE,col.names=c('cell','sample','label','MT_reads'),data.table=FALSE)
stopifnot(nrow(info)==5251,!anyDuplicated(info$cell),all(info$MT_reads>1000))
donor <- info$cell[info$label=='MC_T_cell']
positive <- info$cell[info$label=='CC_cancer_cell']
negative <- info$cell[info$label=='MC_cancer_cell']
stopifnot(length(donor)==1557,length(positive)==1151,length(negative)==2543)
RNGkind('Mersenne-Twister','Inversion','Rejection')
set.seed(202610291L)
reference <- sort(sample(negative,1271L));test <- sort(setdiff(negative,reference))
stopifnot(length(reference)==1271L,length(test)==1272L,!length(intersect(reference,test)))
roles <- list(donor=donor,reference=reference,test=test,seed=202610291L,RNG=RNGkind())
saveRDS(roles,file.path(args['output'],'roles.rds'),compress=FALSE)
cat('roles_fixed donor',length(donor),'reference',length(reference),'test',length(test),'\n');flush.console()

read_wide <- function(path) {
  header <- strsplit(readLines(path,n=1),'\t',fixed=TRUE)[[1]]
  table <- fread(path,skip=1,header=FALSE,data.table=FALSE,nThread=1)
  ids <- as.character(table[[1]]);if(length(header)==ncol(table))header <- header[-1]
  stopifnot(length(header)==ncol(table)-1,!anyDuplicated(ids),!anyDuplicated(header))
  matrix <- as.matrix(table[-1]);dimnames(matrix) <- list(ids,header);matrix
}
coverage <- read_wide(args['coverage']);cells <- sort(info$cell)
coverage <- coverage[,cells,drop=FALSE]
stopifnot(is.integer(coverage),all(coverage>=0))
cat('coverage_loaded',nrow(coverage),ncol(coverage),'\n');flush.console()

con <- file(args['expression'],'rt')
header <- strsplit(readLines(con,n=1),'\t',fixed=TRUE)[[1]]
header <- gsub('.','-',header,fixed=TRUE)
stopifnot(length(header)==5251,!anyDuplicated(header),setequal(header,cells))
selected <- list()
repeat {
  chunk <- readLines(con,n=256);if(!length(chunk))break
  gene <- sub('\t.*','',chunk)
  for(j in which(gene %in% markers)) {
    values <- strsplit(chunk[j],'\t',fixed=TRUE)[[1]]
    stopifnot(length(values)==length(header)+1,!values[1] %in% names(selected))
    selected[[values[1]]] <- as.numeric(values[-1])
  }
}
close(con);stopifnot(setequal(names(selected),markers))
expr <- do.call(rbind,selected[markers]);colnames(expr) <- header;expr <- expr[,cells,drop=FALSE]
stopifnot(all(is.finite(expr)),all(colSums(expr)>0))
cat('expression_loaded',nrow(expr),ncol(expr),'\n');flush.console()

v <- fread(args['variants'],data.table=FALSE)
names(v) <- c('Cell','ID','pos','s_reads','t_reads','AF');v <- v[v$Cell %in% cells,]
ids <- sort(unique(v$ID));ri <- match(v$ID,ids);ci <- match(v$Cell,cells)
positions <- vapply(strsplit(ids,'_',fixed=TRUE),`[`,'',2)
pos <- match(positions,rownames(coverage));stopifnot(!anyNA(pos))
dims <- c(length(ids),length(cells));dn <- list(ids,cells)
af <- sparseMatrix(i=ri,j=ci,x=v$AF,dims=dims,dimnames=dn)
alt <- sparseMatrix(i=ri,j=ci,x=v$s_reads,dims=dims,dimnames=dn)
present <- sparseMatrix(i=ri,j=ci,x=1,dims=dims,dimnames=dn)
keep <- coverage[cbind(pos[ri],ci)]>=5 & v$AF>0
mutated <- sparseMatrix(i=ri[keep],j=ci[keep],x=1,dims=dims,dimnames=dn)
observed <- coverage>=5
technical <- data.frame(cell=cells,MT_reads=info$MT_reads[match(cells,info$cell)],
  mean_MT_coverage=colMeans(coverage),observed_MT_positions=colSums(observed),nuclear_UMI=NA_real_)
d <- list(cells=cells,ids=ids,pos=pos,coverage=coverage,observed=observed,af=af,alt=alt,
          raw_present=present,mutated=mutated,expression=expr,variants=v,technical=technical,info=info)
donor_index <- cols(d,donor);reference_index <- cols(d,reference);candidate_index <- cols(d,test)
visibility <- unique(c(candidate_index,reference_index))
pool <- unique(c(donor_index,visibility))
feature <- select_features(d,donor_index,reference_index,visibility,pool)
if(!length(feature$features))stop('No source-enriched features for the fixed case')
backgrounds <- background(d,feature$features,reference_index)
dna <- neff(d,feature$features,candidate_index,backgrounds)
cat('features_selected',length(feature$features),'\n');flush.console()

pdf(NULL)
rna <- MERCI_MT_est(d$expression,test,donor,reference,organism='mouse')
stopifnot(identical(rownames(rna),test),all(is.finite(rna$Donor_MT_ind)))
rna_rank <- unique_rank(rna$Donor_MT_ind)
stopifnot(identical(unname(rna_rank),unname(rna$W_Donor_rank)))
dna_rank <- unique_rank(dna$count)
dn <- data.frame(MTvar_rank=dna_rank,row.names=test)
rna_input <- data.frame(Donor_MT_ind=rna$Donor_MT_ind,W_Donor_rank=rna_rank,row.names=test)
cutoffs <- seq(100L,10L,-10L)
prediction <- sapply(cutoffs,function(cut)MERCI_ReceiverPre(dn,rna_input,top_rank=cut)$prediction=='Receiver')
counts <- colSums(prediction)
dev.off()

fast_random <- function(n,cutoff,b) {
  k <- ceiling(quantile(seq_len(n),cutoff/100))
  vapply(seq_len(b),function(i)sum(sample.int(n)>n-k & sample.int(n)>n-k),0L)
}
set.seed(2026102901L)
random <- sapply(cutoffs,function(cut)fast_random(length(test),cut,10000L))
maxima <- apply(random,2,max)
saveRDS(random,file.path(args['output'],'random_nulls.rds'),compress=FALSE)
write.csv(data.frame(N=length(test),cutoff=cutoffs,maximum=maxima),
          file.path(args['output'],'thresholds.csv'),row.names=FALSE)
tech <- d$technical[candidate_index,]
score <- data.frame(cell=test,truth=FALSE,Neff=dna$count,DNA_rank=dna_rank,
  RNA_score=rna$Donor_MT_ind,RNA_rank=rna_rank,observed_source_features=dna$observed,
  source_marker_observable=dna$observed>0,MT_reads=tech$MT_reads,mean_MT_coverage=tech$mean_MT_coverage,
  observed_MT_positions=tech$observed_MT_positions,nuclear_UMI=NA_real_,raw_top50=prediction[,cutoffs==50L])
write.table(score,file.path(args['output'],'score.tsv'),sep='\t',quote=FALSE,row.names=FALSE,na='NA')
saveRDS(list(feature=feature,feature_ids=d$ids[feature$features],all_variant_ids=d$ids),
        file.path(args['output'],'feature_selection.rds'),compress=FALSE)
receipt <- list(status='generated_from_four_explicit_processed_inputs',case='bench2_split1_original_regular_negative',
  cells=length(cells),donor=length(donor),reference=length(reference),target=length(test),
  markers=markers,variants=length(ids),observable_variants=feature$observable_variants,
  tested_variants=nrow(feature$stats),selected_features=length(feature$features),
  raw_top50=sum(score$raw_top50),gate=any(counts>maxima),cutoffs=cutoffs,observed=as.integer(counts),
  maxima=as.integer(maxima),random_draws=10000L,split_seed=202610291L,null_seed=2026102901L,
  RNG=RNGkind(),RNA='original MERCI_MT_est on fixed references; negative case only',
  seconds=unname((proc.time()-started)['elapsed']))
write_json(receipt,file.path(args['output'],'generation_receipt.json'),pretty=TRUE,auto_unbox=TRUE)
cat(toJSON(receipt,auto_unbox=TRUE),'\n')
