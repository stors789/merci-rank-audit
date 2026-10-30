# Pinned Neff / human LOO / Rcm on an already qualified full-library source.
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=5L)stop('Expected input output upstream R-library threshold-table')
input <- normalizePath(args[1],mustWork=TRUE);out <- args[2]
upstream_path <- normalizePath(args[3],mustWork=TRUE);userlib <- normalizePath(args[4],mustWork=TRUE)
threshold_table <- normalizePath(args[5],mustWork=TRUE)
.Library.site <- character();.libPaths(c(userlib,.Library));suppressPackageStartupMessages(library(jsonlite))
stopifnot(identical(system2('git',c('-C',upstream_path,'rev-parse','HEAD'),stdout=TRUE),'79e105aaf13159bee325c3ceceb6ee679e430341'))
upstream <- new.env(parent=globalenv())
for(f in c('Cell_Neff_cal.r','w_cal2.r','mut_inf_extr.r','MERCI_LOO_MT_est.r','MERCI_ReceiverPre.r'))sys.source(file.path(upstream_path,'MERCI/R',f),envir=upstream)
stopifnot(!file.exists(out));dir.create(out,recursive=TRUE)
ready <- fromJSON(file.path(input,'ready.json'));stopifnot(identical(ready$status,'ready'))
m <- as.matrix(read.delim(file.path(input,'VAF.tsv'),row.names=1,check.names=FALSE,na.strings='NA'))
coverage <- as.matrix(read.delim(file.path(input,'Coverage.tsv'),row.names=1,check.names=FALSE))
v <- read.delim(file.path(input,'MT_variants.tsv'),check.names=FALSE)
expression <- as.matrix(read.delim(file.path(input,'MT_expression.tsv'),row.names=1,check.names=FALSE))
donors <- readLines(file.path(input,'donors.txt'));targets <- readLines(file.path(input,'targets.txt'));features <- readLines(file.path(input,'features.txt'))
stopifnot(length(donors)>=2L,length(targets)>=3L,length(features)>=1L,all(features%in%rownames(m)),all(c(donors,targets)%in%colnames(expression)))
stopifnot(!anyDuplicated(v[c('Cell','ID')]),all(v$s_reads>0),all(v$t_reads>=v$s_reads),all(abs(v$AF-v$s_reads/v$t_reads)<1e-12))
warnings_seen <- list();times <- list()
timed <- function(name,value) {
  started <- proc.time();cat(name,'\n');flush.console()
  answer <- withCallingHandlers(force(value),warning=function(w) {
    msg <- conditionMessage(w);warnings_seen[[msg]] <<- if(is.null(warnings_seen[[msg]]))1L else warnings_seen[[msg]]+1L
    invokeRestart('muffleWarning')
  })
  times[[name]] <<- as.list((proc.time()-started)[c('user.self','sys.self','elapsed')]);answer
}
pdf(file.path(out,'canonical_MERCI_diagnostic_plots.pdf'))
result <- tryCatch({
  dna <- timed('Neff_LR',upstream$Cell_Neff_cal(m,v,coverage,donors,targets,mutFeatures=features,adjust=TRUE,adjust.method='LR',frac_tran=FALSE,e=.001))
  rna <- timed('RNA_human_LOO',upstream$MERCI_LOO_MT_est(expression,targets,donors,organism='Human',epsilon=.1))
  stopifnot(setequal(rownames(dna),targets),setequal(rownames(rna),targets))
  dna <- dna[targets,,drop=FALSE];rna <- rna[targets,,drop=FALSE]
  stopifnot(all(is.finite(dna$DMTvar_count)),all(is.finite(dna$MTvar_rank)),all(is.finite(rna$Donor_MT_ind)),all(is.finite(rna$W_Donor_rank)))
  prediction <- upstream$MERCI_ReceiverPre(dna,rna,top_rank=50)
  cutoffs <- seq(100L,10L,by=-10L)
  supplied <- read.csv(threshold_table,stringsAsFactors=FALSE)
  supplied <- supplied[supplied$case=='P2' & supplied$null_id=='original',]
  stopifnot(nrow(supplied)==10L,!anyDuplicated(supplied$cutoff),setequal(supplied$cutoff,cutoffs))
  supplied <- supplied[match(cutoffs,supplied$cutoff),]
  expected <- c(3265L,2667L,2125L,1652L,1229L,872L,573L,339L,167L,56L)
  stopifnot(all(supplied$N==length(targets)),!anyNA(supplied$maximum),all(supplied$maximum==expected))
  predictions <- lapply(cutoffs,function(cut)upstream$MERCI_ReceiverPre(dna,rna,top_rank=cut))
  observed_count <- vapply(predictions,function(x)sum(x$prediction=='Receiver'),0L)
  curve <- data.frame(cutoff=cutoffs,N=length(targets),observed=observed_count,maximum=supplied$maximum,Rcm=observed_count/supplied$maximum)
  qualified <- any(observed_count>supplied$maximum)
  observed <- colSums(!is.na(m[features,targets,drop=FALSE]))
  scores <- data.frame(barcode=targets,dna,rna,raw_prediction=prediction$prediction,observed_source_feature_sites=as.integer(observed),source_marker_observable=observed>0,default_gated_candidate=qualified&prediction$prediction=='Receiver',check.names=FALSE)
  write.table(scores,file.path(out,'cells.tsv'),sep='\t',quote=FALSE,row.names=FALSE,na='NA')
  write.table(curve,file.path(out,'gate_curve.tsv'),sep='\t',quote=FALSE,row.names=FALSE,na='NA')
  saveRDS(list(dna=dna,rna=rna,features=features,curve=curve),file.path(out,'scores.rds'))
  list(status=if(qualified)'sample_signal_exploratory_candidates' else 'no_sample_level_signal_under_fixed_gate',default_qualified=qualified,raw_top50_candidates=sum(prediction$prediction=='Receiver'),default_gated_candidates=sum(scores$default_gated_candidate),default_gated_candidates_with_observable_source_feature=sum(scores$default_gated_candidate&scores$source_marker_observable),targets_without_observed_source_feature=sum(observed==0),Neff_unique_values=length(unique(dna$DMTvar_count)),largest_Neff_tie_fraction=max(table(dna$DMTvar_count))/length(targets))
},error=function(e)list(status='run_failed',error=conditionMessage(e)))
dev.off()
write_json(list(result=result,parameters=ready$parameters,source_features=features,actual_warnings=warnings_seen,phase_times=times,patient_replicates=1,scope='one untreated P2 hepatic-compatible target population, not independently malignant; RNA-derived variants not independent DNA; fixed maxima without new null draws',session=capture.output(sessionInfo())),file.path(out,'summary.json'),pretty=TRUE,auto_unbox=TRUE,na='null',digits=12)
if(identical(result$status,'run_failed'))stop(result$error,call.=FALSE)
