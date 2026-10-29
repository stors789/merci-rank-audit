#!/usr/bin/env Rscript
# Score the one qualified source pairing using its supplied matched-N thresholds.
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=7L) stop('Expected upstream variants qualification RNA threshold-table output R-library')
up <- normalizePath(args[1],mustWork=TRUE);input <- normalizePath(args[2],mustWork=TRUE)
qual <- normalizePath(args[3],mustWork=TRUE);rna_input <- normalizePath(args[4],mustWork=TRUE)
threshold_table <- normalizePath(args[5],mustWork=TRUE);output <- args[6];out <- output
userlib <- normalizePath(args[7],mustWork=TRUE)
.Library.site <- character();.libPaths(c(userlib,.Library));suppressPackageStartupMessages(library(jsonlite))
if (dir.exists(output) || file.exists(file.path(out,'score_receipt.json'))) stop('Protect first scores')
stopifnot(identical(system2('git',c('-C',up,'rev-parse','HEAD'),stdout=TRUE),
                    '79e105aaf13159bee325c3ceceb6ee679e430341'))
qualification <- fromJSON(file.path(qual,'source_qualification.json'))
stopifnot(identical(qualification$result$status,'source_features_identified'))
data <- readRDS(file.path(input,'qualification.rds'))
e <- new.env(parent=globalenv())
for (f in c('Cell_Neff_cal.r','w_cal2.r','mut_inf_extr.r','MERCI_LOO_MT_est.r','MERCI_ReceiverPre.r')) sys.source(file.path(up,'MERCI/R',f),envir=e)
# Original sapply loses the variant row name for a single feature; retain the input feature names.
lines <- readLines(file.path(up,'MERCI/R/mut_inf_extr.r'))
before <- 't.results <- data.frame(t(t.results))'
after <- 't.results <- data.frame(t(t.results)); rownames(t.results) <- muts'
stopifnot(sum(grepl(before,lines,fixed=TRUE))==1L)
lines <- gsub(before,after,lines,fixed=TRUE)
con <- textConnection(lines); source(con,local=e); close(con)
expression <- as.matrix(read.delim(file.path(rna_input,'MT_expression.tsv'),row.names=1,check.names=FALSE))
counts <- as.matrix(read.delim(file.path(rna_input,'MT_counts.tsv'),row.names=1,check.names=FALSE))
stopifnot(all(c(data$donors,data$targets)%in%colnames(expression)),length(data$features)>0L,
          identical(dimnames(expression),dimnames(counts)),all(is.finite(expression)))
rna_nonzero <- colSums(counts[,data$targets,drop=FALSE])>0
targets <- data$targets[rna_nonzero]
if (length(data$donors)<2L || length(targets)<3L) stop('Insufficient common cells for RNA LOO')
dir.create(output)
warnings_seen <- list();started <- proc.time();phase <- list()
timed <- function(label,expr) {
  t <- proc.time(); value <- force(expr)
  phase[[label]] <<- unname((proc.time()-t)['elapsed'])
  value
}
pdf(file.path(output,'MERCI_diagnostic_plots.pdf'))
result <- withCallingHandlers(tryCatch({
  dna <- timed('Neff_LR',e$Cell_Neff_cal(data$matrix,data$variants,data$coverage,data$donors,targets,
                                        mutFeatures=data$features,adjust=TRUE,adjust.method='LR',
                                        frac_tran=FALSE,e=.001))
  rna <- timed('RNA_Human_LOO',e$MERCI_LOO_MT_est(expression,targets,data$donors,organism='Human',epsilon=.1))
  stopifnot(setequal(rownames(dna),targets),setequal(rownames(rna),targets))
  dna <- dna[targets,,drop=FALSE];rna <- rna[targets,,drop=FALSE]
  stopifnot(all(is.finite(dna$DMTvar_count)),all(is.finite(dna$MTvar_rank)),
            all(is.finite(rna$Donor_MT_ind)),all(is.finite(rna$W_Donor_rank)))
  n <- length(targets); cutoffs <- seq(100L,10L,by=-10L)
  supplied <- read.csv(threshold_table,stringsAsFactors=FALSE)
  supplied <- supplied[supplied$case=='GBML181' & supplied$null_id=='original',]
  stopifnot(nrow(supplied)==10L,!anyDuplicated(supplied$cutoff),setequal(supplied$cutoff,cutoffs))
  supplied <- supplied[match(cutoffs,supplied$cutoff),]
  expected_maxima <- c(1343L,1103L,883L,687L,522L,372L,249L,151L,80L,26L)
  stopifnot(all(supplied$N==n),!anyNA(supplied$maximum),all(supplied$maximum==expected_maxima))
  maximum <- as.integer(supplied$maximum)
  predictions <- lapply(cutoffs,function(x)e$MERCI_ReceiverPre(dna,rna,top_rank=x))
  observed <- vapply(predictions,function(x)sum(x$prediction=='Receiver'),0L)
  gate <- data.frame(N=n,cutoff=cutoffs,observed=observed,maximum=maximum,Rcm=observed/maximum)
  passed <- any(observed>maximum)
  pred50 <- predictions[[which(cutoffs==50L)]]
  feature_observed <- colSums(!is.na(data$matrix[data$features,targets,drop=FALSE]))
  scores <- data.frame(cell=targets,Neff=as.integer(dna$DMTvar_count),DNA_rank=as.integer(dna$MTvar_rank),
                       RNA_score=rna$Donor_MT_ind,RNA_rank=as.integer(rna$W_Donor_rank),
                       RNA_fraction=rna$Donor_MT_frac,source_feature_sites_observed=feature_observed,
                       observable=feature_observed>0,original_receiver_50=pred50$prediction=='Receiver',
                       check.names=FALSE)
  stopifnot(setequal(scores$DNA_rank,seq_len(n)),setequal(scores$RNA_rank,seq_len(n)))
  write.table(scores,file.path(output,'score.tsv'),sep='\t',quote=FALSE,row.names=FALSE,na='NA')
  write.csv(gate,file.path(out,'matched_N_gate.csv'),row.names=FALSE,na='NA')
  write.csv(gate[,c('N','cutoff','maximum')],file.path(output,'thresholds.csv'),row.names=FALSE)
  saveRDS(list(dna=dna,rna=rna,gate=gate),file.path(output,'score_details.rds'),compress=FALSE)
  list(status='scored', N=n, source_features=length(data$features), qualified_donors=length(data$donors),
       DNA_qualified_targets=length(data$targets), DNA_qualified_RNA_zero_targets=sum(!rna_nonzero),
       original_gate_passed=passed, raw_fixed50_receivers=sum(scores$original_receiver_50),
       default_fixed50_receivers=if(passed)sum(scores$original_receiver_50) else 0L,
       without_observed_source_feature=sum(!scores$observable), unique_Neff=length(unique(scores$Neff)),
       largest_Neff_tie=max(table(scores$Neff)), RNA_score_ties=sum(duplicated(scores$RNA_score)),
       undefined_postprocessed_RNA_fractions=sum(!is.finite(scores$RNA_fraction)),
       max_Rcm_cutoff=cutoffs[which.max(gate$Rcm)], max_Rcm=max(gate$Rcm))
},error=function(err)list(status='run_failed',error=conditionMessage(err))),warning=function(w) {
  msg <- conditionMessage(w)
  warnings_seen[[msg]] <<- if(is.null(warnings_seen[[msg]]))1L else warnings_seen[[msg]]+1L
  invokeRestart('muffleWarning')
})
dev.off()
report <- list(result=result, warnings=warnings_seen, phase_seconds=phase,
               total_seconds=unname((proc.time()-started)['elapsed']),
               threshold_source='requested table; original GBML181 N=1343 maxima checked before score output; no null draws generated',
               fixed_audit_cutoff=50, author_cutoff='maximum Rcm; not reselected in the fixed-50 audit')
write_json(report,file.path(out,'score_receipt.json'),pretty=TRUE,auto_unbox=TRUE,na='null')
cat(toJSON(report,auto_unbox=TRUE,na='null'), '\n')
if (identical(result$status,'run_failed')) quit(status=1L)
