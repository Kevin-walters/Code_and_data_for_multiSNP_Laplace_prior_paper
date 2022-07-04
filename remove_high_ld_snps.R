#' @export

remove_high_ld_SNPs <- function(geno_4array, ld_thresh){
  causal_snps <- geno_4array[[2]]
  kept_snps <- geno_4array[[3]]
  rsq_all <- geno_4array[[4]]
  num_kept <- length(kept_snps)
  num_datasets <- dim(rsq_all)[3]
  rsq_kept <- array(NA, dim = c(num_kept, num_kept, num_datasets))
  for(i in 1:num_datasets){
    rsq_kept[,, i] <- rsq_all[kept_snps, kept_snps, i]
  }
  ident <- diag(rep(1, dim(rsq_kept)[1]))
  ldsnps <- vector("list", length = num_datasets)
  for(i in 1:num_datasets){
    is_high_ld <- apply(rsq_kept[,,i] - ident, MARGIN = 1,
                        function(x) any(x > ld_thresh))
    #ldsnps[[i]] <- (1:num_kept)[is_high_ld]
    ldsnps[[i]] <- which(is_high_ld)
  }
  snps_to_go <- setdiff(unique(unlist(ldsnps)), causal_snps)
  low_ld_snps <- setdiff(1:num_kept, snps_to_go)
  cat("There are", length(low_ld_snps), "common SNPS with r2 <", ld_thresh, ";")
  geno_by_dataset <- array(NA, c(dim(geno_4array[[1]])[1], length(low_ld_snps),
                                 num_datasets))
  for(dataset in 1:num_datasets){
    geno_by_dataset[,, dataset] <- geno_4array[[1]][, low_ld_snps, dataset]
  }
  new_causal <- sapply(causal_snps, function(x) which(low_ld_snps == x))
  cat("causal SNPs are now", new_causal,"\n")
  return(list(geno_by_dataset, new_causal, rsq_kept[low_ld_snps, low_ld_snps,]))
}

