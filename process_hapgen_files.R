################################################################################
############## Process Hapgen files & do Multivar log regression ###############
################################################################################

do_prep <- function(num_datasets, w, causal_snps = c(40, 118),
                    maf_thresh = 0.01, ld_thresh = 0.99, b_thresh = 1,
                    allow_no_causal = allow_no_causal, inpath) {
  remove_rare_or_high_ld <- converthapgenr::create_geno_from_haplo_data(
    c(paste0(inpath, "scen_"),".cases.haps"),
    c(paste0(inpath, "scen_"), ".controls.haps"),
    num_datasets, maf_thresh, causal_snps) %>%
    remove_high_ld_SNPs(ld_thresh)
  #cat("Number SNPs =", dim(remove_rare_or_high_ld[[1]])[2],"\n")
  #cat("Causal SNPs are now", remove_rare_or_high_ld[[2]],"\n")
  causal_snps2 <- remove_rare_or_high_ld[[2]]
  
  high_bf_ind <- matrix(0, nrow = dim(remove_rare_or_high_ld[[1]])[2],
                        ncol = num_datasets)
  for(i in 1:num_datasets){
    high_bf_snps <- select_snps_unibf(
      geno_mat = remove_rare_or_high_ld[[1]][,,i],
      prior_form = "Gaussian", b_thresh = b_thresh, w = w)
    high_bf_ind[high_bf_snps, i] <- 1
  }
  
  cat("Finished finding high BF SNPs\n")
  high_bf_in_any_dataset <- apply(high_bf_ind, 1, function(x) any(x > 0))
  geno_for_multi_fm <- remove_rare_or_high_ld[[1]][,high_bf_in_any_dataset,]
  causal_snps3 <- sapply(causal_snps2,
                         function(x) which(which(high_bf_in_any_dataset) == x))
  cat("There are", sum(high_bf_in_any_dataset),
      "final SNPs; final causal SNPs are", causal_snps3,"\n")
  
  regout <- vector("list", num_datasets)
  cat("Starting MV logistic regression\n")
  for(dataset in 1:num_datasets){
    cat("dataset", dataset, "\n")
    regout[[dataset]] <- multivar_logreg(geno_for_multi_fm[,, dataset])
  }
  
  return(list(remove_rare_or_high_ld, high_bf_in_any_dataset, regout,
              geno_for_multi_fm, causal_snps3))
}
