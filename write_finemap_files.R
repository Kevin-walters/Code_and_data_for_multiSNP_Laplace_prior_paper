################################################################################
########################### write Finemap files ################################
################################################################################

write_fm_files <- function(finemap_path, remove_rare_or_high_ld,
                           snp_ind, regout, geno_for_multi_fm,
                           num_datasets, prior_probs){
  # write finemap ld files
  for(dataset in 1:num_datasets){
    ld_for_finemap <- remove_rare_or_high_ld[[3]][snp_ind, snp_ind, dataset]
    write.table(round(ld_for_finemap, digits = 3),
                file = paste0(finemap_path, "dataset", dataset, ".ld"),
                row.names = F, col.names = F, append = F)
  }
  
  # write finemap z files
  num_cases <- num_controls <- dim(geno_for_multi_fm)[1] / 2
  final_num_snps <- dim(geno_for_multi_fm)[2]
  col1 <- NULL
  for(i in 1:final_num_snps){
    col1 <- rbind(col1, paste0("rs", i))
  }
  col2 <- rep(15, final_num_snps)
  col3 <- 1:final_num_snps
  col4 <- rep("T", final_num_snps)
  col5 <- rep("C", final_num_snps)
  col6 <- round(runif(final_num_snps, 0.01, 0.50), digits = 2)
  first6cols <- cbind(col1, col2, col3, col4, col5 ,col6)
  
  for(dataset in 1:num_datasets){
    geno_array3d <- array(geno_for_multi_fm[,, dataset],
                          dim = c(num_cases * 2, final_num_snps, 1))
    unireg <- Bayesfinemap::run_logistic_regression(geno_array3d, num_cases,
                                                    num_controls,
                                                    final_num_snps,
                                                    num_datasets = 1)
    beta_hat <- round(abs(unireg[, 1, 1]), digits = 4)
    v <- unireg[, 2, 1]
    sdev <- round(sqrt(v), digits = 4)
    final_mat <- cbind(first6cols, beta_hat, sdev)
    colnames(final_mat) <- NULL
    colnames(final_mat) <- c("rsid", "chromosome", "position", "allele1",
                             "allele2", "maf", "beta", "se")
    write.table(final_mat, file = paste0(finemap_path, "dataset", dataset, ".z"),
                row.names = F, col.names = T, append = F, quote = F)
    
  }
  
  # write finemap master file
  write.table("z;ld;k;snp;config;cred;log;n_samples",
              file = paste0(finemap_path, "master"),
              row.names = F, col.names = F, append = F, quote = F)
  
  for(i in 1:num_datasets){
    line1 <- paste0("dataset",i,".z;dataset",i,".ld;dataset.k;outfiles/dataset"
                    ,i, ".snp;outfiles/dataset",i, ".config;outfiles/dataset",i,
                    ".cred;outfiles/dataset",i,".log;",2*num_cases)
    write.table(line1,
                file = paste0(finemap_path, "master"), append = T,
                row.names = F, col.names = F, quote = F)
  }
  
  # write finemap prior probs file
  write.table(round(prior_probs, digits = 3), file = paste0(finemap_path, 
                                                            "dataset.k"),
              row.names = F, col.names = F, append = F, eol = " ")
}
