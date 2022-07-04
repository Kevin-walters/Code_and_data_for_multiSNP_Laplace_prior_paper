#' @export

select_snps_unibf <- function(geno_mat, prior_form, target_num_snps = NULL,
                              b_thresh = NULL, w = NULL, lambda = NULL){
  num_cases <- num_controls <- dim(geno_mat)[1] / 2
  num_snps <- dim(geno_mat)[2]
  #geno_array3d <- array(data = NA, dim = c(2 * num_cases, num_snps, 1 ))
  #geno_array3d[,, 1] <- geno_mat
  geno_array3d <- array(geno_mat, dim = c(dim(geno_mat)[1], dim(geno_mat)[2], 1))

  # perform univar log reg
  unireg <- Bayesfinemap::run_logistic_regression(geno_array3d, num_cases,
                                    num_controls, num_snps, num_datasets = 1)
  if(prior_form == "Gaussian"){
    bf <- Bayesfinemap::calc_wbf(beta_hat = unireg[, 1, 1], v = unireg[, 2, 1], w)
  }
  if(prior_form == "Laplace"){
    bf <- Bayesfinemap::calc_lbf(beta_hat = unireg[, 1, 1], v = unireg[, 2, 1], lambda)
  }
  if(is.null(target_num_snps) & is.null(b_thresh)){
    cat("Error: must specify one of target_num_snps or b_thresh")
    return()
  }
  if(!is.null(target_num_snps) & !is.null(b_thresh)){
    cat("Error: cannot specify both target_num_snps and b_thresh")
    return()
  }
  if(!is.null(b_thresh)) {
    #new_geno <- geno_array3d[, bf > b_thresh, 1]
    #new_geno <- geno_mat[, bf > b_thresh]
    retained_snps <- sort(which(bf > b_thresh))
    #new_causal <- which(retained_snps == causal_snp)
  }
  if(!is.null(target_num_snps)){
    retained_snps <- sort(order(bf, decreasing = T)[1:target_num_snps])
    #new_geno <- geno_array3d[, retained_snps, 1]
    #new_geno <- geno_mat[, retained_snps]
    #new_causal <- which(retained_snps == causal_snp)
  }
  #return(list(new_geno, new_causal))
  return(retained_snps)
}
