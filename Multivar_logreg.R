#' @export

multivar_logreg <- function(geno_mat){
  num_cases <- num_controls <- dim(geno_mat)[1] / 2
  num_snps <- dim(geno_mat)[2]
  #cat("num_snps=", num_snps, "\n")
  geno_and_pheno_data <- create_pheno_geno(geno_mat, num_cases, num_controls)
  multiGLM_output <- run_multiGLM(geno_and_pheno_data, num_snps)
  beta_hats <- matrix(multiGLM_output[[1]], ncol = 1)
  V_inv <- multiGLM_output[[2]]
  return(list(num_snps, beta_hats, V_inv))
}


