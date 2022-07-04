create_pheno_geno <- function(geno_data, num_cases, num_controls){
  pheno <- c(rep(1, num_cases), rep(0, num_controls))
  genopheno <- cbind(pheno, geno_data[, ])
  return(genopheno)
}
