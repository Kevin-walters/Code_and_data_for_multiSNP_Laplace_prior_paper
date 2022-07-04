################################################################################
############################# Run MVlap/norm ppi ###############################
################################################################################

calc_ppi <- function(regout, snp_ind, num_datasets, maxk, w,
                     lambda, prior_prob, causal_snps, allow_no_causal_snps,
                     one_both){
  if(!is.element(one_both, c("l","g","both"))){
    stop("Argument one_both must be g, l or both\n")
  }  
  final_num_snps <- sum(snp_ind)
  gauss_ppi <- matrix(NA, nrow = final_num_snps, ncol = num_datasets)
  lap_ppi <- matrix(NA, nrow = final_num_snps, ncol = num_datasets)
  
  prior_probs <- calc_prior_model(p = final_num_snps, maxk = maxk,
                                  prior_prob = prior_prob,
                                  allow_no_causal = allow_no_causal_snps)
  
  for(dataset in 1:num_datasets){
    cat("dataset = ", dataset, "\n")
    if(one_both == "g" || one_both == "both"){
      out_gauss <- lapmapr::gauss_mv(num_snps = regout[[dataset]][[1]],
                                     beta_hats = regout[[dataset]][[2]],
                                     V_inv = regout[[dataset]][[3]],
                                     maxk = maxk,
                                     w = w,
                                     prior_probs = prior_probs,
                                     allow_no_causal = allow_no_causal_snps)
      gauss_ppi[, dataset] <- out_gauss[[1]]
      cat("Finished Gaussian\n")
    }
    if(one_both == "l" || one_both == "both"){
      out_lap <- lapmapr::laplace_mv(num_snps = regout[[dataset]][[1]],
                                     beta_hats = regout[[dataset]][[2]],
                                     V_inv = regout[[dataset]][[3]],
                                     maxk = maxk,
                                     lambda = lambda,
                                     prior_probs = prior_probs,
                                     allow_no_causal = allow_no_causal_snps)
      lap_ppi[, dataset] <- out_lap[[1]]
      cat("Finished Laplace\n")
    }
  }
  #cat("Finished Gaussian & Laplace prior analysis\n")
  if(one_both == "g"){
    return(list(gauss_ppi, prior_probs, causal_snps))
  }  
  if(one_both == "l"){
    return(list(lap_ppi, prior_probs, causal_snps))
  }  
  if(one_both == "both"){
    return(list(gauss_ppi, lap_ppi, prior_probs, causal_snps))
  }  
}
