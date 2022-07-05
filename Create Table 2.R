table1 <- matrix(NA, nrow = 12, ncol = 8)
for(k in 1:4){
  load(paste0("G:/My Drive/Data/ppi_from_MVLap_papers/or115or125/lowld/",
  "ss4000/vary_max_num_causal_snps/laplace/ppi",k,".RData"))
  causal_snps = ppi[[3]]
  prior_exp <- sum(ppi[[2]] * 1:k)
  ppi <- ppi[[1]]
  post_exp <- sum(ppi)/50
  cs1ppi <- ppi[causal_snps[1], ]
  cs2ppi <- ppi[causal_snps[2], ]
  cs1rank <- apply(-ppi, 2, FUN = rank)[causal_snps[1], ]
  cs2rank <- apply(-ppi, 2, FUN = rank)[causal_snps[2], ]
  table1[1, (2*k-1)] <- prior_exp
  table1[2, (2*k-1)] <- post_exp
  table1[3, (2*k-1)] <- quantile(cs1ppi, 0.50)
  table1[4, (2*k-1)] <- quantile(cs1ppi, 0.10)
  table1[5, (2*k-1)] <- quantile(cs1rank, 0.50)
  table1[6, (2*k-1)] <- quantile(cs1rank, 0.90)
  table1[7, (2*k-1)] <- max(cs1rank)
  table1[8, (2*k-1)] <- quantile(cs2ppi, 0.50)
  table1[9, (2*k-1)] <- quantile(cs2ppi, 0.10)
  table1[10, (2*k-1)] <- quantile(cs2rank, 0.50)
  table1[11, (2*k-1)] <- quantile(cs2rank, 0.90)
  table1[12, (2*k-1)] <- max(cs2rank)
  
  load(paste0("G:/My Drive/Data/ppi_from_MVLap_papers/or115or125/lowld/",
              "ss4000/vary_max_num_causal_snps/gaussian/ppi",k,".RData"))
  causal_snps = ppi[[3]]
  prior_exp <- sum(ppi[[2]] * 1:k)
  ppi <- ppi[[1]]
  post_exp <- sum(ppi)/50
  cs1ppi <- ppi[causal_snps[1], ]
  cs2ppi <- ppi[causal_snps[2], ]
  cs1rank <- apply(-ppi, 2, FUN = rank)[causal_snps[1], ]
  cs2rank <- apply(-ppi, 2, FUN = rank)[causal_snps[2], ]
  table1[1, (2*k)] <- prior_exp
  table1[2, (2*k)] <- post_exp
  table1[3, (2*k)] <- quantile(cs1ppi, 0.50)
  table1[4, (2*k)] <- quantile(cs1ppi, 0.10)
  table1[5, (2*k)] <- quantile(cs1rank, 0.50)
  table1[6, (2*k)] <- quantile(cs1rank, 0.90)
  table1[7, (2*k)] <- max(cs1rank)
  table1[8, (2*k)] <- quantile(cs2ppi, 0.50)
  table1[9, (2*k)] <- quantile(cs2ppi, 0.10)
  table1[10, (2*k)] <- quantile(cs2rank, 0.50)
  table1[11, (2*k)] <- quantile(cs2rank, 0.90)
  table1[12, (2*k)] <- max(cs2rank)
  table1 <- round(table1, digits = 2)
}
knitr::kable(table1, "latex")
# not used
#non_c_meanppi <- meanppi[-c(1, 65)]
#uq <- round(quantile(non_c_meanppi, 0.95), digits = 3)
#top3_meanppi <- mean(sort(non_c_meanppi, decreasing = T)[1:3])
#table1[3:7, k] <- sort(meanppi, decreasing = T)[1:5]
