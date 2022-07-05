devtools::load_all("G:/My Drive/R_packages/multilapmapr")
load("G:\\My Drive\\Data\\hapgen_data\\mvlap\\or115or125\\lowld\\ss4000\\prelims.RData")

ldtext <- "low"
ss <- "4000"

################################################################################
############################### one causal snp #################################
################################################################################

maxk <- 1

################################## laplace #####################################

ppi <- calc_ppi(regout = prelims[[3]],
                snp_ind = prelims[[2]],
                num_datasets = 50,
                maxk = maxk,
                w = 0.04,
                lambda = 7.1,
                prior_prob = 0.04,
                causal_snps = prelims[[5]],
                allow_no_causal_snps = F,
                one_both = "l")
ppi_writefile_path <- paste0("G:/My Drive/Data/ppi_from_MVLap_papers/",
                             "or115or125/", ldtext,"ld/ss",ss, 
                             "/vary_max_num_causal_snps/laplace/ppi",maxk,".RData")
save(ppi, file = ppi_writefile_path)

################################## gaussian ####################################

ppi <- calc_ppi(regout = prelims[[3]],
                snp_ind = prelims[[2]],
                num_datasets = 50,
                maxk = maxk,
                w = 0.04,
                lambda = 7.1,
                prior_prob = 0.04,
                causal_snps = prelims[[5]],
                allow_no_causal_snps = F,
                one_both = "g")
ppi_writefile_path <- paste0("G:/My Drive/Data/ppi_from_MVLap_papers/",
                             "or115or125/", ldtext,"ld/ss",ss, 
                             "/vary_max_num_causal_snps/gaussian/ppi",maxk,".RData")
save(ppi, file = ppi_writefile_path)

################################################################################
############################### two causal snps #################################
################################################################################

maxk <- 2

################################## laplace #####################################

ppi <- calc_ppi(regout = prelims[[3]],
                snp_ind = prelims[[2]],
                num_datasets = 50,
                maxk = maxk,
                w = 0.04,
                lambda = 7.1,
                prior_prob = 0.04,
                causal_snps = prelims[[5]],
                allow_no_causal_snps = F,
                one_both = "l")
ppi_writefile_path <- paste0("G:/My Drive/Data/ppi_from_MVLap_papers/",
                             "or115or125/", ldtext,"ld/ss",ss, 
                             "/vary_max_num_causal_snps/laplace/ppi",maxk,".RData")
save(ppi, file = ppi_writefile_path)

################################## gaussian ####################################

ppi <- calc_ppi(regout = prelims[[3]],
                snp_ind = prelims[[2]],
                num_datasets = 50,
                maxk = maxk,
                w = 0.04,
                lambda = 7.1,
                prior_prob = 0.04,
                causal_snps = prelims[[5]],
                allow_no_causal_snps = F,
                one_both = "g")
ppi_writefile_path <- paste0("G:/My Drive/Data/ppi_from_MVLap_papers/",
                             "or115or125/", ldtext,"ld/ss",ss, 
                             "/vary_max_num_causal_snps/gaussian/ppi",maxk,".RData")
save(ppi, file = ppi_writefile_path)

################################################################################
############################# three causal snps ################################
################################################################################

maxk <- 3

################################## laplace #####################################

ppi <- calc_ppi(regout = prelims[[3]],
                snp_ind = prelims[[2]],
                num_datasets = 50,
                maxk = maxk,
                w = 0.04,
                lambda = 7.1,
                prior_prob = 0.04,
                causal_snps = prelims[[5]],
                allow_no_causal_snps = F,
                one_both = "l")
ppi_writefile_path <- paste0("G:/My Drive/Data/ppi_from_MVLap_papers/",
                             "or115or125/", ldtext,"ld/ss",ss, 
                             "/vary_max_num_causal_snps/laplace/ppi",maxk,".RData")
save(ppi, file = ppi_writefile_path)

################################## gaussian ####################################

ppi <- calc_ppi(regout = prelims[[3]],
                snp_ind = prelims[[2]],
                num_datasets = 50,
                maxk = maxk,
                w = 0.04,
                lambda = 7.1,
                prior_prob = 0.04,
                causal_snps = prelims[[5]],
                allow_no_causal_snps = F,
                one_both = "g")
ppi_writefile_path <- paste0("G:/My Drive/Data/ppi_from_MVLap_papers/",
                             "or115or125/", ldtext,"ld/ss",ss, 
                             "/vary_max_num_causal_snps/gaussian/ppi",maxk,".RData")
save(ppi, file = ppi_writefile_path)

################################################################################
############################# four causal snps #################################
################################################################################

# Laplace PPIs are calculated on Sharc (Uni PC cluster)

