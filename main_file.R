# The three packages below are available at 
# https://github.com/Kevin-walters/
devtools::load_all("G:/My Drive/R_packages/Bayesfinemap")
devtools::load_all("G:/My Drive/R_packages/converthapgen")
devtools::load_all("G:/My Drive/R_packages/lapmapr")
source("calc_model_prior_probs.R")
source("Multivar_logreg.R")
source("remove_high_ld_snps.R")
source("select_snps_using_unibf.R")
source("create_pheno_geno.R")
source("run_multiGLM.R")
source("process_hapgen_files.R")
source("calculate_ppi.R")
source("write_finemap_files.R")
library(magrittr)

############################# Function arguments ###############################
maxk <- 1
num_datasets <- 50
allow_no_causal <- F
w <- 0.04
lambda <- 7.1
prior_prob <- 0.04
ldtext <- "low"
ss <- "4000"
causal_snps <- c(1, 312)
#causal_snps <- c(33, 95)

# path to hapgen .haps files
hapgen_readfile_path <- paste0("G:/My Drive/Data/hapgen_data/mvlap",
                               "/or115or125/",ldtext,"ld/ss",ss,"/")


# path to where ppi.RData should be stored
ppi_writefile_path <- paste0("G:/My Drive/Data/ppi_from_MVLap_papers/",
                             "or115or125/", ldtext,"ld/ss",ss, "/ppi.RData")

# path to where files needed to run finemap should be stored
fm_writefile_path <- paste0("G:/My Drive/Data/finemap/infiles/or115or125/",
                            ldtext,"ld/ss",ss, "/")




prelims <- do_prep(num_datasets = num_datasets,
                   w = w,
                   causal_snps = causal_snps,
                   maf_thresh = 0.01,
                   ld_thresh = 0.99,
                   b_thresh = 1,
                   allow_no_causal = allow_no_causal,
                   inpath = hapgen_readfile_path)
save(prelims, file = paste0(hapgen_readfile_path,"prelims.RData"))

ppi <- calc_ppi(regout = prelims[[3]],
                snp_ind = prelims[[2]],
                num_datasets = num_datasets,
                maxk = maxk,
                w = w,
                lambda = lambda,
                prior_prob = prior_prob,
                causal_snps = prelims[[5]],
                allow_no_causal_snps = allow_no_causal,
                one_both = "both")

save(ppi, file = ppi_writefile_path)

write_fm_files(finemap_path = fm_writefile_path,
               remove_rare_or_high_ld = prelims[[1]],
               snp_ind = prelims[[2]],
               regout = prelims[[3]],
               geno_for_multi_fm = prelims[[4]],
               num_datasets = num_datasets,
               prior_probs = ppi[[3]])

## Need to run finemap before plotting ROCs
