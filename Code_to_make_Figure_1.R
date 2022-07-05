################################################################################
################### function needed to calculate finemap ppi ###################
################################################################################

is_snp_present <- function(char_str, element){
  a <- unlist(strsplit(char_str, ","))
  is.element(element, a)
}


################################################################################
############################## write ROCs ######################################
################################################################################

plot_roc <- function(causal_snps, gppi, lppi, xlim, ylim, fm_stub, 
                     show_xlab = T, legend_x = 0.02, legend_y = 0.5, lwd = 1.7){
  num_snps <- dim(gppi)[1]
  num_datasets <- dim(gppi)[2]
  status_vec <- rep(0, dim(gppi)[1])
  status_vec[causal_snps] <- 1

  status_list <- vector("list", num_datasets)
  gauss_ppi_list <- vector("list", num_datasets)
  lap_ppi_list <- vector("list", num_datasets)
  for(dataset in 1: num_datasets){
    gauss_ppi_list[[dataset]] <- gppi[, dataset]
    lap_ppi_list[[dataset]] <- lppi[, dataset]
    status_list[[dataset]] <- status_vec
  }

  xlab <- ifelse(show_xlab, "False Positive Rate", "")
  gauss_pred <- ROCR::prediction(gauss_ppi_list, labels = status_list)
  perf <- ROCR::performance(gauss_pred, measure="tpr", x.measure="fpr")
  ROCR::plot(perf, xlim = xlim, ylim = ylim,
             xlab = xlab,
             ylab = "True Positive Rate",
             lty = 1, cex.lab =1.5, cex.axis = 1.8,
             lwd = lwd, xaxis.cex.axis = 1.5, yaxis.cex.axis = 1.5,
             xaxis.lwd = 1, yaxis.lwd = 1, avg = "vertical")
  #auc <- unlist(ROCR::performance(gauss_pred,"auc")@y.values)
  #(gauss_auc <- round(auc, 3))
  #mean(gauss_auc)
  gauss_part_auc <- unlist(ROCR::performance(gauss_pred, "auc",
                                             fpr.stop = 0.05)@y.values)
  #mean(gauss_part_auc)
  gpauc <- mean(gauss_part_auc)/0.05
  cat("Gauss pauc = ", round(gpauc, 3), "\t")
  cat("Plotted Gaussian ROC\n")

  lap_pred <- ROCR::prediction(lap_ppi_list, labels = status_list)
  perf <- ROCR::performance(lap_pred, measure="tpr", x.measure="fpr")
  ROCR::plot(perf, lty = 2, lwd = lwd, avg = "vertical", add = T)
  auc <- unlist(ROCR::performance(lap_pred,"auc")@y.values)
  #(lap_auc <- round(auc, 3))
  #mean(lap_auc)
  lap_part_auc <- unlist(ROCR::performance(lap_pred, "auc",
                                           fpr.stop=0.05)@y.values)
  lpauc <- mean(lap_part_auc)/ 0.05
  cat("Lap pauc = ", round(lpauc, 3), "\t")
  cat("Plotted Laplacian ROC\n")

  # fine map ROC
  # calculate PIP
  fm_ppi_list <- vector("list", num_datasets)
  for(dataset in 1:num_datasets){
    data <- read.table(paste0(fm_stub, dataset, ".config"), header = T)
    fm_pip <- rep(NA, num_snps)
    for(i in 1:num_snps){
      snp_name <- paste0("rs",i)
      snps_to_sum <- sapply(data$config, is_snp_present, element = snp_name)
      fm_pip[i] <- sum(data$prob[snps_to_sum])
    }
    fm_ppi_list[[dataset]] <- fm_pip
  }

  fm_pred <- ROCR::prediction(fm_ppi_list, labels = status_list)
  perf <- ROCR::performance(fm_pred, measure="tpr", x.measure="fpr")
  ROCR::plot(perf, lty = 3, lwd = lwd, avg = "vertical", add = T)
  #auc <- unlist(ROCR::performance(fm_pred,"auc")@y.values)
  #(fm_auc <- round(auc, 3))
  #mean(fm_auc)
  fm_part_auc <- unlist(ROCR::performance(fm_pred, "auc",
                                          fpr.stop=0.05)@y.values)
  #mean(fm_part_auc)
  fmpauc <- mean(fm_part_auc)/0.05
  cat("FM pauc = ", round(fmpauc, 3), "\t")
  cat("Plotted Finemap ROC\n")

  gper <- 100 * round(gpauc, 2)
  lper <- 100 * round(lpauc, 2)
  fmper <- 100 * round(fmpauc, 2)
  legend(legend_x, legend_y, 
         c(paste0("Gaussian: ", gper), 
         paste0("Laplace: ", lper),
         paste0("Finemap: ", fmper)),
         lty = c(1, 2, 3), cex = 1.5, bty = "n")
}

# top left plot
# path to finemap output config files
fm_stub<- paste0("G:/My Drive/Data/finemap/outfiles/or115or125/",
                          "lowld/ss4000/dataset")
load("G:/My Drive/Data/ppi_from_MVLap_papers/or115or125/lowld/ss4000/ppi.RData")
par(mfcol = c(3, 2))
par(mar = c(3,4.2,1,0.3))
plot_roc(causal_snps = ppi[[4]],
         gppi = ppi[[1]],
         lppi = ppi[[2]],
         xlim = c(0, 0.05),
         ylim = c(0, 1),
         fm_stub = fm_stub,
         show_xlab = F)



# middle left plot
fm_stub<- paste0("G:/My Drive/Data/finemap/outfiles/or115or125/",
                 "lowld/ss7000/dataset")
load("G:/My Drive/Data/ppi_from_MVLap_papers/or115or125/lowld/ss7000/ppi.RData")
par(mar = c(3,4.2,0.3,0.3))
plot_roc(causal_snps = ppi[[4]],
         gppi = ppi[[1]],
         lppi = ppi[[2]],
         xlim = c(0, 0.05),
         ylim = c(0, 1),
         fm_stub = fm_stub,
         show_xlab = F)



# bottom left plot
fm_stub<- paste0("G:/My Drive/Data/finemap/outfiles/or115or125/",
                 "lowld/ss9000/dataset")
load("G:/My Drive/Data/ppi_from_MVLap_papers/or115or125/lowld/ss9000/ppi.RData")
par(mar = c(4.2,4.2,0.3,0.3))
plot_roc(causal_snps = ppi[[4]],
         gppi = ppi[[1]],
         lppi = ppi[[2]],
         xlim = c(0, 0.05),
         ylim = c(0, 1),
         fm_stub = fm_stub,
         show_xlab = T)


# top right plot
fm_stub<- paste0("G:/My Drive/Data/finemap/outfiles/or115or125/",
                 "highld/ss4000/dataset")
load("G:/My Drive/Data/ppi_from_MVLap_papers/or115or125/highld/ss4000/ppi.RData")
par(mar = c(3, 2.5, 1, 2))
plot_roc(causal_snps = ppi[[4]],
         gppi = ppi[[1]],
         lppi = ppi[[2]],
         xlim = c(0, 0.05),
         ylim = c(0, 1),
         fm_stub = fm_stub,
         show_xlab = F,
         legend_y = 1.03)



# middle right plot
fm_stub<- paste0("G:/My Drive/Data/finemap/outfiles/or115or125/",
                 "highld/ss7000/dataset")
load("G:/My Drive/Data/ppi_from_MVLap_papers/or115or125/highld/ss7000/ppi.RData")
par(mar = c(3, 2.5, 0.3, 2))
plot_roc(causal_snps = ppi[[4]],
         gppi = ppi[[1]],
         lppi = ppi[[2]],
         xlim = c(0, 0.05),
         ylim = c(0, 1),
         fm_stub = fm_stub,
         show_xlab = F, 
         legend_y = 0.76)



# bottom right plot
fm_stub<- paste0("G:/My Drive/Data/finemap/outfiles/or115or125/",
                 "highld/ss9000/dataset")
load("G:/My Drive/Data/ppi_from_MVLap_papers/or115or125/highld/ss9000/ppi.RData")
par(mar = c(4.2, 2.5, 0.3, 2))
plot_roc(causal_snps = ppi[[4]],
         gppi = ppi[[1]],
         lppi = ppi[[2]],
         xlim = c(0, 0.05),
         ylim = c(0, 1),
         fm_stub = fm_stub,
         show_xlab = T, 
         legend_y = 0.84)

