################################################################################
############################## write ROCs ######################################
################################################################################

plot_roc <- function(xlim, ylim, show_xlab = T, ldtext = "low", ss= "4000",
                     legend_x = 0.02, legend_y = 0.5, lwd = 1.7, maxk = 3,
                     linetype = 1:maxk, auc_lim, cex.axis = 1.5, cex.lab = 1.5){
  lper <- rep(NA, maxk)
  for(k in 1:maxk){
    load(paste0("G:/My Drive/Data/ppi_from_MVLap_papers/",
                   "or115or125/", ldtext,"ld/ss",ss, 
                   "/vary_max_num_causal_snps/laplace/ppi",k,".RData"))
    lppi <- ppi
    #cat(length(lppi))
    causal_snps = lppi[[3]]
    num_snps <- dim(lppi[[1]])[1]
    num_datasets <- dim(lppi[[1]])[2]
    status_vec <- rep(0, dim(lppi[[1]])[1])
    status_vec[causal_snps] <- 1
    
    status_list <- vector("list", num_datasets)
    lap_ppi_list <- vector("list", num_datasets)
    for(dataset in 1: num_datasets){
      lap_ppi_list[[dataset]] <- lppi[[1]][, dataset]
      status_list[[dataset]] <- status_vec
    }
    
    xlab <- ifelse(show_xlab, "False Positive Rate", "")
    
    lap_pred <- ROCR::prediction(lap_ppi_list, labels = status_list)
    perf <- ROCR::performance(lap_pred, measure="tpr", x.measure="fpr")
    par(cex.axis = cex.axis)
    if(k==1) {ROCR::plot(perf, lty = linetype[k], lwd = lwd, 
                        avg = "vertical", xlim = xlim, xlab = xlab,
                        cex.lab = cex.lab, 
                        ylab = "True Positive Rate")} else
                          
    ROCR::plot(perf, lty = linetype[k], lwd = lwd, avg = "vertical", add = T,
               cex.axis = cex.axis, cex.lab = cex.lab, xlab = xlab,
               ylab = "True Positive Rate")
        auc <- unlist(ROCR::performance(lap_pred,"auc")@y.values)
    lap_part_auc <- unlist(ROCR::performance(lap_pred, "auc",
                                             fpr.stop = auc_lim)@y.values)
    lpauc <- mean(lap_part_auc)/ auc_lim
    #cat("Lap pauc = ", round(lpauc, 3), "\t")
    lper[k] <- 100 * round(lpauc, 2)
  } # end of k loop
  leg_text <- paste0("K = 1; pAUC= ", lper[1])
  for(k in 2:maxk){
    leg_text <- c(leg_text, paste0("K = ", k,"; pAUC= ", lper[k]))
  }
  legend(legend_x, legend_y, 
         leg_text,
         lty = 1:maxk, cex = 1.5, bty = "n")
}

plot_roc(xlim = c(0, 0.05), 
         ylim = c(0, 1),
         legend_x = 0.025, 
         legend_y = 0.35,
         maxk= 4,auc_lim = 0.05)

