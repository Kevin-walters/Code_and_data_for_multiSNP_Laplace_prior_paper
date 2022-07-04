#' Calculate the prior probability distribution of the number of causal SNPs
#'
#' @param p The number of SNPs in the dataset.
#' @param maxk The maximum number of causal SNPs allowed in the model.
#' @param prior_prob The prior probability that a given SNP is causal.
#' @param allow_no_causal Are you allowing for the possibility that there are no
#'   causal SNPs? I.e. is the first element of \code{prior_probs} the
#'   probability that there are no causal SNPs a priori?
#' @return The prior probability distribution of the number of causal SNPs.
#' @export

calc_prior_model<-function(p, maxk, prior_prob, allow_no_causal){
  if(allow_no_causal == F){
    un_norm <- dbinom(x = 1:maxk, size = p, prob = prior_prob)
  } else un_norm <- dbinom(x = 0:maxk, size = p, prob = prior_prob)
  un_norm / sum(un_norm)
}
