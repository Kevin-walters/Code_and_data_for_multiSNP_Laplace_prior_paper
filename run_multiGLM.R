run_multiGLM <- function(genopheno, p){
  PredictorVariables <- paste("genopheno[,",2:(p + 1), "]", sep = "")
  Formula <- as.formula(paste("genopheno[, 1] ~ ",
                              paste(PredictorVariables, collapse = " + ")))
  model <- glm(Formula, family = binomial, maxit = 100)
  beta_hat <- coefficients(model)[-1]
  beta_hats <- na.omit(beta_hat)
  if(any(is.na(beta_hat))) cat("some effect sizes are NAs")
  var_beta_hats <- vcov(model)[-1, -1]
  return(list(beta_hats, solve(var_beta_hats)))
}
