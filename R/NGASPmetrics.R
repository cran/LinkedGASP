NGASPmetrics <-
function(GASP,true_output,ref_output) {
  #### validation of the GASP ####
  #Mahal_dist <- (true_output - GASP$mu)%*%solve(GASP$var,(true_output - GASP$mu))
  MSPE_base = mean((ref_output-true_output)^2)
  RMSPE_base = sqrt(MSPE_base)
  MSPE = mean((GASP$mu-true_output)^2)
  RMSPE = sqrt(MSPE) 
  #PMCC = mean((GASP$mu-true_output)^2 + diag(GASP$var))
  #RPMCC = sqrt(PMCC)
  ql = as.vector(GASP$mu+sqrt(diag(GASP$var))*qnorm(0.025))
  qu = as.vector(GASP$mu+sqrt(diag(GASP$var))*qnorm(0.975))
  coverage = true_output <= qu & true_output >= ql
  if (length(coverage) > 1) {
    if (any(coverage == FALSE)) {
      for (i in which(coverage == FALSE)) {
        coverage[i] = all.equal(true_output[i],ql[i])==TRUE || all.equal(true_output[i],qu[i])==TRUE
      }
    }
  }
  emp_cov = mean(coverage)
  length_CIs=mean(qu-ql)
  return (list(RMSPE_base = RMSPE_base,RMSPE = RMSPE, emp_cov = emp_cov, length_CIs = length_CIs, ratio = RMSPE_base/RMSPE, CIs = qu-ql))
}
