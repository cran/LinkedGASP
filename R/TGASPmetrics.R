TGASPmetrics <-
  function (TGASP, true_output, ref_output)
  {
    MSPE_base = mean((ref_output - true_output)^2)
    RMSPE_base = sqrt(MSPE_base)
    MSPE = mean((TGASP$mu - true_output)^2)
    RMSPE = sqrt(MSPE)
    #PMCC = mean((TGASP$mu - true_output)^2 + diag(TGASP$sigma) * TGASP$df/(TGASP$df - 2))
    #RPMCC = sqrt(PMCC)
    ql = as.vector(TGASP$mu + sqrt(diag(TGASP$sigma)) * qt(0.025,df = TGASP$df))
    qu = as.vector(TGASP$mu + sqrt(diag(TGASP$sigma)) * qt(0.975,df = TGASP$df))
    coverage = true_output <= qu & true_output >= ql
    if (length(coverage) > 1) {
      if (any(coverage == FALSE)) {
        for (i in which(coverage == FALSE)) {
          coverage[i] = all.equal(true_output[i], ql[i]) == TRUE || all.equal(true_output[i], qu[i]) == TRUE
        }
      }
    }
emp_cov = mean(coverage)
length_CIs = mean(qu - ql)
return(list(RMSPE_base = RMSPE_base, RMSPE = RMSPE,
            emp_cov = emp_cov, length_CIs = length_CIs, ratio = RMSPE_base/RMSPE,
            CIs = qu - ql))
  }