link <-
function(f1_MLEs, f2_MLEs, test_input) {
  ## f1_MLEs should be a list of MLEs of different coordinates, f2_MLEs - just MLEs themselves
  f1_MLEs = list(f1_MLEs)
  s = as.matrix(mat.or.vec(dim(test_input)[1],length(f1_MLEs)))
  s.var = as.matrix(mat.or.vec(dim(test_input)[1],length(f1_MLEs)))
  #source('Emulators.R')
  for (c in 1:length(f1_MLEs)) { 
    em = eval_type2_GASP(test_input,f1_MLEs[[c]]) # inputs are emulated with the second type GASP
    s[,c] <- em$mu
    s.var[,c] <- diag(em$var)
  } 
  
  my_input = as.matrix(f2_MLEs$data$training)
  betaM = exp(f2_MLEs$delta)
  integr = 1
  
  corr.fun <- function(training,delta) {exp(-exp(delta)*(iDist(training))^2)}
  cor = matrix(1,dim(my_input)[1],dim(my_input)[1])
  
  for (i in 1:length(betaM)) {
    integr = integr*exp(-betaM[i]*t(outer(my_input[,i],s[,i],"-"))^2/(1 + 2*s.var[,i]*betaM[i]))/sqrt(1+2*s.var[,i]*betaM[i])
    cor = cor*corr.fun(my_input[,i],f2_MLEs$delta[i])  
  }
  nu = exp(f2_MLEs$eta)/(1 + exp(f2_MLEs$eta))
  cor = cor + nu*diag(dim(my_input)[1])
  
  q <- length(f2_MLEs$basis)
  if (q > 0) {
    H = NULL
    for (i in 1:q) {
      H = cbind(H,f2_MLEs$basis[[i]](my_input))
    }
    mu = H%*%f2_MLEs$beta
  } else mu = rep(0,dim(my_input)[1])
  
  a=as.vector(solve(cor,f2_MLEs$data$fD-mu))
  mu = cbind(rep(1,dim(s)[1]),s[,1])%*%f2_MLEs$beta + integr%*%a
  
  g = solve(cor)
  
  v = rep(0,length(mu))
  for (i in 1:length(v)) {
    mat = 1
    for (j in 1:length(f1_MLEs)) {
      mat = mat*exp(-2*betaM[j]*(outer(my_input[,j],my_input[,j],"+")/2-s[i,j])^2/(1 + 4*s.var[i,j]*betaM[j])-betaM[j]*outer(my_input[,j],my_input[,j],"-")^2/2)/
        sqrt(1+4*s.var[i,j]*betaM[j])
    }
    v[i] = f2_MLEs$sigma.sq*(1 + nu - sum(g*mat)) + sum(outer(a,a,"*")*mat)
    v[i] = v[i] + 2*f2_MLEs$beta[2]*sum(integr[i,]*a*(2*betaM[1]*s.var[i,1]*my_input[,1] + s[i,1])/(1 + 2*s.var[i,1]*betaM[1]))
  }
  v = v + f2_MLEs$beta[1]^2 + 2*prod(f2_MLEs$beta)*s[,1] + f2_MLEs$beta[2]^2*(s.var[,1] + s[,1]^2) +
    2*f2_MLEs$beta[1]*integr%*%a - mu^2
  v = as.vector(v*as.numeric(v>=0))
  
  comb.em1 = list(x = test_input, mu = as.vector(mu), var = diag(as.vector(v)))
  
  #analytical form for variance of the second N-GASP
  A = cor
  AH = solve(A,H)
  HAH = t(H)%*%AH
  m = solve(HAH)
  M = solve(HAH,t(AH))
  g = AH%*%M - solve(A)
  
  sigma = rep(0,length(mu))
  sigma2 = rep(0,length(mu))
  for (i in 1:length(sigma)) {
    mat = 1
    for (j in 1:length(f1_MLEs)) {
      mat = mat*exp(-2*betaM[j]*(outer(my_input[,j],my_input[,j],"+")/2-s[i,j])^2/(1 + 4*s.var[i,j]*betaM[j])-betaM[j]*outer(my_input[,j],my_input[,j],"-")^2/2)/
        sqrt(1+4*s.var[i,j]*betaM[j])
    }
    vsp1 = exp(-betaM[1]*(my_input[,1]-s[i,1])^2/(1 + 2*s.var[i,1]*betaM[1]))/sqrt(1 + 2*s.var[i,1]*betaM[1])
    vsp3 = (2*betaM[1]*s.var[i,1]*my_input[,1] + s[i,1])*exp(-betaM[1]*(my_input[,1]-s[i,1])^2/(1 + 2*s.var[i,1]*betaM[1]))/sqrt((1 + 2*s.var[i,1]*betaM[1])^3)
    if (length(betaM)>1) {
      for (j in 2:length(betaM)) {
        vsp = exp(-betaM[j]*(my_input[,j]-s[i,j])^2/(1 + 2*s.var[i,j]*betaM[j]))/sqrt(1 + 2*s.var[i,j]*betaM[j])
        vsp1 = vsp1*vsp
        vsp3 = vsp3*vsp
      }
    }
    sigma[i] = 1 + m[1,1] + (m[1,2] + m[2,1])*s[i,1] + m[2,2]*(s.var[i,1] + s[i,1]^2) + 
      sum(g*mat) -  2*(sum(M[1,]*vsp1) + sum(M[2,]*vsp3))
    sigma2[i] = as.numeric(-t(f2_MLEs$data$fD)%*%g%*%f2_MLEs$data$fD)/(dim(H)[1]-dim(H)[2]-2)*(sigma[i]+ nu) + sum(outer(a,a,"*")*mat)
    
    sigma[i] = f2_MLEs$sigma.sq*(sigma[i]+ nu) + sum(outer(a,a,"*")*mat)
    
    sigma[i] = sigma[i] + 2*f2_MLEs$beta[2]*sum(integr[i,]*a*(2*betaM[1]*s.var[i,1]*my_input[,1] + s[i,1])/(1 + 2*s.var[i,1]*betaM[1]))
    
    sigma2[i] = sigma2[i] + 2*f2_MLEs$beta[2]*sum(integr[i,]*a*(2*betaM[1]*s.var[i,1]*my_input[,1] + s[i,1])/(1 + 2*s.var[i,1]*betaM[1]))
  }
  sigma = sigma + f2_MLEs$beta[1]^2 + 2*prod(f2_MLEs$beta)*s[,1] + f2_MLEs$beta[2]^2*(s.var[,1] + s[,1]^2) +
    2*f2_MLEs$beta[1]*integr%*%a - mu^2
  sigma = as.vector(sigma)
  sigma = sigma*as.numeric(sigma>=0)
  comb.em2 = list(x = test_input, mu = as.vector(mu), var = diag(as.vector(sigma)))
  
  sigma2 = sigma2 + f2_MLEs$beta[1]^2 + 2*prod(f2_MLEs$beta)*s[,1] + f2_MLEs$beta[2]^2*(s.var[,1] + s[,1]^2) +
    2*f2_MLEs$beta[1]*integr%*%a - mu^2
  sigma2 = as.vector(sigma2)
  sigma2 = sigma2*as.numeric(sigma2>=0)
  comb.em3 = list(x = test_input, mu = as.vector(mu), var = diag(as.vector(sigma2)))
  comb.emT = list(x = test_input, mu = as.vector(mu), sigma = diag(as.vector(sigma2)),df = dim(H)[1]-dim(H)[2])
  return(list(em1 = comb.em1, em2 = comb.em2, em3 = comb.em3, emT = comb.emT))
}
