eval_TGASP <-
function(input,GASPparams) {
  # This is the third type of the emulator, which is a T-distribution
  beta <- GASPparams$beta
  #delta <- GASPparams$delta
  basis <- GASPparams$basis
  corr.cols <- GASPparams$corr.cols
  eta <- GASPparams$eta
  nu <- exp(eta)/(1+exp(eta))
  training <- as.matrix(GASPparams$data$training)
  n <- dim(training)[1]
  p = length(corr.cols)
  q <- length(basis)
  xi.sq <- GASPparams$sigma.sq*(n-q+2)/(n-q)
  data<-GASPparams$data
  smooth <- data$smooth
  fD<-data$fD
  
  if (!("delta" %in% names(GASPparams))) {
    if (!("beta" %in% names(GASPparams))) {
      return(list(x = input, mu = rep(0,dim(input)[1]), sigma = xi.sq*(diag(dim(input)[1])), df = n - q, xi.sq = xi.sq))
    } else {
      H = NULL
      Hp = NULL
      for (i in 1:q) {
        H = cbind(H,basis[[i]](training))
        Hp = cbind(Hp,basis[[i]](input))
      }
      return(list(x = input, mu = Hp%*%beta, sigma = xi.sq*(diag(dim(input)[1]) + Hp%*%solve(t(H)%*%H,t(Hp))), df = n - q, xi.sq = xi.sq))
    }
    
    
    H = NULL
    Hp = NULL
    for (i in 1:q) {
      H = cbind(H,basis[[i]](training))
      Hp = cbind(Hp,basis[[i]](input))
    }
    
    beta <- GASPparams$beta
    return(list(x = input, mu = Hp%*%beta,var = GASPparams$sigma.sq*(diag(dim(input)[1]) + Hp%*%solve(t(H)%*%H,t(Hp)))))     
  } else {
    delta <- GASPparams$delta
    corr.fun <- function(training1,training2,delta,smooth) {if (missing(training2)) exp(-exp(delta)*(iDist(training1))^smooth) else
      exp(-exp(delta)*(iDist(training1,training2))^smooth)}
    
    
    A <- corr.fun(training[,corr.cols[1]],delta = delta[1],smooth = smooth[1])
    corr.tr.input = corr.fun(input[,corr.cols[1]],training[,corr.cols[1]],delta[1],smooth[1])
    cor = corr.fun(input[,corr.cols[1]],delta = delta[1],smooth = smooth[1])
    
    if (p > 1) for (i in 2:p) {
      A = A*corr.fun(training[,corr.cols[i]],delta = delta[i],smooth = smooth[i])
      corr.tr.input = corr.tr.input*corr.fun(input[,corr.cols[i]],training[,corr.cols[i]],delta[i],smooth[i])
      cor = cor*corr.fun(input[,corr.cols[i]],delta = delta[i],smooth = smooth[i])
    }  
    A = A + nu*diag(n)
    
    if (!("beta" %in% names(GASPparams))) {
      sigma = cor - corr.tr.input%*%solve(A,t(corr.tr.input))
      sigma = (sigma + t(sigma))/2
      mu = as.vector(corr.tr.input%*%solve(A,fD))
    } else {
      H = NULL
      Hp = NULL
      for (i in 1:q) {
        H = cbind(H,basis[[i]](training))
        Hp = cbind(Hp,basis[[i]](input))
      }
      
      AH <- solve(A,H)
      HAH <- crossprod(H,AH)
      
      addM = Hp - corr.tr.input%*%AH
      sigma = addM%*%solve(HAH,t(addM))+cor- 
        corr.tr.input%*%solve(A,t(corr.tr.input))
      sigma = (sigma + t(sigma))/2
      mu = as.vector(corr.tr.input%*%solve(A,fD - H%*%beta) + Hp%*%beta)
    }
    diag(sigma) = diag(sigma)*(diag(sigma)>=0)
    sigma = sigma + nu*diag(dim(input)[1])
    return(list(x = input, mu = mu, sigma = xi.sq*sigma, df = n - q, xi.sq = xi.sq))
  }
}
