eval_type2_GASP <-
function(input,GASPparams) {
  if (!("beta" %in% names(GASPparams))) return(eval_type1_GASP(input,GASPparams))
  else{
    beta <- GASPparams$beta
    #delta <- GASPparams$delta
    basis <- GASPparams$basis
    corr.cols <- GASPparams$corr.cols
    eta <- GASPparams$eta
    nu <- exp(eta)/(1+exp(eta))
    sigma.sq <- GASPparams$sigma.sq
    training <- as.matrix(GASPparams$data$training)
    q <- length(basis)
    H = NULL
    Hp = NULL
    for (i in 1:q) {
      H = cbind(H,basis[[i]](training))
      Hp = cbind(Hp,basis[[i]](input))
    }
    n <- dim(training)[1]
    p = length(corr.cols)
    
    if (!("delta" %in% names(GASPparams))) {
      H = NULL
      Hp = NULL
      for (i in 1:q) {
        H = cbind(H,basis[[i]](training))
        Hp = cbind(Hp,basis[[i]](input))
      }
      
      beta <- GASPparams$beta
      return(list(x = input, mu = Hp%*%beta,var = sigma.sq*(diag(dim(input)[1]) + Hp%*%solve(t(H)%*%H,t(Hp)))))     
    } else {
      delta <- GASPparams$delta
      
      data<-GASPparams$data
      smooth <- data$smooth
      fD<-data$fD
      
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
      
      AH <- solve(A,H)
      HAH <- crossprod(H,AH)
      
      addM = Hp - corr.tr.input%*%AH
      cor = addM%*%solve(HAH,t(addM))+cor-corr.tr.input%*%solve(A,t(corr.tr.input))
      cor = (cor + t(cor))/2
      cor = cor + nu*diag(dim(input)[1])
      
      diag(cor) = diag(cor)*(diag(cor)>=0)
      mu = corr.tr.input%*%solve(A,fD - H%*%beta) + Hp%*%beta
      
      return(list(x = input, mu = mu, sigma.sq = sigma.sq, cor = cor, var = sigma.sq*cor))
    }
  }
}
