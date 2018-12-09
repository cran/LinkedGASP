eval_type1_GASP <-
function(input,GASPparams) {
  delta <- GASPparams$delta
  basis <- GASPparams$basis
  corr.cols <- GASPparams$corr.cols
  sigma.sq <- GASPparams$sigma.sq
  data <- GASPparams$data
  eta <- GASPparams$eta
  nu <- exp(eta)/(1+exp(eta))
  smooth <- data$smooth
  fD <- data$fD
  training <- as.matrix(GASPparams$data$training)
  corr.fun <- function(training1,training2,delta,smooth) {if (missing(training2)) exp(-exp(delta)*(iDist(training1))^smooth) else
    exp(-exp(delta)*(iDist(training1,training2))^smooth)}
  
  n <- dim(training)[1]
  p <- length(corr.cols)
  q <- length(basis)
  
  if (!("delta" %in% names(GASPparams))) {
    if (q > 0) {
      H <- NULL
      Hp <- NULL
      for (i in 1:q) {
        H <- cbind(H,basis[[i]](training))
        Hp <- cbind(Hp,basis[[i]](input))
      }
      
      beta <- GASPparams$beta
      return(list(x = input, mu = Hp%*%beta,var = sigma.sq*diag(dim(input)[1])))
    } else return(list(x = input, mu = rep(0,dim(input)[1]),var = sigma.sq*diag(dim(input)[1])))
  } else {
    delta <- GASPparams$delta
    
    A <- corr.fun(training[,corr.cols[1]],delta = delta[1],smooth = smooth[1])
    corr.tr.input = corr.fun(input[,corr.cols[1]],training[,corr.cols[1]],delta[1],smooth[1])
    cor = corr.fun(input[,corr.cols[1]],delta = delta[1],smooth = smooth[1])
    
    if (p > 1) for (i in 2:p) {
      A = A*corr.fun(training[,corr.cols[i]],delta = delta[i],smooth = smooth[i])
      corr.tr.input = corr.tr.input*corr.fun(input[,corr.cols[i]],training[,corr.cols[i]],delta[i],smooth[i])
      cor = cor*corr.fun(input[,corr.cols[i]],delta = delta[i],smooth = smooth[i])
    }
    A = A + nu*diag(n)
    
    cor = cor + nu*diag(dim(input)[1]) - corr.tr.input%*%solve(A,t(corr.tr.input))
    cor = (cor + t(cor))/2
    diag(cor) = diag(cor)*(diag(cor)>=0)
    
    if (q > 0) {
      H = NULL
      Hp = NULL
      for (i in 1:q) {
        H = cbind(H,basis[[i]](training))
        Hp = cbind(Hp,basis[[i]](input))
      }
      
      beta <- GASPparams$beta
      return(list(x = input, mu = Hp%*%beta + corr.tr.input%*%solve(A,fD-H%*%beta),
                  var = sigma.sq*cor))
    } else return(list(x = input, mu = corr.tr.input%*%solve(A,fD),var = sigma.sq*cor))
  }
}
