eval_GASP_RFP <-
function(data,basis,corr.cols,nugget = F) {
  training <- as.data.frame(data$training)
  smooth <- data$smooth
  fD <- data$fD
  q <- length(basis)
  if (q > 0) H <- mapply(function(cnt) {basis[[cnt]](data$training)}, as.list(1:q), SIMPLIFY = T)

  corr.fun <- function(tr,d,s) {exp(-exp(d)*(iDist(tr))^s)}

 # library(spBayes)
 # library(nloptr)

  n <- dim(training)[1]
  p <- length(corr.cols)
  if (p > 0) {
    if (q > 0) {
      loglik <- function(par, n, q, nugget) {
        if (nugget) {
          delta <- par[1:(length(par)-1)]
          eta <- par[length(par)]
        } else {
          delta <- par
          eta <- -Inf
        }
        nu <- exp(eta)/(1+exp(eta))

        A <- mapply(corr.fun, training[,corr.cols, drop = F], delta, smooth, SIMPLIFY = FALSE)
        A <- Reduce("*", A)
        A <- A + nu*diag(n)

        AH <- solve(A,H)
        HAH <- crossprod(H,AH)
        Q <- solve(A) - AH%*%solve(HAH,t(AH))

        if (nugget) {
          Wb <- mapply(function(tr,s) {-((A- nu*diag(n))*iDist(tr)^s)%*%Q}, training[,corr.cols, drop = F], smooth,SIMPLIFY = FALSE)
          Wb[[p+1]]<-Q

          IR <- matrix(n-q, p+2, p+2)
          IR[1,2:(p+2)] <- IR[2:(p+2),1]<- mapply(function(l) {sum(diag(l))}, Wb, SIMPLIFY = TRUE)

          for(i in 1:(p+1)) {
            IR[i+1,2:(p+2)] <- IR[2:(p+2),i+1] <- mapply(function(l) {sum(diag(Wb[[i]]%*%l))}, Wb,SIMPLIFY = TRUE)
          }

          return(sum(delta) + log(nu) -.5*determinant(A)$mod[1] -.5*determinant(HAH)$mod[1] - (n - q)/2 * log(crossprod(fD,(Q%*%fD)))+.5*determinant(IR)$mod[1])

        } else {
          Wb <- mapply(function(tr,s) {-(A*iDist(tr)^s)%*%Q}, training[,corr.cols, drop = F], smooth,SIMPLIFY = FALSE)

          IR <- matrix(n-q, p+1, p+1)
          IR[1,2:(p+1)] <- IR[2:(p+1),1]<- mapply(function(l) {sum(diag(l))}, Wb, SIMPLIFY = TRUE)

          for(i in 1:p) {
            IR[i+1,2:(p+1)] <- IR[2:(p+1),i+1] <- mapply(function(l) {sum(diag(Wb[[i]]%*%l))}, Wb,SIMPLIFY = TRUE)
          }

          return(sum(delta) -.5*determinant(A)$mod[1] -.5*determinant(HAH)$mod[1] - (n - q)/2 * log(crossprod(fD,(Q%*%fD)))+.5*determinant(IR)$mod[1])
        }
      }
    } else {
      loglik <- function(par, n, q, nugget) {
        if (nugget) {
          delta <- par[1:(length(par)-1)]
          eta <- par[length(par)]
        } else {
          delta <- par
          eta <- -Inf
        }

        nu <- exp(eta)/(1+exp(eta))

        A <- mapply(corr.fun, training[,corr.cols, drop = F], delta, smooth, SIMPLIFY = FALSE)
        A <- Reduce("*", A)
        A <- A + nu*diag(n)

        Q <- solve(A)

        if (nugget) {
          Wb <- mapply(function(tr,s) {-((A- nu*diag(n))*iDist(tr)^s)%*%Q}, training[,corr.cols, drop = F], smooth,SIMPLIFY = FALSE)
          Wb[[p+1]]<-Q

          IR <- matrix(n, p+2, p+2)
          IR[1,2:(p+2)] <- IR[2:(p+2),1]<- mapply(function(l) {sum(diag(l))}, Wb, SIMPLIFY = TRUE)

          for(i in 1:(p+1)) {
            IR[i+1,2:(p+2)] <- IR[2:(p+2),i+1] <- mapply(function(l) {sum(diag(Wb[[i]]%*%l))}, Wb,SIMPLIFY = TRUE)
          }

          return(sum(delta) + log(nu) -.5*determinant(A)$mod[1] - n/2 * log(crossprod(fD,(Q%*%fD)))+.5*determinant(IR)$mod[1])

        } else {
          Wb <- mapply(function(tr,s) {-(A*iDist(tr)^s)%*%Q}, training[,corr.cols, drop = F],smooth,SIMPLIFY = FALSE)

          IR <- matrix(n, p+1, p+1)
          IR[1,2:(p+1)] <- IR[2:(p+1),1]<- mapply(function(l) {sum(diag(l))}, Wb, SIMPLIFY = TRUE)

          for(i in 1:p) {
            IR[i+1,2:(p+1)] <- IR[2:(p+1),i+1] <- mapply(function(l) {sum(diag(Wb[[i]]%*%l))}, Wb,SIMPLIFY = TRUE)
          }

          #return(sum(delta) -.5*log(det(A)) - n/2 * log(crossprod(fD,(Q%*%fD)))+.5*log(det(IR)))
          return(sum(delta) -.5*determinant(A)$mod[1] - n/2 * log(crossprod(fD,(Q%*%fD)))+.5*determinant(IR)$mod[1])
        }
      }
    }
    start.range <- apply(apply(training[, corr.cols, drop = F], 2,range), 2, diff)/2
    start.delta <- -smooth * log(start.range) #starting values for optimization of integrated likelihood are half the range of the domains of inputs
    if (nugget) {
      v <- var(lm(as.data.frame(cbind(fD, training[, corr.cols, drop = F])))$residuals)
      start.nu <- v/(var(lm(fD~1)$residuals) - v)
      if (start.nu>1) start.nu <- 1/start.nu
      if (!start.nu) start.nu <- .00001
      start.eta <- log(start.nu/(1-start.nu))
      max.loglik <- try(optim(c(start.delta, start.eta), function(par) {-loglik(par, n, q, nugget)}), silent = T)
      while (!("par" %in% names(max.loglik)) & start.nu<1) {
        max.loglik <- try(optim(c(start.delta, start.eta), function(par) {-loglik(par, n, q, nugget)}), silent = T)
        start.nu = start.nu*10
        start.eta <- log(start.nu/(1-start.nu))
      }
      delta <- max.loglik$par[1:(length(max.loglik$par)-1)]
      eta <- max.loglik$par[length(max.loglik$par)]
    } else {
      if (length(start.delta)<2) {
        #max.loglik <- try(optim(start.delta, lower = start.delta - 100, upper = start.delta + 100,function(par) {-loglik(par, n, q, nugget)},method = "Brent"), silent = T)

        max.loglik <- nloptr(start.delta,eval_f = function(par) {-loglik(par, n, q, nugget)},opts = list( "algorithm" = "NLOPT_LN_COBYLA","xtol_rel"  = 1.0e-8))
      } else {
        #max.loglik <- try(optim(start.delta, function(par) {-loglik(par, n, q, nugget)}), silent = T)
        max.loglik <- nloptr(start.delta,eval_f = function(par) {-loglik(par, n, q, nugget)},opts = list( "algorithm" = "NLOPT_LN_COBYLA","xtol_rel"  = 1.0e-8))

      }

      delta <- max.loglik$solution
      eta <- -Inf
    }
    nu <- exp(eta)/(1+exp(eta))

    A <- mapply(corr.fun, training[, corr.cols, drop = F], delta, smooth, SIMPLIFY = FALSE)
    A <- Reduce("*", A)
    A <- A + nu*diag(n)

    if (q>0) {
      AH <- solve(A, H)
      HAH <- crossprod(H, AH)
      sigma.sq <- as.numeric(crossprod(fD, solve(A) - AH%*%solve(HAH, t(AH)))%*%fD/(n-q+2))
      beta <- solve(HAH, t(AH))%*%fD
      return(list(beta = beta, delta = delta, eta = eta, sigma.sq = sigma.sq, data = data, nugget = nugget, basis=basis, corr.cols = corr.cols))
    } else {
      sigma.sq <- as.numeric(crossprod(fD,(solve(A,fD)))/(n+2))
      return(list(delta = delta, eta = eta, sigma.sq = sigma.sq, data = data, nugget = nugget, basis=basis, corr.cols = corr.cols))
    }
  } else {
    if (q > 0) {
      HH <- crossprod(H)
      sigma.sq <- as.numeric(crossprod(fD,(diag(n) - H%*%solve(HH,t(H))))%*%fD/(n-q+2))
      beta <- solve(HH,t(H))%*%fD
      return(list(beta = beta, eta = -Inf, sigma.sq = sigma.sq, data = data, nugget = FALSE, basis = basis, corr.cols = corr.cols))
    } else {
      sigma.sq <- as.numeric(crossprod(fD))/(n+2)
      return(list(eta = -Inf, sigma.sq = sigma.sq, data = data, nugget = FALSE, basis=basis, corr.cols = corr.cols))
    }
  }
}
