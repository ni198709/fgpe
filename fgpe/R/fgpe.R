#' The graphical projection estimator (GPE) based on the partially penalized regression
#' @title simple function
#' @description today,I create my first function,a very usrful function.
#' @param n n is the number of observations
#' @param p p is the dimension of graphical model
#' @param mean mean is a mean vector,
#' @param sigma sigma is a covariance matrix
#' @param length length is a tuning parameter of HOLP or SIS
#' @return A list of
#' \itemize{
#'  \item{\code{gamma_hat}: estimator for regression coefficient matrix}
#'  \item{\code{p_valu}: p_value matrix of Hypothesis testing}
#' }
#' @export
fgpe <-
  function(n, p, mean, sigma, length) {
    colsum = c(rep(0, p))
    I = diag(1, n, n)
    p_valu = matrix(0, p - 1, p)
    gamma_hat = matrix(0, p - 1, p)
    X = matrix(0, n, p)
    X = mvrnorm(n = n, mean, sigma)
    
    for (j in 1:p) {
      Xjbu = X[, -j]
      Xj   = X[, j]
      gujilingling = scalreg(Xjbu, Xj, lam0 = lambda.univ)
      gammajini = coef(gujilingling)
      tausj = gujilingling$hsigma
      taus = tausj
      
      #SIS
      #omegaj=t(Xjbu)%*%Xj
      #Holp
      omegaj = t(Xjbu) %*% (solve(Xjbu %*% t(Xjbu))) %*% Xj
      BICj1 = c(rep(0, length))
      min = 2
      for (sj in 2:length) {
        Sj = sort(order(abs(omegaj), decreasing = TRUE)[1:sj])
        xj = Xjbu[, Sj] %*% solve(t(Xjbu[, Sj]) %*% Xjbu[, Sj]) %*% t(Xjbu[, Sj]) %*%
          Xj
        residualj = Xj - xj
        resj = sum(residualj^2)
        
        BICj1[sj] = resj / (n) + sj * (taus^2) * log(p) * log(log(n)) / n
        if (BICj1[sj] <= BICj1[min]) {
          min = sj
        }
      }
      #Min[k] = min
      Sj.hat = sort(order(abs(omegaj), decreasing = TRUE)[1:min])
      Xj2 = Xjbu[, Sj.hat]
      Pj = Xj2 %*% solve(t(Xj2) %*% Xj2) %*% t(Xj2)
      
      Sjbu = setdiff(1:(p - 1), Sj.hat)
      lamdazheng = (2 / n * log(length(Sjbu)))^0.5
      
      psij = matrix(0, n, p - 1)
      zj = matrix(0, n, p - 1)
      zetaj = c(rep(0, p - 1))
      gamma_hatj = c(rep(0, p - 1))
      p_valuj = c(rep(0, p - 1))
      
      for (i in 1:(p - 1)) {
        if (i %in% Sj.hat) {
          Xj3 = Xjbu[, setdiff(Sj.hat, i)]
          #psij= (I-Xj3%*%solve(t(Xj3)%*%Xj3)%*%t(Xj3))%*%Xjbu
          psij[, Sjbu] = (I - Xj3 %*% solve(t(Xj3) %*% Xj3) %*% t(Xj3)) %*% Xjbu[, Sjbu]
          psij[, i] = (I - Xj3 %*% solve(t(Xj3) %*% Xj3) %*% t(Xj3)) %*% Xjbu[, i]
          fitj = glmnet(psij[, Sjbu], psij[, i], intercept = 0)
          keaij = coef(fitj)[-1, ]
          meihaoj = ncol(keaij)
          BIC2j = c(rep(0, meihaoj))
          min1 = 1
          for (sj1 in 1:meihaoj) {
            residualj = psij[, i] - psij[, Sjbu] %*% keaij[, sj1]
            resj = sum(residualj^2)
            BIC2j[sj1] = log(resj) + length(which(keaij[, sj1] != 0)) * log(p) *
              log(log(n)) / n
            if (BIC2j[sj1] <= BIC2j[min1]) {
              min1 = sj1
            }
          }
          
          zj[, i] = psij[, i] - psij[, Sjbu] %*% keaij[, min1]
        }
        else if (i %in% Sjbu) {
          psij[, Sjbu] = (I - Pj) %*% Xjbu[, Sjbu]
          fitj = glmnet(psij[, setdiff(Sjbu, i)], psij[, i], intercept = 0)
          keaij = coef(fitj)[-1, ]
          meihaoj = ncol(keaij)
          BIC2j = c(rep(0, meihaoj))
          min1 = 1
          for (sj1 in 1:meihaoj) {
            residualj = psij[, i] - psij[, setdiff(Sjbu, i)] %*% keaij[, sj1]
            resj = sum(residualj^2)
            BIC2j[sj1] = log(resj) + length(which(keaij[, sj1] != 0)) * log(p) *
              log(log(n)) / n
            if (BIC2j[sj1] <= BIC2j[min1]) {
              min1 = sj1
            }
          }
          
          zj[, i] = psij[, i] - psij[, setdiff(Sjbu, i)] %*% keaij[, min1]
        }
        
        zetaj[i] = ((t(zj[, i]) %*% zj[, i])^0.5) / abs(t(Xjbu[, i]) %*% zj[, i])
        
        gamma_hatj[i] = (t(zj[, i]) %*% Xj) / (t(Xjbu[, i]) %*% zj[, i])
        
        p_valuj[i] = 2 * (1 - pnorm(abs(gamma_hatj[i]) / (zetaj[i] * taus)))
        
      }
      gamma_hat[, j] = t(gamma_hatj)
      p_valu[, j] = t(p_valuj)
      
    }
    return(p_valu)
  }



#' The false positive rates and the true positive rates based on the graphical projection estimator (GPE)
#' @title simple function
#' @param M M is the number of parallelism
#' @param p p is the dimension of graphical model
#' @param resultmcp resultmcp is a list of matrices generated based on the results of function parallelism
#' @return A list of 
#' \itemize{
#'  \item{\code{EPs}: the true positive rates of Hypothesis testing}
#'  \item{ETs: the false positive rates of Hypothesis testing}
#' }
#' @export
tes.fgpe <-
  function(M, p, resultmcp) {
    for (t in 1:M) {
      for (jj in 1:p) {
        a[, jj, t] = resultmcp [[t]][, jj]
      }
    }
    s1 = matrix(0, p - 1, p)
    for (jj in 1:p) {
      for (ii in 1:(p - 1)) {
        s1[ii, jj] = sum(I(a[ii, jj, ] < 0.05))
      }
    }
    s = matrix(0, p, p)
    for (jj in 1:p) {
      s[-jj, jj] = s1[, jj]
    }
    EP = matrix(0, p, p)
    for (ii in 1:p) {
      for (jj in 1:p) {
        if (abs(ii - jj) == 1) {
          EP[ii, jj] = s[ii, jj]
        }
        if (abs(ii - jj) == 2) {
          EP[ii, jj] = s[ii, jj]
        }
      }
    }
    
    ET = s - EP
    
    ETs = sum(ET) / ((p * p - p - 4 * p + 6) * M)
    
    EPs = sum(EP) / ((p * 4 - 6) * M)
    return(list(ETs = ETs, EPs = EPs))
  }