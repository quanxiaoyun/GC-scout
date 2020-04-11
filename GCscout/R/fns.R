#' The rank-based estimate of the correlation matrix
#'
#' \code{cov_mat} computes the estimate of the correlation matrix between the columns of \code{dat}.
#' @param dat a numeric matrix or data frame, of dimensions n by p. If matrix, each row is a length-p vector representing one observation on p variables.
#' @param ordinal_index a numeric or logical vector. The indices of columns in \code{dat} that are of ordinal values (rather than continuous).
#' @examples
#' dat <- cbind(c(0,1,0,0,1), c(1,0,1,1,0), rnorm(5), rnorm(5))
#' cov_mat(dat, c(1,2))
#' @importFrom pbivnorm pbivnorm
#' @importFrom mvtnorm pmvnorm
#' @export
cov_mat = function(dat, ordinal_index){
  n = nrow(dat)
  p = ncol(dat)
  M = diag(p)
  cts_ind = (1:p)[-ordinal_index]
  M[ordinal_index, ordinal_index] = ordinal_block_cov(dat[,ordinal_index])
  M[ordinal_index, cts_ind] = ord_cts_cov(dat[,ordinal_index], dat[,cts_ind])
  M[cts_ind, ordinal_index] = M[ordinal_index, cts_ind]
  M[cts_ind, cts_ind] = cts_block_cov(dat[,cts_ind])
  return(M)
}


#' Expectation of sample Kendall's tau.
#'
#' @param r Correlation coefficient between the two latent variables.
#' @param dj1 The cutoff \eqn{\Delta_j^1}.
#' @param dj2 The cutoff \eqn{\Delta_j^2}.
#' @param dk1 The cutoff \eqn{\Delta_k^1}.
#' @param dk2 The cutoff \eqn{\Delta_k^2}.
#' @return The population version of Kendall's tau.
#' @examples
#' G_ordinal(.5, -0.3, 0.3, 0, Inf)
#' @importFrom pbivnorm pbivnorm
#' @importFrom stats pnorm
#' @export
G_ordinal <- function(r,dj1,dj2,dk1,dk2){
  if (is.infinite(dj2) && is.infinite(dk2)){
    # both variables are of binary data
    (pbivnorm(dj1,dk1,r)-pnorm(dj1)*pnorm(dk1))*2
  }else if (is.infinite(dj2)){
    (pnorm(dk2)*(1-pnorm(dk1)-pnorm(dj1)+pbivnorm(dj1,dk1,r))-(pnorm(dk1) - 1)*(pbivnorm(dj1,dk2,r)- pnorm(dk2)))*2
  }else if (is.infinite(dk2)){
    (pnorm(dj2)*(1-pnorm(dk1)-pnorm(dj1)+pbivnorm(dj1,dk1,r))-(pbivnorm(dj2,dk1,r) - pnorm(dj2))*(pnorm(dj1)- 1))*2
  }else{
    (pbivnorm(dj2,dk2,r)*(1-pnorm(dk1)-pnorm(dj1)+pbivnorm(dj1,dk1,r))-(pbivnorm(dj2,dk1,r) - pnorm(dj2))*(pbivnorm(dj1,dk2,r)- pnorm(dk2)))*2
  }
}

#' Moment estimate of the cutoffs for discretization.
#'
#' @param x a numeric vector. Observed values of a variable. Must be of (ordinal) discrete values.
#' @return The moment estimate of the cutoffs.
#' @examples
#' x <- c(0, 1, 0, 0, 1)
#' d_hat_vec(x)
#' @importFrom stats qnorm
#' @export
d_hat_vec <- function(x){
  if(length(levels(as.factor(x))) == 2){
    c(qnorm(sum(x==min(x))/length(x)), Inf)
  }else{
    c(qnorm(sum(x==min(x))/length(x)), qnorm(1-sum(x==max(x))/length(x)))
  }
}



#' @importFrom stats uniroot
inverse = function(f){
  function(y) uniroot((function(x) f(x)-y),c(-1,1))$root
}

#' Kendall's tau-a.
#'
#' Compute the Kendall rank correlation coefficient (Tau-a statistic) between two variables.
#' @param x A numeric matrix.
#' @param y A numeric matrix.
#' @param j an integer. Indicating which column of \code{x} is to be considered for this function.
#' @param k an integer. Indicating which column of \code{y} is to be considered.
#' @param n an integer. The number of observations. Has to be same for \code{x} and \code{y}.
#' @return The Kendall's tau-a correlation coefficient between the \code{j}-th column of \code{x} and the \code{k}-th column of \code{y}.
#' @examples
#' x <- cbind(c(0,1,0,0,1), c(1,1,0,1,0))
#' y <- cbind(rnorm(5), rnorm(5))
#' tau_hat(x, x, 1, 2, 5)
#' tau_hat(x, y, 1, 1, 5)
#' @export
tau_hat = function(x,y,j,k,n){
  # x,y are both binary or both continuous
  # this function compute the tau between x[,j] and y[,k]
  tau = 0;
  for(i in 1:n){
    for(i_prime in i:n){
      tau = tau + sign((x[i,j]- x[i_prime,j])*
                         (y[i,k] - y[i_prime,k]))
    }
  }
  return(tau/(n*(n-1)/2))
}


ordinal_block_cov = function(x){
  x = as.matrix(x)
  n = nrow(x)
  a = ncol(x)
  M = diag(a)

  if(a>1){
    for(k in 2:a){
      dkhat = d_hat_vec(x[,k])
      for(j in 1:(k-1)){
        djhat = d_hat_vec(x[,j])
        tau = tau_hat(x,x,j,k,n)
        dj1 = djhat[1]
        dj2 = djhat[2]
        dk1 = dkhat[1]
        dk2 = dkhat[2]
        if(tau > G_ordinal(1,dj1,dj2,dk1,dk2)){
          M[j,k] = 1
          M[k,j] = 1
        }else if(tau < G_ordinal(-1,dj1,dj2,dk1,dk2)){
          M[j,k] = -1
          M[k,j] = -1
        }else{
          G_inverse = inverse(function(z) G_ordinal(z,dj1,dj2,dk1,dk2))
          r = G_inverse(tau)
          M[j,k] = r
          M[k,j] = r
        }}}}
  return(M)
}

#' @importFrom stats cor
cts_block_cov = function(x){
  x = as.matrix(x)
  n = nrow(x)
  b = ncol(x)
  M = diag(b);
  tau_mat = cor(x, method = 'kendall')
  if(b>1){
    for(k in 2:b){
      for(j in 1:(k-1)){
        tau = tau_mat[j,k]
        r = sin(pi*tau/2)
        M[j,k] = r
        M[k,j] = r
      }
    }
  }
  return(M)
}

ord_cts_cov = function(x,y){
  x = as.matrix(x) # ordinal
  y = as.matrix(y) # continuous
  a = ncol(x); b = ncol(y); n = nrow(x)
  M = matrix(0,nrow=a,ncol=b)
  for(j in 1:a){
    djhat = d_hat_vec(x[,j])
    for(k in 1:b){
      tau = tau_hat(x,y,j,k,n)
      if(tau < G_mixed_k(-1, djhat)){
        M[j,k] = -1
      }else if(tau > G_mixed_k(1, djhat)){
        M[j,k] = 1
      }else{
        G_inverse = inverse(function(z) G_mixed_k(z, djhat))
        a = G_inverse(tau)
        M[j,k] = a
      }}}
  return(M)
}

G_mixed_k = function(r, dj){
  sig = diag(3)
  sig[1,3] = sig[3,1] = r/sqrt(2)
  sig[2,3] = sig[3,2] = -r/sqrt(2)
  p = length(dj)
  sum = 4*pbivnorm(dj[p], 0, r/sqrt(2)) - 2* pnorm(dj[p])
  for (i in 1: (p-1)){
    upr1 = c(dj[i], dj[i+1], 0)
    upr2 = c(dj[i + 1], dj[i], 0)

    sum = sum + 2*pmvnorm(lower=-Inf,upper = upr1,sigma=sig)[1] -
      2*pmvnorm(lower=-Inf,upper = upr2,sigma=sig)[1]
  }
  return(sum)
}

#' Perform cross-validation for the GC-Scout regression
#'
#' This function returns cross-validation errors for a range of tuning parameter
#' values (lam1 and lam2). Error measures can be specified as either MAPE (mean absolute
#' percentage error) or MSE (Mean Squared Error). Note that this function allows GC_Scout(0,1),
#' GC_Scout(1,1), and GC_Scout(2,1) models to be either fit individually or altogether in one run.
#' @param K An integer. Number of cross-validation folds to be performed; default is 5.
#' @param MAPE Logical. If TRUE, MAPE (mean absolute percentage error) will be
#' used to measure cross-validation error; if FALSE, MSE (Mean Squared Error)
#' will be used instead. Default is FALSE.
#' @inheritParams GC_scout
#' @return A list with the following elements:
#' \itemize{
#'   \item mse_gc_scout_01 - the beta estimates for GC_Scout(0,1) model (only if \code{gc_scout_01 = TRUE}; NULL otherwise)
#'   \item mse_gc_scout_11 - the beta estimates for GC_Scout(2,1) model (only if \code{gc_scout_21 = TRUE}; NULL otherwise)
#'   \item mse_gc_scout_21 - the beta estimates for GC_Scout(1,1) model (only if \code{gc_scout_11 = TRUE}; NULL otherwise)
#'   \item bestlam1_11 - the best \eqn{\lamdba_1} associated with lowest cv error (only if \code{gc_scout_11 = TRUE}; NULL otherwise)
#'   \item bestlam1_21 - the best \eqn{\lamdba_1} associated with lowest cv error (only if \code{gc_scout_21 = TRUE}; NULL otherwise)
#'   \item bestlam2_01 - the best \eqn{\lamdba_2} associated with lowest cv error (only if \code{gc_scout_01 = TRUE}; NULL otherwise)
#'   \item bestlam2_11 - the best \eqn{\lamdba_2} associated with lowest cv error (only if \code{gc_scout_11 = TRUE}; NULL otherwise)
#'   \item bestlam2_21 - the best \eqn{\lamdba_2} associated with lowest cv error (only if \code{gc_scout_21 = TRUE}; NULL otherwise)
#' }
#' @examples
#' data(Sigma_witten)
#' n <- 40; d <- ncol(Sigma_witten); d3 <- 4;
#' delta = c(qnorm(0.5), qnorm(0.75))
#' dat <- mvtnorm::rmvnorm(n, sigma = Sigma_witten)
#' y_tilde = dat[,d]
#' y = y_tilde; y = exp(y_tilde);
#' x_tilde = dat[,-d]; x = x_tilde^5
#' x[,1:d3] = as.numeric(cut(x[,1:d3], breaks = c(-Inf,(delta)^5,Inf)))-1
#' cv.GC_scout(x, y, lam1s = seq(0.01,0.2,length.out = 10), K=3,
#'  lam2s = seq(0.01,0.2,length.out = 10), d2 = seq(1:4),
#'  gc_scout_01 = FALSE, gc_scout_21 = TRUE, gc_scout_11 = FALSE)

#' @seealso \code{\link{GC_scout}}
#' @export
cv.GC_scout = function(x,y, d2, lam1s = seq(0.01,0.25, length.out = 10), lam2s = seq(0.01,0.25,length.out = 10),
                       K=5, gc_scout_01 = TRUE, gc_scout_21 = TRUE,
                       gc_scout_11 = TRUE, MAPE=FALSE){

  all.folds = cv.folds(length(y), K)

  cv_error_gc_scout_01 = matrix(0, length(lam2s), K)
  cv_error_gc_scout_21 = array(0, dim=c(length(lam1s), length(lam2s), K))
  cv_error_gc_scout_11 = array(0, dim=c(length(lam1s), length(lam2s), K))

  for (k in seq(K)){
    vldn = all.folds[[k]]

    # get betahat for each folds
    cv.result = GC_scout(x[-vldn,], y[-vldn], d2=d2, lam1s, lam2s,
                         gc_scout_01 = gc_scout_01,
                         gc_scout_21 = gc_scout_21,
                         gc_scout_11 = gc_scout_11)

    if (gc_scout_01 == TRUE){
      # beta_gc_scout_01
      beta_gc_scout_01_mat = cv.result$beta_gc_scout_01_mat
      cv_error_gc_scout_01[,k] = pred_error(x_new = x[vldn, ], y_new = y[vldn],
                                             x = x[-vldn, ], y = y[-vldn], beta_gc_scout_01_mat, MAPE)
    }

    if (gc_scout_21){
      # beta_gc_scout_21
      beta_scout_21_mat = cv.result$beta_gc_scout_21_mat
      cv_error_gc_scout_21[,,k] = pred_error(x[vldn, ], y[vldn], x[-vldn, ], y[-vldn],beta_scout_21_mat, MAPE)
    }

    if (gc_scout_11){
      # beta_gc_scout_11
      beta_scout_11_mat = cv.result$beta_gc_scout_11_mat
      cv_error_gc_scout_11[,,k] = pred_error(x[vldn, ], y[vldn], x[-vldn, ], y[-vldn],beta_scout_11_mat, MAPE)
    }

  }

  if (gc_scout_21){
    cv_error_21 = apply(cv_error_gc_scout_21, c(1,2), mean)
    bestlam1_21 = lam1s[which.min(rowMeans(cv_error_21))]
    bestlam2_21 = lam2s[which.min(colMeans(cv_error_21))]
  }else{
    bestlam1_21 = NULL
    bestlam2_21 = NULL

  }
  if (gc_scout_11==TRUE){
    cv_error_11 = apply(cv_error_gc_scout_11, c(1,2), mean)
    bestlam1_11 = lam1s[which.min(rowMeans(cv_error_11))]
    bestlam2_11 = lam2s[which.min(colMeans(cv_error_11))]
  }else{
    bestlam1_11 = NULL
    bestlam2_11 = NULL
  }
  if (gc_scout_01 == TRUE){
    cv_error_01 = apply(cv_error_gc_scout_01, 1, mean)
    bestlam2_01 = lam2s[which.min(cv_error_01)]
  }else{
    bestlam2_01 = NULL
  }
  return(list(cv_error_gc_scout_01 = cv_error_gc_scout_01, cv_error_gc_scout_21 = cv_error_gc_scout_21,
              cv_error_gc_scout_11 = cv_error_gc_scout_11,
              bestlam1_11 = bestlam1_11, bestlam2_11 = bestlam2_11,
              bestlam1_21 = bestlam1_21, bestlam2_21 = bestlam2_21,
              bestlam2_01 = bestlam2_01))
}

#' Fit the covariance-regularized regression under Gaussian
#' copula model (GC-Scout).
#'
#' This function implements Algorithm 1 from Quan et al (2019).
#' Specifically, for observed pairs (\code{x},\code{y}), it estimates
#' the regression coefficients \eqn{\beta} for the latent linear model under
#' latent Gaussian Copula model: \eqn{\tilde{y} = \tilde{x}\beta + \epsilon}
#' @param x A numeric matrix. Used for training.
#' @param y A numeric vector. Used for training.
#' @param lam1s A numeric number or vector. The tuning parameters for
#' regularization of the covariance matrix. If
#' \code{gc_scout_01 = TRUE}, then no covariance regularization
#' is taking place and \code{lam1s} will be ignored. If
#' \code{gc_scout_11 = TRUE},  as graphical lasso can be uncomfortably slow.
#' @param lam2s A numeric number or vector. The tuning parameters
#' for the \eqn{L_1} regularization of the regression coefficients,
#' using the regularized covariance matrix.
#' @param d2 A numeric vector. The indices of discrete variables in \code{x}.
#' @param gc_scout_01 A logical value. If TRUE, GC_Scout(0,1) (aka GC Lasso) will be performed. Default is TRUE.
#' @param gc_scout_21 A logical value. If TRUE, GC_Scout(2,1) will be performed. Default is TRUE.
#' @param gc_scout_11 A logical value. If TRUE, GC_Scout(1,1) will be performed. Default is TRUE.
#' @return
#' \itemize{
#'   \item beta_gc_scout_01_mat - the beta estimates for GC_Scout(0,1) model (\code{gc_scout_01 = TRUE})
#'   \item beta_gc_scout_21_mat - the beta estimates for GC_Scout(2,1) model (\code{gc_scout_21 = TRUE})
#'   \item beta_gc_scout_11_mat - the beta estimates for GC_Scout(1,1) model (\code{gc_scout_11 = TRUE})
#' }
#' @examples
#' data(Sigma_witten)
#' n <- 40; d <- ncol(Sigma_witten); d3 <- 4;
#' delta = c(qnorm(0.5), qnorm(0.75))
#' dat <- mvtnorm::rmvnorm(n, sigma = Sigma_witten)
#' y_tilde = dat[,d]
#' y = y_tilde; y = exp(y_tilde);
#' x_tilde = dat[,-d]; x = x_tilde^5
#' x[,1:d3] = as.numeric(cut(x[,1:d3], breaks = c(-Inf,(delta)^5,Inf)))-1
#' GC_scout(x, y, lam1s = seq(0.01,0.2,length.out = 10),
#'  lam2s = seq(0.01,0.2,length.out = 10), d2 = seq(1:4),
#'  gc_scout_01 = FALSE, gc_scout_21 = TRUE, gc_scout_11 = FALSE)
#' @importFrom scout crossProdLasso
#' @importFrom glasso glasso
#' @importFrom Matrix nearPD
#' @export
GC_scout = function(x, y, lam1s, lam2s, d2,
                    gc_scout_01 = TRUE, gc_scout_21 = TRUE, gc_scout_11 = TRUE){
  # estimate the correlation matrix for latent variabbles (x,y)
  cov_m = cov_mat(cbind(x,y), d2)
  d = ncol(x) + 1
  cov_xy = cov_m[-d,d]
  cov_x = cov_m[-d,-d]

  beta_gc_scout_01_mat = matrix(0, ncol(x), length(lam2s))
  beta_gc_scout_21_mat = array(0, dim = c(ncol(x), length(lam1s), length(lam2s)))
  beta_gc_scout_11_mat = array(0, dim = c(ncol(x), length(lam1s), length(lam2s)))

  for (j in 1:length(lam2s)){
    if(gc_scout_01 == TRUE){
      beta_gc_scout_01 = crossProdLasso(cov_x, cov_xy, rho = lam2s[j])$beta
      beta_gc_scout_01_mat[,j] = as.vector(beta_gc_scout_01)
    }

    for (i in 1:length(lam1s)){
      if (gc_scout_21 == TRUE){
        g.out = gc_ridge(cov_x, rho = lam1s[i])
        beta_gc_scout_21 = lasso_one(g.out, cov_xy, rho = lam2s[j])$beta
        beta_gc_scout_21_mat[,i,j] = as.vector(beta_gc_scout_21)
      }

      if (gc_scout_11 == TRUE){
        g.out <- glasso(as.matrix(nearPD(cov_x)$mat), rho=lam1s[i])
        beta_gc_scout_11 <- lasso_one(g.out$w, cov_xy, rho=lam2s[j])$beta
        beta_gc_scout_11_mat[,i,j] = as.vector(beta_gc_scout_11)
      }
    }
  }
  return(list(beta_gc_scout_01_mat = beta_gc_scout_01_mat,
              beta_gc_scout_21_mat = beta_gc_scout_21_mat,
              beta_gc_scout_11_mat = beta_gc_scout_11_mat))
}


gc_ridge = function(cov_x, rho){
  e = eigen(cov_x)
  d_tilde= -sqrt(2*rho) + (-e$values + sqrt(e$values^2+8*rho))/2
  cov_x_ridge = sqrt(2*rho) * diag(ncol(cov_x)) + e$vectors %*% diag(e$values + d_tilde) %*% t(e$vectors)
  return(cov_x_ridge)
}

cv.folds <- function(n, folds = 10){
  split(sample(1:n), rep(1:folds, length = n))
}

cut_x = function(x, delta){
  y = x
  if(length(delta) > 1 &&  delta[1] == delta[2]){
    delta[1] = as.numeric(names(table(x))[1])
    delta[2] = as.numeric(names(table(x))[2])
  }
  k = length(delta)
  delta = c(-Inf, delta, Inf)
  for (i in 1:(k + 1)){
    y[which(x <= delta[i + 1]  & x > delta[i])] = i-1
  }
  return(y)
}



#' Predict y for given x using GC-Scout method.
#'
#' @param x_new A numeric matrix. Must have the same number of columns as \code{x}
#' @param x A numeric matrix. Used for training.
#' @param y A numeric vector. Used for training.
#' @param betahat A numeric matrix or vector. If vector, its length must be equal to \code{x}'s number of columns. If matrix, each column represents one beta.
#' @return Given the observed pairs (\code{x},\code{y}), and the
#' \eqn{\beta} estimate \code{betahat}, one can predict new response
#' values for new covariates \code{x_new}, under latent Gaussian Copula model.
#' Specifically, this function implements Algorithm 2 from Quan et al (2019).
#' @examples
#' n <- 200; d <- ncol(Sigma_witten); d3 <- 4; delta = c(qnorm(0.5), qnorm(0.75))
#' dat <- mvtnorm::rmvnorm(n, sigma = Sigma_witten)
#' y_tilde = dat[,d]
#' y = y_tilde; y = exp(y_tilde);
#' x_tilde = dat[,-d]; x = x_tilde^5
#' x[,1:d3] = as.numeric(cut(x[,1:d3], breaks = c(-Inf,(delta)^5,Inf)))-1
#' x_train = x[1:40,]; y_train = y[1:40]
#' x_test = x[41:n,]; y_test = y[41:n]
#' predict(x_test, x_train, y_train, c(0.5,0.3,0,0,0,0,-0.2,-0.2))
#' @export
predict = function(x_new, x, y, betahat){
  # compute x_tilde
  x_tilde = x
  for(i in 1:nrow(x)){
    x_tilde[i,] = unlist(lapply(1:ncol(x), function(j) f_hat(x[i,j], x[,j])))
  }

  # compute y_tilde
  y_tilde = unlist(lapply(1:length(y), function(j) f_hat(y[j], y)))

  # compute x_new_tilde
  x_new_tilde = x_new
  for(i in 1:nrow(x_new_tilde)){
    x_new_tilde[i,] = unlist(lapply(1:ncol(x), function(j) f_hat(x_new_tilde[i,j], x[,j])))
  }

  # compute intercept first then y_new_tilde_hat
  if(length(dim(betahat)) == 2){
    interceptmat = apply(betahat, 2, function(beta) as.numeric(mean(y_tilde) - apply(x_tilde,2, mean) %*% beta))
    y_new_tilde = sweep(apply(betahat, 2, function(beta) x_new_tilde %*% beta),1, interceptmat, '+')
    y_new_hat = apply(y_new_tilde, 2, f_0_inv, y)
  }else if(length(dim(betahat)) == 3){
    interceptmat = apply(betahat, c(2,3), function(beta) as.numeric(mean(y_tilde) - apply(x_tilde,2, mean) %*% beta))
    y_new_tilde = sweep(apply(betahat, c(2,3), function(beta) x_new_tilde %*% beta), c(2,3), interceptmat, "+")
    y_new_hat = apply(y_new_tilde, c(2,3), f_0_inv, y)
  }else{
    interceptmat = as.numeric(mean(y_tilde) - apply(x_tilde,2, mean) %*% betahat)
    y_new_tilde = as.matrix(x_new_tilde) %*% as.matrix(betahat) + interceptmat
    y_new_hat = f_0_inv(y_new_tilde, y)
  }

  return(y_new_hat)
}

#' @importFrom stats ecdf
f_hat = function(t, y){
  # compute f_hat(x) using winsorized F_tilde -- winsorized empirical cdf
  n = length(y)
  F_tilde = ecdf(y)(t)
  if(F_tilde < 1/n^2){
    F_hat = 1/n^2
  }else if (F_tilde <= 1-(1/n)^2){
    F_hat = F_tilde
  }else{
    F_hat = 1- (1/n)^2
  }
  return(qnorm(F_hat))
}



#' @importFrom stats quantile
f_0_inv = function(t,y){

  phi_t = pnorm(t)
  n = length(y)
  # winsorized empirical CDF
  phi_t = phi_t *(phi_t > 1/n^2 & phi_t < 1 - 1/n^2) + 1/n^2 *(phi_t < 1/n^2) +  (1-1/n^2)* (phi_t > 1 - 1/n^2)
  return(as.numeric(quantile(y, phi_t, type = 1)))
}

lasso_one=function(w,ss, rho, thr=1.0e-4, maxit=100,trace=F, beta.init=NULL){
  # does lasso fit of a single response variable  on predictors,
  # via coordinate descent. Data is inner products w and ss

  n=length(ss)

  if(length(rho)==1){rho=rep(rho,n)}

  itrace=1*trace

  if(is.null(beta.init)){ beta.init=rep(0,n)}
  mode(rho)="single"
  mode(ss)="single"
  mode(w)="single"
  mode(n)="integer"
  mode(maxit)="integer"
  mode(itrace)="integer"
  mode(thr)="single"
  mode(beta.init)="single"

  junk<-.Fortran("lasso7", rho, n, as.matrix(w), ss, thr, xx=beta.init, PACKAGE="scout")

  return(list(beta=junk$xx))
}


#' Compute prediction errors of GC-scout method.
#'
#' @param y_new A numeric vector. The length must be equal to the number of rows of \code{x}.
#' @param MAPE Logical. If TRUE, MAPE (mean absolute percentage error) will be
#' used to measure cross-validation error; if FALSE, MSE (Mean Squared Error)
#' will be used instead. Default is FALSE.
#' @inheritParams predict
#' @return The prediction errors following the Algorithm 2 in Quan et al (2020),
#' that can be either Mean Absolute Percentage Error (MAPE) if \code{MAPE = TRUE}
#' or Mean Squared Error if \code{MAPE = FALSE}.
#'
#' @seealso \code{\link{predict}}
#' @examples
#' n <- 200; d <- ncol(Sigma_witten); d3 <- 4; delta = c(qnorm(0.5), qnorm(0.75))
#' dat <- mvtnorm::rmvnorm(n, sigma = Sigma_witten)
#' y_tilde = dat[,d]
#' y = y_tilde; y = exp(y_tilde);
#' x_tilde = dat[,-d]; x = x_tilde^5
#' x[,1:d3] = as.numeric(cut(x[,1:d3], breaks = c(-Inf,(delta)^5,Inf)))-1
#' x_train = x[1:40,]; y_train = y[1:40]
#' x_test = x[41:n,]; y_test = y[41:n]
#' pred_error(x_test,y_test, x_train, y_train, c(0.5,0.3,0,0,0,0,-0.2,-0.2),MAPE = TRUE)
#' @export
pred_error = function(x_new, y_new, x, y, betahat, MAPE=FALSE){
  # betahat is matrix of beta's (Each column is one beta)

  yhats = predict(x_new,x,y,betahat)

  if(is.vector(betahat)){
    if(MAPE){
      error = mean(abs(yhats - y_new)/y_new)
    }else{
      error = mean((yhats - y_new)^2)
    }

  }else{
    if(MAPE){
      error = colMeans(abs(yhats-y_new)/y_new)
    }else{
      error = colMeans((yhats-y_new)^2)
    }
  }
  return(error)
}

