#' Internal Matrix Helper Functions
#' @noRd
cormax.ind <- function(n) {
  diag(n)
}

cormax.exch <- function(n, alpha) {
  mat <- matrix(alpha, n, n)
  diag(mat) <- 1
  mat
}

cormax.ar1 <- function(n, alpha) {
  exponent <- abs(matrix(1:n - 1, nrow = n, ncol = n, byrow = TRUE) -
                    matrix(1:n - 1, nrow = n, ncol = n, byrow = FALSE))
  alpha^exponent
}

#' Kauermann & Carroll (KC) Variance Estimator
#' @noRd
GEE.var.kc.r <- function(formula, id, family = gaussian, data, corstr = "independence") {
  environment(formula) <- environment()

  if (!id %in% names(data)) stop("ID variable not found.")
  data <- data[order(data[[id]]), ]
  id_vec <- data[[id]]

  invisible(utils::capture.output({
    suppressMessages({
      gee.fit <- gee::gee(formula, data = data, id = id_vec, family = family, corstr = corstr, silent = TRUE)
    })
  }))

  beta_est <- gee.fit$coefficient
  alpha <- gee.fit$working.correlation[1, 2]
  len <- length(beta_est)

  unique_ids <- unique(id_vec)
  cluster_n <- as.vector(table(id_vec)[as.character(unique_ids)])
  ncluster <- max(cluster_n)
  size <- length(unique_ids)

  m <- model.frame(formula, data)
  data$response <- model.response(m, "numeric")
  mat <- as.data.frame(model.matrix(formula, m))
  mat$subj <- rep(unique_ids, cluster_n)

  if (is.character(corstr)) {
    var <- switch(corstr,
                  "independence" = cormax.ind(ncluster),
                  "exchangeable" = cormax.exch(ncluster, alpha),
                  "AR-M" = cormax.ar1(ncluster, alpha),
                  "unstructured" = summary(gee.fit)$working.correlation
    )
  }

  if (is.character(family)) {
    family <- switch(family, "gaussian" = "gaussian", "binomial" = "binomial", "poisson" = "poisson")
  } else if (is.function(family)) {
    family <- family()[[1]]
  }

  step11 <- matrix(0, nrow = len, ncol = len)
  for (i in 1:size) {
    covariate <- as.matrix(subset(mat[, -length(mat[1, ])], mat$subj == unique_ids[i]))
    var_i <- var[1:cluster_n[i], 1:cluster_n[i]]
    if (family == "poisson") {
      D <- geesmv::mat.prod(covariate, exp(covariate %*% beta_est))

      Vi <- diag(sqrt(c(exp(covariate %*% beta_est))), cluster_n[i]) %*% var_i %*% diag(sqrt(c(exp(covariate %*% beta_est))), cluster_n[i])
      xx <- t(D) %*% solve(Vi) %*% D
      step11 <- step11 + xx
    } else if (family == "gaussian") {
      xx <- t(covariate) %*% solve(var_i) %*% covariate
      step11 <- step11 + xx
    }
  }

  step12 <- matrix(0, nrow = len, ncol = len)
  step13 <- matrix(0, nrow = len^2, ncol = 1)
  step14 <- matrix(0, nrow = len^2, ncol = len^2)
  p <- matrix(0, nrow = len^2, ncol = size)

  for (i in 1:size) {
    y <- as.matrix(data$response[id_vec == unique_ids[i]])
    covariate <- as.matrix(subset(mat[, -length(mat[1, ])], mat$subj == unique_ids[i]))
    var_i <- var[1:cluster_n[i], 1:cluster_n[i]]

    if (family == "poisson") {
      D <- geesmv::mat.prod(covariate, exp(covariate %*% beta_est))
      Vi <- diag(sqrt(c(exp(covariate %*% beta_est))), cluster_n[i]) %*% var_i %*% diag(sqrt(c(exp(covariate %*% beta_est))), cluster_n[i])

      resid_term <- cormax.ind(cluster_n[i]) - D %*% solve(step11) %*% t(D) %*% solve(Vi)
      inv_term <- geesmv::mat.sqrt.inv(resid_term)

      xy <- t(D) %*% solve(Vi) %*% inv_term %*% (y - exp(covariate %*% beta_est))
      step12 <- step12 + xy %*% t(xy)
      step13 <- step13 + matrixcalc::vec(xy %*% t(xy))
      p[, i] <- matrixcalc::vec(xy %*% t(xy))
    } else if (family == "gaussian") {
      resid_term <- cormax.ind(cluster_n[i]) - covariate %*% solve(step11) %*% t(covariate) %*% solve(var_i)
      inv_term <- geesmv::mat.sqrt.inv(resid_term)

      xy <- t(covariate) %*% solve(var_i) %*% inv_term %*% (y - covariate %*% beta_est)
      step12 <- step12 + xy %*% t(xy)
      step13 <- step13 + matrixcalc::vec(xy %*% t(xy))
      p[, i] <- matrixcalc::vec(xy %*% t(xy))
    }
  }
  for (i in 1:size) {
    dif <- (p[, i] - step13 / size) %*% t(p[, i] - step13 / size)
    step14 <- step14 + dif
  }
  cov.beta <- solve(step11) %*% (step12) %*% solve(step11)
  return(list(cov.beta = cov.beta))
}

#' Mancl & DeRouen (MD) Variance Estimator
#' @noRd
GEE.var.md.r <- function(formula, id, family = gaussian, data, corstr = "independence") {
  environment(formula) <- environment()

  if (!id %in% names(data)) stop("ID variable not found.")
  data <- data[order(data[[id]]), ]
  id_vec <- data[[id]]

  invisible(utils::capture.output({
    suppressMessages({
      gee.fit <- gee::gee(formula, data = data, id = id_vec, family = family, corstr = corstr, silent = TRUE)
    })
  }))

  beta_est <- gee.fit$coefficient; alpha <- gee.fit$working.correlation[1,2]; len <- length(beta_est)
  unique_ids <- unique(id_vec)
  cluster_n <- as.vector(table(id_vec)[as.character(unique_ids)])
  ncluster <- max(cluster_n)
  size <- length(unique_ids)

  m <- model.frame(formula, data); data$response <- model.response(m, "numeric"); mat <- as.data.frame(model.matrix(formula, m))
  mat$subj <- rep(unique_ids, cluster_n)

  if(is.character(corstr)){
    var <- switch(corstr,
                  "independence"=cormax.ind(ncluster),
                  "exchangeable"=cormax.exch(ncluster, alpha),
                  "AR-M"=cormax.ar1(ncluster, alpha),
                  "unstructured"=summary(gee.fit)$working.correlation)
  }
  if(is.character(family)){ family <- switch(family, "gaussian"="gaussian", "binomial"="binomial", "poisson"="poisson") } else if(is.function(family)){ family <- family()[[1]] }

  step11<-matrix(0,nrow=len,ncol=len)
  for (i in 1:size){
    covariate<-as.matrix(subset(mat[,-length(mat[1,])], mat$subj==unique(id_vec)[i]))
    var_i=var[1:cluster_n[i],1:cluster_n[i]]
    if (family=="poisson") {
      D <- geesmv::mat.prod(covariate, exp(covariate %*% beta_est))

      Vi <- diag(sqrt(c(exp(covariate %*% beta_est))),cluster_n[i])%*%var_i%*%diag(sqrt(c(exp(covariate%*%beta_est))),cluster_n[i])
      xx<-t(D)%*%solve(Vi)%*%D
      step11<-step11+xx
    } else if (family=="gaussian"){ xx<-t(covariate)%*%solve(var_i)%*%covariate; step11<-step11+xx }
  }

  step12<-matrix(0,nrow=len,ncol=len); step13<-matrix(0,nrow=len^2,ncol=1); step14<-matrix(0,nrow=len^2,ncol=len^2); p<-matrix(0,nrow=len^2,ncol=size)
  for (i in 1:size){
    y<-as.matrix(data$response[id_vec==unique(id_vec)[i]])
    covariate<-as.matrix(subset(mat[,-length(mat[1,])], mat$subj==unique(id_vec)[i]))
    var_i=var[1:cluster_n[i],1:cluster_n[i]]
    if (family=="poisson") {
      D <- geesmv::mat.prod(covariate, exp(covariate %*% beta_est))

      Vi <- diag(sqrt(c(exp(covariate %*% beta_est))),cluster_n[i])%*%var_i%*%diag(sqrt(c(exp(covariate%*%beta_est))),cluster_n[i])

      xy<-t(D)%*%solve(Vi)%*%solve(cormax.ind(cluster_n[i])-D%*%solve(step11)%*%t(D)%*%solve(Vi))%*%(y-exp(covariate%*%beta_est))
      step12<-step12+xy%*%t(xy); step13<-step13+matrixcalc::vec(xy%*%t(xy)); p[,i]<-matrixcalc::vec(xy%*%t(xy))
    } else if (family=="gaussian"){
      xy<-t(covariate)%*%solve(var_i)%*%solve(cormax.ind(cluster_n[i])-covariate%*%solve(step11)%*%t(covariate)%*%solve(var_i))%*%(y-covariate%*%beta_est)
      step12<-step12+xy%*%t(xy); step13<-step13+matrixcalc::vec(xy%*%t(xy)); p[,i]<-matrixcalc::vec(xy%*%t(xy))
    }
  }
  for (i in 1:size){ dif<-(p[,i]-step13/size)%*%t(p[,i]-step13/size); step14<-step14+dif }
  cov.beta<-solve(step11)%*%(step12)%*%solve(step11)
  return(list(cov.beta=cov.beta))
}
