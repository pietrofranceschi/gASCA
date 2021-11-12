#' @importFrom Matrix rankMatrix
#' @importFrom parallel mclapply
#' @import stats 
#' @import graphics

# ------------------------------------------------------------------------
## Auxiliary Functions for singular value decomposition
# -----------------------------------------------------------------------

ASCA_svd_engine <- function(a) {
  rk <- Matrix::rankMatrix(a)

  ## Ok now we have to handle the pathological case of a matrix of rank one
  ## with equal values on the columns (that ASCA invariably shows up as the mu term)
  
  if ((rk[1] == 1) & (diff(range(a[,1])) < .Machine$double.eps ^ 0.5)) {
    return(NULL)
  }
  sv <- svd(scale(a,scale = FALSE), nv = rk[1], nu = rk[1])
  ## here we have to handle the construction of the diagonal matrix 
  
  scores <- sv$u %*% diag(sv$d[1:rk[1]], nrow = rk[1])
  loadings <-  sv$v
  variance <-  sv$d[1:rk[1]]^2
 
  ## assign names to rows and columns
  dimnames(loadings) <- list(colnames(a),paste0("eig",1:rk[1]))
  dimnames(scores) <- list(rownames(a),paste0("eig",1:rk[1]))
  
  ## pack everything into a list
  return(list(scores = scores,
                loadings = loadings,
                variance = variance)
    )
}



#' ASCA_svd
#'
#' @param arr the multidimensional (2D or 3D) array. In 2D case rows are samples and columns variables. 
#' In the 3D case (the output of `ASCA_decompose`) organized as samples x terms x variables 
#' the SVD is performed along the 2nd dimension, aka one SVD per term.
#' If the matrix has rank 1 the decompoition is set to NULL
#'
#' @return a list with the results of the SVD decomposition. 
#' 
#'  * scores: a list with the scores for each term, in the linear predictor space
#'  * loadings: a list with the scores for each term
#'  * variances: a list with the variance of each eigenvector for each term 
#'  
#'  
#'  
#'@examples
#' ## load the data
#' data("synth_count_data")
#' 
#' ## perform the ASCA decomposition
#' decomposition_test <- ASCA_decompose(
#' d = synth_count_data$design,
#' x = synth_count_data$counts, 
#' f = "time + treatment + time:treatment",
#' glm_par = list(family = poisson())
#' )
#' 
#' ## Run the SVD of the terms
#' svd_test <- ASCA_svd(decomposition_test$decomposition)
#'  
#'  
#' @export
#'
#' 

ASCA_svd <- function(arr){
  ## check the size of the arr
  if (length(dim(arr)) == 2) {
    return(ASCA_svd_engine(arr))
  } else {
    svs <- apply(arr,2, function(mat) ASCA_svd_engine(mat))
    ## add the term name to the output list
    names(svs) <- dimnames(arr)[[2]]
    ## format the output
    return(list(scores = lapply(svs, function(x) x$scores),
                loadings = lapply(svs, function(x) x$loadings),
                variances = lapply(svs, function(x) x$variance)
    )
    )
  }
}


#' ASCA_svd_response
#'
#' @param arr the multidimensional (2D or 3D) array. In 2D case rows are samples and columns variables. 
#' In the 3D case (the output of `ASCA_decompose`) organized as samples x terms x variables 
#' the SVD is performed along the 2nd dimension, aka one SVD per term.
#' If the matrix has rank 1 the decompoition is set to NULL
#' @param invlink the function to be applied to the scores of the terms that are not NULL or error
#'
#' @return a list with the results of the SVD decomposition. 
#'  * scores: a list with the scores for each term, in the response
#'  * loadings: a list with the scores for each term
#'  * varimp: a list with the importance of each variable for each term
#'  * variances: a list with the variance of each eigenvector for each term 
#' @export
#'
#' 



ASCA_svd_response <- function(arr, invlink = identity) {
  svdout <- ASCA_svd(arr)
  unchangedid <- c(names(which(sapply(svdout$scores, is.null))), "error")
  term_s <- names(svdout$scores)
  changeids <- term_s[!(term_s %in% unchangedid)]
  scores_response <- svdout$scores
  scores_response[changeids] <- lapply(scores_response[changeids], invlink)
  svdout$scores <- scores_response
  return(svdout)
}



# ------------------------------------------------------------------------
## functions for glm decomposition
# -----------------------------------------------------------------------



#' ASCA_decompose
#' 
#' @description 
#' This function performs the full ASCA decomposition
#'
#' @details 
#' The ASCA decomposition of a data matrix is performed by using Generalized Linear Models to perform 
#' the estimation of the expected values. The use of GLM's allows the extension of the method to non normal data 
#' and unbalanced designs. 
#' This function performs only the decomposition without the SVD, which have to be performed by `ASCA_svd`
#' The `pseudoR2` is calculated as 1-(residual deviance/null deviance). The variable importance list contains 
#' two data.frames. The `wterm_varimp` stores measure of the importance of each variable \eqn{c} for each term calculated as 
#' the squared norm of each column vector divided by the squared Frobenius norm of the term 
#'  
#' \deqn{Varimp(A,c) = |c|^{2}/\|A|_{F}^{2}}
#' 
#' Since the terms are mean centered, the previous quantity is equivalent to the variance of variable \eqn{c} (for that term) 
#' normalized by the sum of the variances of all the variables.
#' The `bterm_varimp` is instead a measure the weight of each variable across the different decomposition terms. 
#' Also in this case, the squared norm of the columns vectors relative to \eqn{c} is considered, but now it is normalized 
#' to the sum of the squared norms of the column vectors associated to \eqn{c} across the different terms. 
#' It is important to remember that the second quantity should be combined with a measure of model fit to obtain a meaningful measure of
#' the variable importance across the different terms.   
#' 
#' 
#' @md
#' @param d a data.frame/matrix with the design
#' @param x a data.frame/matrix of numeric values to be decomposed 
#' @param f a string holding the formula of the decomposition
#' @param glm_par a list with the parameters to be passed to the `glm` call
#' @param res_type the types of GLM residuals
#'
#' @return a list with the full outcomes of the decomposition with the following elements
#' 
#' * decomposition: the 3D array holding the results of the decomposition
#' * mu: a vector with the constant terms of the univariate models
#' * residuals: a matrix holding the model residuals
#' * prediction:  the matrix with the predicted values in the linear predictor space
#' * pseudoR2: a parameter to assess the goodness of fit for the model on each variable. 
#' * glm_par: a list with the parameters used for modeling
#' * res_type: the type of residuals
#' * varimp: a list composed by two data.frames `wterm_varimp` and  `bterm_varimp`. see details
#' 
#' 
#' 
#' @examples
#' ## load the data
#' data("synth_count_data")
#' 
#' ## perform the ASCA decomposition
#' decomposition_test <- ASCA_decompose(
#' d = synth_count_data$design,
#' x = synth_count_data$counts, 
#' f = "time + treatment + time:treatment",
#' glm_par = list(family = poisson())
#' )
#' 
#' 
#' @export
#'
#' 


ASCA_decompose <- function(d,x,f, 
                           glm_par = vector(mode = "list", length = 0),
                           res_type = "pearson"){
  ## contrasts should be changed so first I save the old ones
  oldcontrasts <- options("contrasts")
  ## then I set the new ones
  options(contrasts = rep("contr.sum", 2))
  ## calculate the models for all variables
  ## return the terms calculated
  myform = paste0("y ","~",f)
  
  ## prepare the list of model parameters
  glm_par_run <- glm_par
  glm_par_run$formula <-  myform
  
  
  ## transform x and d in data frames ...
  d <-as.data.frame(d)
  x <-as.data.frame(x)
  
  ## check if d has numeric columns
  if (sum(vapply(d,is.numeric, FUN.VALUE = vector("logical",1)))) {
    stop("One of the columns of d is numeric, transform it to a factor")
  }
  
  

  ## run the models and format the decomposition
  models <- parallel::mclapply(1:ncol(x), function(v){
    df <- cbind(y = x[,v],d)
    glm_par_run$data <- df
    ## run the modeling
    m <- do.call(stats::glm,glm_par_run)
    ## the terms in form of matrix
    terms <- stats::predict(m, type = "terms")
    ## organize the constant term
    mu <-  attr(terms,"constant")
  
    ## the matrix of terms is centered, so we need add the constant
    ## the mu term represents the "grand mean"
    
    ##  here the type of residuals should be set, now is the default (the deviance residuals)
    residuals <- stats::residuals(m, type = res_type) 
    
    
    return(list(dec = terms,
                mu = mu,
                pseudoR2 = 1-(m$deviance/m$null.deviance),
                residuals = residuals,
                prediction = stats::predict(m, type = "link")))
  })
  
  ## organize the terms in a 3d array
  termarray <- vapply(models, function(m) m$dec, FUN.VALUE = models[[1]]$dec)
  dimnames(termarray) <- list(rownames(x),
                              colnames(models[[1]]$dec),
                              colnames(x))
  resid_mat <- vapply(models, function(m) m$residuals, FUN.VALUE = models[[1]]$residuals)
  predict_mat <- vapply(models, function(m) m$prediction, FUN.VALUE = models[[1]]$prediction)
  
  
  ## pseudoR2
  pseudoR2_vec <- vapply(models,function(m) m$pseudoR2, FUN.VALUE = double(1))
  names(pseudoR2_vec) <- colnames(x)
  
  ## mu
  mu_vec <- vapply(models,function(m) m$mu, FUN.VALUE = double(1))
  names(mu_vec) <- colnames(x)
  
  ## add the variable names to the residuals
  colnames(resid_mat) <- colnames(x)
  
  ## add the variable names to the predictions
  colnames(predict_mat) <- colnames(x)
  
  ## calculate the variable importance from the decomposition
  
  
  ## reset the contrasts!
  options("contrasts" = oldcontrasts$contrasts)
  return(list(decomposition = termarray,
              mu = mu_vec,
              residuals = resid_mat,
              prediction = predict_mat,
              pseudoR2 = pseudoR2_vec,
              glm_par = glm_par,
              res_type = res_type,
              varimp = get_varimp(termarray)
              ))
}


# ------------------------------------------------------------------------
## functions for bootstrap and validations
# -----------------------------------------------------------------------



#' ASCA_do_bootstrap
#'
#' @description The function performs a "design" aware bootstrap of a data matrix.
#'
#' @param d a dataframe with the design
#' @param x the data matrix
#'
#' @return a list containing the initial design and the bootstrap data matrix
#' @export
#'
#' @examples
#' 
ASCA_do_bootstrap <- function(d,x){
  ## here we construct a factor to associate each sample to a level of
  ## the design
  
  ## check if d has numeric columns
  if (sum(vapply(d,is.numeric, FUN.VALUE = vector("logical",1)))) {
    stop("One of the columns of d is numeric, transform it to a factor")
  }
  
  if (ncol(d) == 1) {
    factid <- d[,1]
  } else {
    factid <- apply(d,1,paste,collapse = "_")
  }
  ## now we split the data by the factor
  splitx <- split(as.data.frame(x),factid)
  ## create the bootstrap sample
  splitx_bootstrap <- lapply(splitx, function(i) i[sample(1:nrow(i), nrow(i), replace = TRUE),])
  ## rejoin the x
  x_boot <- unsplit(splitx_bootstrap,factid)
  return(list(x_val = x_boot,
              d_val = d))
}


#' ASCA_do_jacknife
#'
#' @description The function performs a "design" aware jacknifing of a data matrix.
#'
#' @param d a dataframe with the design
#' @param x the data matrix
#' @param nout the number of samples for left out for each design level   
#'
#' @return a list containing the initial design and the reduced data matrix
#' @export
#'
#' @examples
#' 

ASCA_do_jacknife <- function(d,x, nout = 1){
  
  ## here we construct a factor to associate each sample to a level of
  ## the design
  ## we need to add an error if the nout is bigger than the minimum size of the classes
  
  
  ## check if d has numeric columns
  if (sum(vapply(d,is.numeric, FUN.VALUE = vector("logical",1)))) {
    stop("One of the columns of d is numeric, transform it to a factor")
  }
  
  
  if (ncol(d) == 1) {
    factid <- d[,1]
  } else {
    factid <- apply(d,1,paste,collapse = "_")
  }
  ## now we split the data by the factor
  splitx <- split(as.data.frame(x),factid)
  splitd <- split(d,factid)
  ## create the jacknife
  splitx_jk <- lapply(splitx, function(i) i[sample(1:nrow(i), nrow(i)-nout),])
  splitd_jk <- lapply(splitd, function(i) i[sample(1:nrow(i), nrow(i)-nout),])
  ## rejoin the x
  ## here is differnt from bootstrap because also d should be smaller
  x_jk <- do.call(rbind,splitx_jk)
  rownames(x_jk) <- NULL
  d_jk <- do.call(rbind,splitd_jk)
  rownames(d_jk) <- NULL
  return(list(x_val = x_jk,
              d_val = d_jk))
}



#' ASCA_combine_terms
#'
#' @param arr3d 
#' @param comb 
#'
#' @return
#' @export
#'
#' @examples
#' 

ASCA_combine_terms <- function(arr3d,comb){
  
  combine_engine <- function(arr2d,comb){
    df_terms <- as.data.frame(arr2d)
    new_term_name <- paste(comb,collapse = "+")
    df_terms[[new_term_name]] <- Reduce('+',df_terms[,comb])
    df_terms <- df_terms[,!(colnames(df_terms) %in% comb)]
    out <- as.matrix(df_terms)
    return(out)
  }
  
  out <- apply(arr3d,3,function(x) combine_engine(x,comb), simplify = FALSE)
  return(simplify2array(out))
}


#' ASCA_fix_sign
#'
#' @description The function flips the sign of B if the scalar product with A is bigger than 90 degrees
#' 
#' @param A the reference vector
#' @param B the second vector
#'
#' @return
#' @export
#'
#' @examples
#' 
ASCA_fix_sign <- function(A,B){
  mult <- sign(diag(t(A) %*% B))
  flipB <- t(apply(B,1, function(r) r*mult))
  return(flipB)
}


# ------------------------------------------------------------------------
## functions to calculate the variable importance
# -----------------------------------------------------------------------


#' get_varimp
#'
#' @description The function calculates the variable importance from the decomposition array. 
#'
#'
#' @param termarray 
#'
#' @return
#' @export
#'
#' @examples
#' 

get_varimp <- function(termarray) {
  termlist <- as.list(dimnames(termarray)[[2]])
  termlist <- lapply(termlist, function(t) scale(termarray[,t,], scale = FALSE))
  names(termlist) <- dimnames(termarray)[[2]]
  varimpmat <- vapply(termlist, function(t) colSums(abs(t)^2), FUN.VALUE = double(10))
  pippo <- colSums(varimpmat)
  pluto <- colSums(t(varimpmat))
  return(list(
    "wterm_varimp" = as.data.frame(sweep(varimpmat,2,pippo,"/")),
    "bterm_varimp" = as.data.frame(sweep(t(varimpmat),2,pluto,"/"))))
}

