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
  ## calculate the variable importance combining loadings and explained variances
  ## variable importance is calculated combining the information contained in the 
  ## relevant components
  varimp <- rowSums(sweep(abs(loadings),2,variance,"*")/sum(variance))
  
  ## assign names to rows and columns
  names(varimp) <- colnames(a)
  dimnames(loadings) <- list(colnames(a),paste0("eig",1:rk[1]))
  dimnames(scores) <- list(rownames(a),paste0("eig",1:rk[1]))
  
  ## pack everything into a list
  return(list(scores = scores,
                loadings = loadings,
                varimp = varimp,
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
#'  * scores: a list with the scores for each term, in the linear predictor space
#'  * loadings: a list with the scores for each term
#'  * varimp: a list with the importance of each variable for each term
#'  * variances: a list with the variance of each eigenvector for each term 
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
                varimp = sapply(svs, function(x) x$varimp),
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
#' * deviance: the array containing the deviance of the models
#' * glm_par: a list with the parameters used for modeling
#' * res_type: the type of residuals
#' 
#' 
#' 
#' @examples
#' ## load the data
#' data("synth_count_data")
#' 
#' ## perform the ASCA decomposition
#' ASCA_test <- ASCA_decompose(
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
    terms <- cbind(mu = stats::coefficients(m)["(Intercept)"], stats::predict(m, type = "terms"))
    ## the matrix of terms is centered, so we need add the constant
    ## the mu term represents the "grand mean"
    
    ##  here the type of residuals should be set, now is the default (the deviance residuals)
    err <- stats::residuals(m, type = res_type) 
    
    
    return(list(dec = cbind(terms = terms, error = err),
                deviance = stats::deviance(m)))
  })
  
  ## organize the terms in a 3d array
  termarray <- vapply(models, function(m) m$dec, FUN.VALUE = models[[1]]$dec)
  dimnames(termarray) <- list(rownames(x),
                              colnames(models[[1]]$dec),
                              colnames(x))
  
  ##here we scale the error parameter
  termarray[,"error",] <- scale(termarray[,"error",])
  
  
  ## reset the contrasts!
  options("contrasts" = oldcontrasts$contrasts)
  return(list(decomposition = termarray,
              deviance = sapply(models,function(m) m$deviance),
              glm_par = glm_par,
              res_type = res_type
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



