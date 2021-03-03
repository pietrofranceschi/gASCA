#' @importFrom Matrix rankMatrix
#' @importFrom parallel mclapply
#' @import stats 
#' @import graphics

# ------------------------------------------------------------------------
## Auxiliary Functions for singular value decomposition
# -----------------------------------------------------------------------

ASCA_svd_engine <- function(a) {
  rk <- Matrix::rankMatrix(a)
  if (rk[1] == 1) {
    stop("One of the components has rank 1")
  }
  sv <- svd(scale(a,scale = FALSE), nv = rk[1], nu = rk[1])
  scores <- sv$u %*% diag(sv$d[1:rk[1]])
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
                variance = sum(variance))
    )
}


combine_terms <- function(terms,comb){
  df_terms <- as.data.frame(terms)
  new_term_name <- paste(comb,collapse = "+")
  df_terms[[new_term_name]] <- Reduce('+',df_terms[,comb])
  df_terms <- df_terms[,!(colnames(df_terms) %in% comb)]
  terms <- as.matrix(df_terms)
  return(terms)
}


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
                variances = t(sapply(svs, function(x) x$variance))
    )
    )
  }
}


# ------------------------------------------------------------------------
## Auxiliary Functions for glm decomposition
# -----------------------------------------------------------------------

## This function uses glm to decompose a data matrix 

ASCA_decompose <- function(d,x,f, comb = NULL, glm_par = vector(mode = "list", length = 0)){
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
    err <- stats::residuals(m) 
    
    ## here we combine (sum) the terms which are present in the "combine" argument
    if (!is.null(comb)) {
      terms <- combine_terms(terms,comb)
    }
    
    return(list(dec = cbind(terms = terms, error = err),
                deviance = stats::deviance(m)))
  })
  
  ## organize the terms in a 3d array
  termarray <- sapply(models, function(m) m$dec, simplify = "array")
  dimnames(termarray) <- list(rownames(x),
                              colnames(models[[1]]$dec),
                              colnames(x))
  
  ## reset the contrasts!
  options("contrasts" = oldcontrasts$contrasts)
  return(list(decomposition = termarray,
              exp_variation = apply(termarray,2,function(t) norm(t, type = "F")),
              deviance = sapply(models,function(m) m$deviance),
              d = d,
              x = x,
              f = f,
              glm_par = glm_par,
              comb = comb,
              validation = NULL
              ))
}


# ------------------------------------------------------------------------
## ASCA S3 class, methods
# -----------------------------------------------------------------------

## This wrapper performs the complete ASCA decomposition and create the ASCA S3 object

#' ASCA
#' @md 
#' @description 
#' This function performs the full ASCA decomposition and acts as a constructor of the ASCA S3 object
#'
#' @details 
#' The ASCA decomposition of a data matrix is performed by using Generalized Linear Models to perform 
#' the estimation of the expected values. the use of GLM's allows the extension of the method to non normal data 
#' and unbalanced designs 
#'
#' @param d a data frame with the experimental design, which associate each sample to the design factors 
#' @param x a matrix with the data measured on the different samples
#' @param f a string with the formula which will be used for the `glm` decomposition
#' @param comb a character vector specifying which terms of the decomposition should be combined
#' @param glm_par a list with the parameters to be passed to the `glm` call
#'
#' @return
#' An object of class `asca`. Actually a list holding the results:
#' 
#'* of the `glm` decomposition
#'* of the `svd` of the factor matrices
#'* of the validation of the object (if present)
#'  
#' @note 
#' There are specific methods to handle and extract information from the `asca` object
#'* `print`
#'* `summary`
#'* [scores] 
#'* [loadings]
#'* [varimp]
#'* [plot.asca]
#'* [validate]
#'
#' @examples
#' ## load the data
#' data("synth_count_data")
#' 
#' ## perform the ASCA decomposition
#' ASCA_test <- ASCA(
#' d = synth_count_data$design,
#' x = synth_count_data$counts, 
#' f = "time + treatment + time:treatment",
#' glm_par = list(family = poisson())
#' )
#' 
#' ## decomposition with a combination of treatment and interaction
#' ASCA_test_comb <- ASCA(
#' d = synth_count_data$design,
#' x = synth_count_data$counts, 
#' f = "time + treatment + time:treatment",
#' comb = c("treatment","time:treatment"),
#' glm_par = list(family = poisson())
#' )
#' 
#' ## Print a summary of the decomposition
#' print(ASCA_test)
#' summary(ASCA_test)
#' 
#' @export
#' 
ASCA <- function(d,x,f, comb = NULL, glm_par = vector(mode = "list", length = 0)) {
  asca <- ASCA_decompose(d = d, 
                        x = x,
                        f = f,
                        comb = comb,
                        glm_par = glm_par)
  
  nterms <- dim(asca$decomposition)[2]
  asca$svd <- ASCA_svd(asca$decomposition[,2:nterms,])
  attr(asca, "class") <- "asca"
  return(asca)
}

## ------------------------------------------------------------------------------------------------

## #' @rdname ASCA
## #' @export
## print <- function(x) {
##   UseMethod("print")
## }


#' @export
#' @rdname ASCA
print.asca <- function(x) {
  cat("ASCA decomposition","\n")
  cat("\n")
  cat("Number of samples ", nrow(x$x), "\n")
  cat("Number of variables", ncol(x$x), "\n")
  cat("\n")
  cat("Formula: ",x$f,"\n")
  cat("Terms: ")
  nterms <- dim(x$decomposition)[2]
  cat(paste(dimnames(x$decomposition)[[2]][2:nterms],collapse = ",")); cat("\n")
}

## #' @rdname ASCA
## #' @export
## summary <- function(x) {
##   UseMethod("summary")
## }


#' @export
#' @rdname ASCA
summary.asca <- function(x) {
  cat("ASCA decomposition","\n")
  cat("\n")
  cat("Number of samples ", nrow(x$x), "\n")
  cat("Number of variables", ncol(x$x), "\n")
  cat("\n")
  cat("Formula: ",x$f,"\n")
  cat("Terms: ")
  nterms <- dim(x$decomposition)[2]
  cat(paste(dimnames(x$decomposition)[[2]][2:nterms],collapse = ",")); cat("\n")
  is_validated <- ifelse(length(x$validation) == 0, FALSE, TRUE)
  cat("Validation: ",is_validated,"\n")
  cat("\n")
  cat("glm parameters: \n")
  print(x$glm_par)
}


#' @title scores
#' @description Extract the scores from an `asca` object
#' @param x an object of class asca
#' @return 
#' A list of matrices with the scores of the decomposition. The number of eigenvectors 
#' is determined by the number of levels for each factor
#' @export
scores <- function(x) {
  UseMethod("scores")
}

#' @export
scores.asca <- function(x){
  return(x$svd$scores)
}

#' @title loadings
#' @description Extracts the loadings of the different terms from the `asca` object 
#' @param x an object of class asca
#' @return A list of matrices with the loadings of the decomposition. The number of eigenvectors 
#' is determined by the number of levels for each factor
#' @export
loadings <- function(x) {
  UseMethod("loadings")
}

#' @export
loadings.asca <- function(x){
  return(x$svd$loadings)
}

#' @title varimp
#' @description  Extracts the loadings of the different terms from the `asca` object 
#' @details 
#' The variable importance is calculated by weighting the absolute value of the loadings of the 
#' singular values by the relative importance of the eigenvalues
#' @param x the asca object
#' @return A matrix with the variable importance for each factor of the decomposition
#' @export
varimp <- function(x) {
  UseMethod("varimp")
}

#' @export
varimp.asca <- function(x){
  return(x$svd$varimp)
}
  


#' @title getValidation.asca
#' @description  Extracts the results of the validation from the `asca` object 
#' @param x the asca object
#' @return A list with the resampling quantiles for the loadings and the variable importance
#' @export
getValidation <- function(x) {
  UseMethod("getValidation")
}

#' @export
getValidation <- function(x){
  if (is.null(x$validation)) {stop("The asca object does not contain validation results")}
  return(x$validation)
}





