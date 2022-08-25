
#' ASCA_svd
#' 
#' @description 
#' This function is designed to perform the term-wise decomposition by R base `prcomp`
#' 
#' @param decomposition the decomposition list produced by  `ASCA_decompose`.
#' @param ... additional parameters to be passed to `prcomp`
#'
#' @return a list with the results of the SVD decomposition for each term. Each element is the output of a  
#' `prcomp` call. 
#' 
#'  
#' @examples
#' ## load the data
#' data("synth_count_data")
#' 
#' ## perform the ASCA decomposition
#' dec_test <- ASCA_decompose(
#' d = synth_count_data$design,
#' x = synth_count_data$counts, 
#' f = "time + treatment + time:treatment",
#' glm_par = list(family = poisson())
#' )
#' 
#' ## Run the SVD of the terms
#' svd_test <- ASCA_svd(dec_test$decomposition)
#'  
#' @export
ASCA_svd <- function(decomposition, ...){
  pcas <- lapply(decomposition,stats::prcomp,...)

  return(pcas)
  }


# ------------------------------------------------------------------------
## functions for glm decomposition
# -----------------------------------------------------------------------

#' ASCA_decompose
#' 
#' @description 
#' This function performs the full ASCA decomposition
#' @md
#' @details 
#' The ASCA decomposition of a data matrix is performed by using Generalized Linear Models 
#' to estimate univariate expected values. The use of GLM's allows the extension of the method to non normal data 
#' and unbalanced designs. 
#' This function performs only the decomposition without the SVD, which have to be performed by `ASCA_svd`.
#' The level of fit for each variable is assessed calculating the `pseudoR2`, which is defined as:
#' 
#' \deqn{1-residual_deviance/null_deviance}
#' 
#' The variable importance element stores a measure of the importance of each variable \eqn{c} for each term calculated 
#' as the norm of each column. L2 norm is also calculated for each decomposition term.
#' 
#' It is important to highlight that in the case of count data with large fraction of zeroes the variable importance and the 
#' term L2 norms cannot be considered as reliable measure of importance because in presence of log links very low expected values are associated to 
#' very large negative values in the linear predictor space.
#' 
#' @param d a data.frame/matrix with the design
#' @param x a data.frame/matrix of numeric values to be decomposed 
#' @param f a string holding the formula of the decomposition
#' @param glm_par a list with the parameters to be passed to the `glm` call
#' @param res_type the types of GLM residuals
#'
#' @md
#' @return a list with the full outcomes of the decomposition with the following elements
#' 
#' * decomposition: a list holding the results of the decomposition
#' * mu: a vector with the constant terms of the univariate models
#' * residuals: a matrix holding the model residuals. Their type is stored in the `res_type` element.
#' * prediction:  the matrix with the predicted values in the linear predictor space
#' * pseudoR2: a parameter to assess the goodness of fit for the model on each variable. 
#' * glm_par: a list with the parameters used for modelling
#' * res_type: the type of residuals
#' * varimp: the importance of the individual variables in the decomposition terms
#' * terms_L2: the L2 norm of the individual terms back transformed in the response space
#' * d: a data.frame with the design
#' * x: a data.frame with the initial data
#' * f: the string defining the decomposition
#' * combined: a vector holding the combined terms ()
#' * invlink: the inverse of the link function used in the glm fitting 
#' 
#' 
#' 
#' @examples
#' ## load the data
#' data("synth_count_data")
#' 
#' ## perform the ASCA decomposition
#' dec_test <- ASCA_decompose(
#' d = synth_count_data$design,
#' x = synth_count_data$counts, 
#' f = "time + treatment + time:treatment",
#' glm_par = list(family = poisson())
#' )
#' 
#' 
#' @export
#' @import parallel
#' 


ASCA_decompose <- function(d,x,f, 
                           glm_par = vector(mode = "list", length = 0),
                           res_type = "response"){
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
    
    ##  here the type of residuals should be set, now is response
    residuals <- stats::residuals(m, type = res_type) 
    
    
    return(list(dec = terms,
                mu = mu,
                pseudoR2 = 1-(m$deviance/m$null.deviance),
                residuals = residuals,
                prediction = stats::predict(m, type = "link"),
                linkinv = m$family$linkinv))
  })
  
  ## organize the terms in a 3d array
  termarray <- vapply(models, function(m) m$dec, FUN.VALUE = models[[1]]$dec)
  dimnames(termarray) <- list(rownames(x),
                              colnames(models[[1]]$dec),
                              colnames(x))
  
  termlist <- lapply(colnames(models[[1]]$dec), function(t) termarray[ ,t,])
  names(termlist) <- colnames(models[[1]]$dec)
  
  
  ## residuals and predictions
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
  
  
  output <- list(decomposition = termlist,
                 mu = mu_vec,
                 residuals = resid_mat,
                 prediction = predict_mat,
                 pseudoR2 = pseudoR2_vec,
                 glm_par = glm_par,
                 res_type = res_type,
                 varimp = ASCA_get_varimp(termlist),                                     
                 terms_L2 = vapply(termlist,function(t) ASCA_norm(t),
                                 FUN.VALUE = double(1)),                      
                 d = d,
                 x = x,
                 f = f,
                 combined = NULL,
                 linkinv = models[[1]]$linkinv
  )
  
  
  return(output)
}


#' ASCA_combine_terms
#'
#' @description The function combines two or more terms of the ASCA decomposition 
#' 
#' @param asca the results of `ASCA_decompose`
#' @param comb a vector of character stings with the names of the terms to combine
#'
#' @return a new ASCA decomposition with updated terms and terms derived quantities. 
#' The `combined` element keep tracks of the combined factors
#' 
#' @export
#'
#' @examples
#' 
#' ## load the data
#' data("synth_count_data")
#' 
#' ## perform the ASCA decomposition
#' dec_test <- ASCA_decompose(
#' d = synth_count_data$design,
#' x = synth_count_data$counts, 
#' f = "time + treatment + time:treatment",
#' glm_par = list(family = poisson())
#' )
#'
#' ## combine two terms
#' comb_terms <- ASCA_combine_terms(dec_test, c("treatment","time:treatment"))
#' 
#' ## validate the combined decomposition
#' comb_test_val <- ASCA_permutation(comb_terms, 10)
#' 
#' 

ASCA_combine_terms <- function(asca,comb){
  
  new_term_name <- paste(comb,collapse = "+")
  combterm <- do.call('+',asca$decomposition[comb])
  
  new_decomposition_list <- asca$decomposition
  new_decomposition_list[[new_term_name]] <- combterm
  
  new_asca <- asca
  
  new_asca$decomposition <- new_decomposition_list
  new_asca$varimp <-  ASCA_get_varimp(new_decomposition_list)
  new_asca$terms_L2 <-  vapply(new_decomposition_list,function(t) ASCA_norm(t),
                               FUN.VALUE = double(1))    
  new_asca$combined <- comb
  
  return(new_asca)
}




# ------------------------------------------------------------------------
## functions to calculate the variable importance
# -----------------------------------------------------------------------


#' ASCA_get_varimp
#'
#' @description The function calculates the variable importance from the decomposition array. 
#' See \link[gASCA]{ASCA_decompose} for additional details on how the variable importance is calculated
#'
#'
#' @param termlist The list of matrices holding the results of an ASCA decomposition
#' @param ftrans The function which will be applied to the input before calculating the norm
#'
#' @return a matrix with the variable importance for each term (see \link[gASCA]{ASCA_decompose})
#' 


ASCA_get_varimp <- function(termlist,ftrans = identity) {
  out <- vapply(termlist, 
                function(t) apply(t,2, function(c) ASCA_norm(c,ftrans = ftrans)), 
                FUN.VALUE = double(ncol(termlist[[1]])))
  return(out)
}



# ------------------------------------------------------------------------
## functions to calculate the matrix norm
# -----------------------------------------------------------------------


#' ASCA_norm
#' 
#' @description The function calculates the norm of a vector (standatd norm) or of a matrix (L2 norm) 
#' after the transformation by f_trans
#' 
#' 
#' @param m the matrix/vector
#' @param ftrans the function which will be applied to the input before calculating the norm
#'
#' @return the norm of the matrix/vector
#' 
#' 


ASCA_norm <- function(m, ftrans = identity) {
  m1 <- ftrans(m)
  if (is.matrix(m1)){
    return(norm(m1, type = "F"))
  } else {
    return(norm(m1, type = "2"))
  }
}





# ------------------------------------------------------------------------
## functions for validation
# -----------------------------------------------------------------------

#' ASCA_permute_design
#' 
#' @description The function performs a random permutation of a design data frame.
#' All columns are permuted independently
#'
#' 
#' @param d a design data.frame 
#'
#' @return a data.frame
#' 


ASCA_permute_design <- function (d){
  perm_d <- apply(d, 2L, sample)
  perm_d
}


#' ASCA_permutation
#'
#' @description the function implements a design permutation strategy to validate the results of an ASCA decomposition
#'
#' @param asca an object of class `ASCA_decomposition` 
#' @param turns the number of permutations
#' @param qt the empirical quantile
#' @md
#' @return A list with the following elements
#' 
#' * `R2_qt`: a matrix with the quantiles for the pseusoR2 of the individual variables
#' * `L2_qt`: a matrix with the quantiles for the L2 norm of the individual terms
#' * `varimp_qt`: a list holding the matrices of quantiles for each variable in each term
#'  
#' @export
#'
#' @examples
#' 
#' ## load the data
#' data("synth_count_data")
#' 
#' ## perform the ASCA decomposition
#' dec_test <- ASCA_decompose( d = synth_count_data$design,
#' x = synth_count_data$counts, 
#' f = "time + treatment + time:treatment",
#' glm_par = list(family = poisson())
#' )
#' 
#' ## validate the outcomes
#' 
#' val_out <- ASCA_permutation(dec_test, 10)
#' 
#' 


ASCA_permutation <- function(asca,
                             turns = 100,
                             qt = 0.95){
  
  ## objects to collect the outcomes
  random_ascas_R2 <- matrix(0,ncol = length(asca$pseudoR2),nrow = turns)
  colnames(random_ascas_R2) <- names(asca$pseudoR2)
  
  random_ascas_L2 <- matrix(0, ncol = length(asca$terms_L2), nrow = turns)
  colnames(random_ascas_L2) <- names(asca$terms_L2)
  
  random_ascas_varimp <- array(0,c(dim(asca$varimp),turns))
  random_ascas_varimp_names <- dimnames(asca$varimp)
  random_ascas_varimp_names[[3]] <- NULL
  dimnames(random_ascas_varimp) <- random_ascas_varimp_names
  
  ## run the battery ascas
  for (i in 1:turns){
    cat("Permutation #",i,"\n")
    d_perm <- ASCA_permute_design(asca$d)
    running_asca <- ASCA_decompose(d = d_perm,
                                   x = asca$x,
                                   f = asca$f,
                                   glm_par = asca$glm_par)
    if (!is.null(asca$combined)) {
      running_asca <- ASCA_combine_terms(running_asca,asca$combined)
    }
    random_ascas_R2[i,] <- running_asca$pseudoR2
    random_ascas_L2[i,] <- running_asca$terms_L2
    random_ascas_varimp[,,i] <- running_asca$varimp
  } 
  
  ## organize the outputs
  
  R2_qt <- apply(random_ascas_R2,2,stats::quantile,probs = qt)
  L2_qt <- apply(random_ascas_L2,2,stats::quantile,probs = qt)
  varimp_qt <- apply(random_ascas_varimp,c(1,2),stats::quantile,probs = qt)
  
  ## variable importance here is an array, let's make it a list ...
  
  varimp_qt_list_names <- as.list(dimnames(varimp_qt)[[1]])
  names(varimp_qt_list_names) <- dimnames(varimp_qt)[[1]]
  
  # varimp_qt_list <- lapply(varimp_qt_list_names, function(n) varimp_qt[n,,])
  
  
  out <- list("R2_qt" = R2_qt,
              "L2_qt" = L2_qt,
              "varimp_qt" = varimp_qt)
  
  return(out)
}


# ------------------------------------------------------------------------
## Subsetting variables
# -----------------------------------------------------------------------


#' ASCA_trim_vars
#'
#' @param decomposition a list of decomposition matrices produced by `ASCA_decompose` or `ASCA_combine_terms`
#' @param keep a vector of numbers or names indicating the variables which should be kept
#'
#' @return a list of reduced matrices 
#' 
#' @export
#' 
#' @examples
#' 
#' ## load the data
#' data("synth_count_data")
#' 
#' ## perform the ASCA decomposition
#' dec_test <- ASCA_decompose( d = synth_count_data$design,
#' x = synth_count_data$counts, 
#' f = "time + treatment + time:treatment",
#' glm_par = list(family = poisson())
#' )
#' 
#' ## keep only the forst six species 
#' 
#' trimmed_dec <- ASCA_trim_vars(dec_test$decomposition, 1:6)
#' 
#' 


ASCA_trim_vars <- function(decomposition, keep = NULL) {
  if (is.null(keep)) {
    return(decomposition)
  } else {
    out <- lapply(decomposition, function(t) t[,keep])
    return(out)
  }
}





