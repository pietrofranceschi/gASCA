
# ------------------------------------------------------------------------
## Auxiliary Functions for rotations
# -----------------------------------------------------------------------



procrustes <- function(A, B){
  # center and normalize A 
  A.centered <- t(scale(t(A), center = TRUE, scale = FALSE))
  A.size <- norm(A.centered, type = "F") / (ncol(A) * nrow(A))
  #A.normalized <- A.centered / A.size
  A.normalized <- A.centered 
  
  # center and normalize B
  B.centered <- t(scale(t(B), center = TRUE, scale = FALSE))
  B.size <- norm(B.centered, type = "F") / (ncol(B) * nrow(B))
  #B.normalized <- B.centered / B.size
  B.normalized <- B.centered 
  
  # Rotation matrix T 
  svd.results <- svd(B.normalized %*% t(A.normalized))
  U <- svd.results$u
  V <- svd.results$v
  T <- V %*% t(U)
  
  # B transformed
  B.transformed <- T %*% B.normalized
  
  # Error after superimposition
  RSS <- norm(A.normalized - B.transformed,  type = "F")
  
  # Return
  return(list(A.normalized = A.normalized, 
              B.normalized = B.normalized, 
              rotation.mtx = T, 
              B.transformed = B.transformed, 
              RSS = RSS))
}



fix_sign <- function(A,B){
  mult <- sign(diag(t(A) %*% B))
  flipB <- t(apply(B,1, function(r) r*mult))
  return(flipB)
}



# ------------------------------------------------------------------------
## Auxiliary Functions for validation
# -----------------------------------------------------------------------
## both methods are returning a list with x_val and d_val elements



do_bootstrap <- function(d,x){
  ## here we construct a factor to associate each sample to a level of
  ## the design
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


do_jacknife <- function(d,x, nout = 1){
  
  ## here we construct a factor to associate each sample to a level of
  ## the design
  ## we need to add an error if the nout is bigger than the minimum size of the classes
  
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



#' @title validate
#' @md
#' @description 
#' This function implements a resampling strategy (based on bootstrap or jacknifing) to validate the results 
#' of an ASCA decomposition in terms of loadings and variable importance
#' 
#' @details 
#' Validation is performed by evaluating the variability of the estimates of the variable importance and 
#' of the loadings upon a bootstrapping/jacknifing which preserve the experimental design. The results are 
#' presented in form of the empirical quantiles of the validated parameters.
#' 
#' @param x An asca object
#' @param method A character vector which specify the type of validation strategy. Either "jacknifing" (the default) or 
#' "bootstrap". 
#' @param turns How many independent resampling should be performed
#' @param nout How many samples will be "left" out in the case of jacknifing
#' @param qt A vector with the quantiles which will be calculated on the resamplings
#'
#' @return
#' an object of class `asca` containing also the results of the validation.
#' 
#' @examples
#' 
#' ## load the data
#' data("synth_count_data")
#' 
#' ## perform the ASCA decomposition
#' ASCA_test <- ASCA(
#' d = synth_count_data$design,
#' x = synth_count_data$counts, 
#' f = "time + treatment + time:treatment",
#' glm_par = list(family = poisson()))
#' 
#' ## perform validation
#' ASCA_test_v <- validate(ASCA_test,
#' turns = 100,
#' )
#' 
#' ## the results are included in the asca object inside the validation element
#' 
#' @export
#' 
validate <- function(x,method,turns,nout,qt) {
  UseMethod("validate")
}

#' @export
validate.asca <- function(x,
                          method = c("jacknifing","bootstrap"),
                          turns = 200, 
                          nout = 1,
                          qt = c(0.05,0.95)){
  ## extract the "reference " loadings from the decomposition of the complete matrix
  refload <- loadings(x)
  terms <- names(refload)
  method <- match.arg(method)
  cat("Validation: ",method, " turns: ", turns )
  ## perform the repeated validations
  turns <- lapply(1:turns, function(b){
    if (method == "jacknifing") {
      val <- do_jacknife(x$d,x$x,nout = nout)
    } else {
      val <- do_bootstrap(x$d,x$x)
    }
    
      asca_v <- ASCA(d = val$d_val,
                x = val$x_val,
                f = x$f,
                comb = x$comb,
                glm_par = x$glm_par)
    
    ## rotate the loadings if needed
    rotated_loadings <- lapply(terms, function(t) {
      out <- fix_sign(refload[[t]],
                      loadings(asca_v)[[t]])
      return(out)
    })
    
    ## we put in the asca objects the "rotated" loadings
    names(rotated_loadings) <- terms
    asca_v$svd$loadings <- rotated_loadings 
    
    return(list(loadings = loadings(asca_v),
                scores = scores(asca_v),
                varimp = varimp(asca_v)))
  })
  
  
  load_array <- lapply(terms, function(x) vapply(turns, function(t) t$loadings[[x]], turns[[1]]$loadings[[x]]))
  names(load_array) <- terms
  
  scores_array <- lapply(terms, function(x) vapply(turns, function(t) t$scores[[x]], turns[[1]]$scores[[x]]))
  names(scores_array) <- terms
  
  varimp_array <- vapply(turns, function(t) t$varimp, turns[[1]]$varimp)
  
  ## now we need to calculate the empitrical quantiles
  # # of the loadings, the scores and the variable importance
  loadings_q1 <- lapply(load_array, function(term) apply(term,c(1,2),function(x) quantile(x,qt[1])))
  loadings_q2 <- lapply(load_array, function(term) apply(term,c(1,2),function(x) quantile(x,qt[2])))
  
  scores_q1 <- lapply(scores_array, function(term) apply(term,c(1,2),function(x) quantile(x,qt[1])))
  scores_q2 <- lapply(scores_array, function(term) apply(term,c(1,2),function(x) quantile(x,qt[2])))
  
  varimp_q1 <- apply(varimp_array,c(1,2),function(x) quantile(x,qt[1]))
  varimp_q2 <- apply(varimp_array,c(1,2),function(x) quantile(x,qt[2]))
  
  ## put the outcomes in the validation list of the asca object
  
  qtnames <- paste0("q",qt)
  validation_list_names <- c("method", "varimp", "loadings","scores")
  x$validation <- vector("list", length(validation_list_names))
  names(x$validation) <- validation_list_names
  
  x$validation$method <- method
  x$validation$loadings <- list(loadings_q1,loadings_q2)
  names(x$validation$loadings) <- qtnames
  x$validation$varimp <- list(varimp_q1,varimp_q2)
  names(x$validation$varimp) <- qtnames
  x$validation$scores <- list(scores_q1,scores_q2)
  names(x$validation$scores) <- qtnames
  return(x)
} 

