#' Decompose the terms ASCA
#'
#'The function performs the svd decomposition of the different ASCA terms
#'
#' @param arr a 3d array containing the terms of the decomposition. See details.
#' @param eigthresh the threshold (relative importance) used to discard low importance components. The default value is 0.1 "\%" of the sum off the eigenvalues
#' @param returnloadings should the loadings of svd be returned
#'
#' @return A list with the following elements:
#' \itemize{
#'  \item{loadings: }{the loadings of the initial variables in each term}
#'  \item{varimp: }{the importance of each initial variable for the decomposition of each term}
#'  \item{variances: }{the variance cexplained by each decomposition term, calculated from the SVD eigenvalues}
#' }
#'
#' @details The input array is 3D. The first dimensions is the row of the initial matrix (samples).
#' The second dimension encodes the ASCA decomposition terms. The third one represents the columns of the intial
#' data matrix, i.e. the original variables.
#'
#' The explained variance is calculated as the squared sum of the "relevant" terms of the d vector containing the singular values obtained from svd
#'
#' @export
#'
#' @examples
#'

ASCA_svd <- function(arr,
                     eigthresh = 0.001,
                     returnloadings = TRUE){

    svs <- apply(arr,2, function(mat) {
      ## perform svd af the separate terms
      sv <- svd(scale(mat,scale = FALSE))
      idin <- which(sv$d/sum(sv$d) > eigthresh)
      variances <-  sv$d[idin]^2
      loadings <-  sv$v[,idin, drop = FALSE]

      ## calculate the variable importance combining loadings and explained variances
      varimp <- rowSums(sweep(abs(loadings),2,variances,"*")/sum(variances))
      if (returnloadings) {
        out <- list(loadings = loadings,
                    varimp = varimp,
                    variances = sum(variances))
      } else {
        out <- list(loadings = NULL,
                    varimp = varimp,
                    variances = sum(variances))
      }
      return(out)
      })

    ## add the term name to the output list
    names(svs) <- dimnames(arr)[[2]]

    ## format the output
    return(list(loadings = lapply(svs, function(x) x$loadings),
                varimp = sapply(svs, function(x) x$varimp),
                variances = t(sapply(svs, function(x) x$variances))
                )
           )
    }




#' Matrix norm for ASCA
#'
#' Calculate the squared norm of a matrix as the square of the Frobenius norm
#'
#' @param m a matrix of numbers
#'
#' @return the squared Frobenius norm
#'
#' @details The squared matrix norm is used to quantify the percentage of variation explained by each ASCA component
#'
#' @export
#'
#'
#'

ASCA_mnorm <- function(m){
  out <- sum(abs(m)^2)
  return(out)
}


#' ASCA decomposition
#'
#' Perform the ASCA decomposition of a data matrix according to a given experimental design
#'
#' @param d  a data.frame/matrix of factors specifying the design
#' @param x  a data.fram/matrix with the response (n samples x m variables)
#' @param f  the formula representing the decomposition
#'
#' @return A list with the following elements:
#' \itemize{
#'  \item{decomposition: }{a 3d decomposition array}
#'  \item{svds: }{the SVD of the ASCA terms}
#'  \item{exp_variation: }{the initial data frame }
#'  \item{d: }{the data.frame/matrix of factors specifying the design}
#'  \item{x: }{the data.fram/matrix with the response (n samples x m variables)}
#'  \item{f: }{the formula representing the decomposition}
#' }
#'
#' @details The decompusition is calculated by recursively applying lm/glm with the correct cotrasts (sum)
#'
#'
#'
#' @export
#'
#' @examples

ASCA_decompose <- function(d,x,f){

  ## contrasts should be changed so first I save the old ones
  oldcontrasts <- options("contrasts")
  ## then I set the new ones
  options(contrasts = rep ("contr.sum", 2))

  ## calculate the models for all variables
  ## return the terms calculated
  models <- mclapply(1:ncol(x), function(v){
    df <- cbind(y = x[,v],d)
    myform = paste0("y ","~",f)
    ## the models
    mod <- lm(formula = myform, data = df)
    ## the terms in frm of matrix
    terms <- predict(mod, type = "terms")
    ## the error of the decomposition
    err <- residuals(mod)
    return(cbind(terms, err = err))
  })

  ## organize thems in a 3d array
  termarray <- sapply(models, function(m) m, simplify = "array")

  dimnames(termarray) <- list(NULL,
                              colnames(models[[1]]),
                              colnames(x))

  ## perform the svd of all the terms excluding the error
  svds <- ASCA_svd(termarray[,1:dim(termarray)[2]-1,,drop = FALSE])

  ## reset the contrasts!
  options("contrasts" = oldcontrasts$contrasts)


  return(list(decomposition = termarray,
              svds = svds,
              exp_variation = apply(termarray,2,ASCA_mnorm),
              d = d,
              x = x,
              f = f))
}





#' ASCA_validate
#'
#' apply a permutation approach to validate the significance of the different terms of the ASCA decomposition
#'
#' @param asca the output of the ASCA_decompose
#' @param npermutations the number of permutations used to construct the null distribution
#'
#' @details The significance of each term of the ASCA decomposition
#' is evaluated by calculating the empirical distribution of the explained variance by recursive permutation of the
#' design factors. Each decomposition term is considered significance if it accounts for a fraction of variance
#' larger than the one that would be obtained by chance by a random assignment of the samples to the design factors
#'
#'
#' @return
#' @export
#'
#' @examples
#'

ASCA_validate <- function(asca,
                          npermutations = 100,
                          upquantile = 0.95){
  rnd_variances <- t(sapply(1:npermutations, function(n){
    ## permute the design factors
    rnd_d <- as.data.frame(lapply(asca$d,sample))
    rnd_asca <- ASCA_decompose(rnd_d,asca$x,asca$f)
    return(rnd_asca$exp_variation)
  }))

  return(list(random_variances = rnd_variances,
              quantiles = apply(rnd_variances, 2, function(x) quantile(x,upquantile))
              )
         )
}






## to do
## - x should be matrix, d data frame
## - evaluate parallel computing for validation
## - multivariate vs univariate testing?







