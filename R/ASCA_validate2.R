# 
#' ASCA validation2
#'
#'Function to validate the ASCA decomposition; apply a permutation approach to validate the significance of the different terms of the ASCA decomposition
#'
#'
#' @param asca The object retured by the ASCA_decompose
#' @param npermutations Number of permutations to be used to construct the null distribution of ASCA decomposition terms
#' @param upquantile The upper quantile limits of the null distribution of ASCA terms, used for testing significance
#'
#' @return Returns a list with the following elements:
#' \itemize{
#'  \item{random_variances: }{a matrix with null values of explained variance by each factor}
#'  \item{upper.quantile.variances: }{the upper quantile of the null distribution of explained variance}
#'  \item{random_varimp.array: }{an array of null values of variable importance; the last dimension of the array equals npermutation}
#'  \item{upper.quantile.varimp: }{the observed variable importance, reporting only variable-factor cases where values > null upperquantile}
#'  }
#' 
#' 
#' @details The significance of each term of the ASCA decomposition
#' is evaluated by calculating the empirical distribution of the explained variance and variable importance 
#' by recursive permutation of the design factors. Each decomposition term is considered significance if it accounts for a fraction of variance
#' larger than the one that would be obtained by chance by a random assignment of the samples to the design factors
#' 
#' @export
#'
#' @examples
#' 
#' 
ASCA_validate2<-function(asca,
                         npermutations=100,
                         upquantile=0.95){
  
  shf_varimp <- function(asca){  #internal function to create a shuffle-based varimp
    
    rnd_d <- as.data.frame(lapply(asca$d,sample)) # shuffle the design factors
    rnd_asca <- ASCA_decompose(rnd_d,asca$x,asca$f) # run ASCA on the randomised data structure
    rnd_varimp<-rnd_asca$svds$varimp # the randomised varimp (variable importance)
    rnd_variation<-rnd_asca$exp_variation # the randomised explained variation
    return(list(rnd_variation=rnd_variation, rnd_varimp=rnd_varimp))
  }
  
  rnd_asca.list<- replicate(npermutations, shf_varimp(asca)) # replicate the shuffle ntimes; produces array of random variation + random varimp
  
  rnd_varimp.array<- # isolate the array of random varimp and assing names
    array(unlist(rnd_asca.list[2,]), dim=c(dim(asca$svds$varimp),npermutations), dimnames=list(NULL,dimnames(asca$svds$varimp)[[2]],NULL))
  
  upperquantile.varimp<- apply(rnd_varimp.array, c(1,2), function(x) quantile(x,upquantile)) # get the upper quant of the varimp table
  obs.varimp<-asca$svds$varimp # extract the observed varimp from the asca decomposition object
  obs.varimp[obs.varimp < upperquantile.varimp] <-NA # only show the variable-factor cases where obs varimp > upperquantile random
  
  rnd_variances.mat<- # isolate the array of random variation and assign names
    t(matrix(unlist(rnd_asca.list[1,]),nrow=length(asca$exp_variation)))
  dimnames(rnd_variances.mat)<-list(NULL, names(asca$exp_variation))
  
  # return all
  return=list(random_variances=rnd_variances.mat, 
              upperquantile.variances=apply(rnd_variances.mat, 2, function(x) quantile(x,upquantile)),
              random_varimp.array=rnd_varimp.array,
              upperquantile.varimp= upperquantile.varimp,
              important.var=obs.varimp)
}




#' Plot ASCA validate
#'
#' @param asca The object retured by the ASCA_decompose
#' @param asca_validate The object returned by the ASCA_validate2
#' @param var The variable that needs to plot (called by number)
#' @param fact The factor for which the variable importance is calculated (called by number)
#'
#' @return
#' @export
#'
#' @details The plot is used to visually explore the deviation of the observed variable importance (for a given factor) relative to its null distribution.
#' Each variable and factor should be called in the function
#' 
#' @examples
#' 
#' 
#' 
plot.validate.var<-function(asca, asca_validate, var, fact){
  hist(asca_validate$random_varimp.array[var,fact,], breaks = 30, xlim=c(0,1),
       main=substitute(paste("Histogram of ", asca_validate, sep=", ", "variable ", var, sep="_", "factor ", fact)))
  abline(v=asca$svds$varimp[var,fact], col="red")
}


