#' Synthetic count data
#'
#' A dataset mimicking the results of a counting experiments of 10 variables over 45 samples.
#' 
#' Design: 
#' 
#' * 5 timepoints; 
#' * 3 treatments ;
#' * 3 replicates
#' 
#' The variables were created from a two latent variable model which were "projected" on the the final variables 
#'
#' @format a list with two components:
#' \describe{
#'   \item{counts}{a matrix with the counts (45x10)}
#'   \item{design}{a data.frame with the designh factors}
#'   ...
#' }

"synth_count_data"