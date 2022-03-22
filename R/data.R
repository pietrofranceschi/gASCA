#' 
#' 
#' Synthetic count data
#' @md
#' @description 
#' A dataset mimicking the results of a counting experiments of 10 variables over 45 samples.
#' 
#' Design: 
#' 
#' * 5 timepoints; 
#' * 3 treatments ;
#' * 10 replicates
#' 
#' The variables were created from a two latent variable model which were "projected" on the the final variables. In details
#' 
#' * s1,s2,s3,s4 respond to the "time" factor. s1 and s2 contribute to the first latent factor, s3 and s4 to the second
#' * s5,s6  respond to the "treatment" factor
#' * s1 and s2 respond also to the interaction term "time:treatment". In other words, for them the effect of treatment depends on time.
#'
#' @format a list with two components:
#' \describe{
#'   \item{counts}{a matrix with the counts (45x10)}
#'   \item{design}{a data.frame with the designh factors}
#'   ...
#' }

"synth_count_data"