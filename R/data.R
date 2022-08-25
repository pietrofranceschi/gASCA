#' 
#' 
#' Synthetic count data
#' @md
#' @description 
#' A dataset mimicking the results of a counting experiments of 15 variables over 150 samples.
#' 15 count variables are obtained as a function of 3 different treatments at 5 points in time. 
#' These variables could represent, for instance, count data from 10 different invertebrate species found on plants at 
#' Time (A)  0, 1, 2, 3, 4 after the plants have been treated with one of three treatments (B1, B2, B3). 
#' For each combination of time and treatment, K = 10 replicates (k = 1...K) were simulated. 
#' 
#' In particular.
#' * Two principal components will be used to describe the variation in matrix H_A (time), he variation described by the first component will be only present in the first two variables, while the variation of the second component will be only present in variables 3 and 4.
#' * One principal component will be used to describe the treatment effect in matrix $H_{B}$ (treatment) and this variation was only present in variables 5 and 6. 
#' * Two components will be used for the interaction effect in matrix HAB: PC1 was only active in variables 7 and 8, while PC2 variation, which was only assigned to found in variables 9 and 10.
#' * Five additional variables not responding to the design were added to the dataset
#' 
#'
#' @format a list with two components:
#' \describe{
#'   \item{counts}{a matrix with the counts (150x15)}
#'   \item{design}{a data.frame with the design factors}
#'   ...
#' }

"synth_count_data"