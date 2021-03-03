

#' @title plot.asca
#' @md
#' @description 
#' General method for plotting an `asca` object
#' 
#' @param x the `asca` object
#' @param type a character string specifying the type of the plot, see details
#' @param term the term of the decomposition
#' @param term2 optional. The second term to be displayed. Used only if the type is `scores`
#' @param eig the eigenvector of the term
#' @param eig2 optional. The eigenvector of `term2`
#' @param sort should de variable be ordered according to their importance?
#' @param validation should the plot include information on the validation? Used only if type is `loadings` or `varimp`
#' @param n.var number of variables to be displayed
#' @param main title of the plot
#' @param ... additional arguments passed to the base plotting functions
#'
#'
#' @details
#' The method allow for three type of plots. `varimp` plot shows the contribution of the initial variables 
#' to one of the terms of the decomposition. When the argument `validation` is true, the plot also shows 
#' the confidence intervals obtained during validation. `loadings` plot allows to visualize the weight 
#' of the initial variables on the eigenvectors of a specific decomposition term. Also here, the `validation` 
#' argument will include the confidence intervals. 
#' The last type of plot (`scores`) displays the position of the samples on a projection plane which 
#' potentially combines the different decomposition terms  
#'
#' @export
#'
#' @examples
#' 
#' ## load the data
#' data("synth_count_data")
#' 
#' ## perform the ASCA decomposition
#' ASCA_test <- ASCA(
#'  d = synth_count_data$design,
#'  x = synth_count_data$counts, 
#'  f = "time + treatment + time:treatment",
#'  glm_par = list(family = poisson())
#' )
#'
#' ## plot the variable importance for the first eigenvector of the "time" term
#' plot(ASCA_test, type = "varimp", term = "time", eig = 1)
#'   
#' ## for the same term, the loadings are instead ...
#' plot(ASCA_test, type = "loadings", term = "time", eig = 1)
#' 
#' ## decomposition with a combination of treatment and interaction
#' ASCA_test_comb <- ASCA(
#'  d = synth_count_data$design,
#'  x = synth_count_data$counts, 
#'   f = "time + treatment + time:treatment",
#'   comb = c("treatment","time:treatment"),
#'   glm_par = list(family = poisson())
#' )
#' 
#' 
#' ## plot the scores of eig1 for "time" and teig1 "treatment+time:treatment"
#' mycol <- c("#081d5860", "#225ea860", "#1d91c060","#7fcdbb60", "#c7e9b460")
#' names(mycol) <- unique(ASCA_test_comb$d$time)
#' mypch <- c(15,1,19)
#' names(mypch) <- unique(ASCA_test_comb$d$treatment)
#' plot(ASCA_test_comb, 
#'     type = "scores", 
#'     term = "time", 
#'     eig = 1, 
#'     term2 = "treatment+time:treatment", 
#'     eig2 =  1,
#'     col = mycol[ASCA_test_comb$d$time],
#'     pch = mypch[ASCA_test_comb$d$treatment],
#'     cex = 1)
#' 
#' ## perform validation
#' ASCA_test_v <- validate(ASCA_test,
#' turns = 100)
#' 
#' ## plot the variable importance with the confidence limits
#' plot(ASCA_test_v, type = "varimp", term = "time", eig = 1, validation = TRUE)
#' 
 
plot.asca <- function(x,
                      term = 1,
                      type = c("varimp","loadings","scores"),
                      term2 = term,
                      eig = 1, 
                      eig2 = 2,
                      sort = NULL,
                      validation = FALSE,
                      n.var = min(30, ncol(x$x)),
                      main = NULL,
                      ...
                      ){
  type <- match.arg(type)
  switch(type,
         varimp =
           plot_varimp(asca = x, term = term, eig = eig, sort = sort, 
                       validation = validation, n.var = n.var, main = main, 
                       ...),
         loadings = 
           plot_loadings(asca = x, term = term, eig = eig, sort = sort, 
                         validation = validation, n.var = n.var, main = main, 
                         ...),
         scores = 
           plot_scores(asca = x, 
                       term1 = term,
                       term2 = term2,
                       eigen1 = eig,
                       eigen2 = eig2, 
                       main = main,
                       ...)
  )
  
}


## Here a generic function to plot ASCA
plot_varimp <- function(asca, term, 
                        eig,
                        sort, 
                        validation,
                        n.var,
                        main,
                        ...){
  
  if (is.null(main)) main <- paste0("Factor: ",term)
  if (is.null(sort)) sort <- TRUE
  ## get the variable importance and find the order for having them descending
  imp <- varimp(asca)[,term]
  ord <- if (sort) rev(order(imp,decreasing=TRUE)[1:n.var]) else 1:n.var
  
  ## plot without validation intervals
  if (!validation) {
    dotchart(imp[ord], xlab= "Variable Importance", ylab="", main=main, ...)
    segments(x0 = 0, x1 = imp[ord], y0 = seq(1:n.var), y1 = seq(1:n.var))
    abline(v = 0, lty = 2, lwd = 1)
  } else {
  ## plot with the validation intervals  
    if (is.null(asca$validation)) {stop("The asca object has not been validated")} 
  ## get out the quantiles of the validations
    q1 <- asca$validation$varimp[[1]][,term]
    q2 <- asca$validation$varimp[[2]][,term]
  ## make the dotchart plot with the quantiles
    dotchart(imp[ord], xlab= "Variable Importance", ylab="", main=main, 
             xlim = c(0,max(q2[ord])), 
             ...)
    points(y = seq(1:n.var), x = q2[ord], pch= "|", col = "#8a3335")
    points(y = seq(1:n.var), x = q1[ord], pch= "|", col = "#8a3335")
    segments(x0 = q1[ord], x1 = q2[ord], y0 = seq(1:n.var), y1 = seq(1:n.var), col = "#8a3335")
    abline(v = 0, lty = 2, lwd = 1)
  }
}


plot_loadings <- function(asca, 
                          term, 
                          eig,
                          sort, 
                          validation,
                          n.var,
                          main, ...){
  
  if (is.null(main)) main <- paste0("Factor: ",term,", Eigenvector: ",eig)
  if (is.null(sort)) sort <- FALSE
  
  
  imp <- loadings(asca)[[term]][,eig]
  ord <- if (sort) rev(order(imp,decreasing=TRUE)[1:n.var]) else 1:n.var
  
  if (!validation) {
    dotchart(imp[ord], xlab= "Variable Loading", ylab="", main=main, ...)
    segments(x0 = 0, x1 = imp[ord], y0 = seq(1:n.var), y1 = seq(1:n.var))
    abline(v = 0, lty = 2, lwd = 1)
    
  } else {
    if (is.null(asca$validation)) {stop("The asca object has not been validated")} 
    q1 <- asca$validation$loadings[[1]][[term]][,eig]
    q2 <- asca$validation$loadings[[2]][[term]][,eig]
    dotchart(imp[ord], xlab= "Variable Loading", ylab="", main=main, 
             xlim = range(c(q1,q2)), 
             ...)
    points(y = seq(1:n.var), x = q2[ord], pch= "|", col = "#8a3335")
    points(y = seq(1:n.var), x = q1[ord], pch= "|", col = "#8a3335")
    segments(x0 = q1[ord], x1 = q2[ord], y0 = seq(1:n.var), y1 = seq(1:n.var), col = "#8a3335")
    abline(v = 0, lty = 2, lwd = 1)
  }
}


plot_scores <- function(asca, 
                        term1,
                        term2,
                        eigen1,
                        eigen2,
                        main,
                        ...) {
  
  if (is.null(main)) main <- paste0("Factors: ",paste0(term1, " vs ", term2),"\n",
                                    "Eigenvectors: ", eigen1,"vs",eigen2)
  
  plot(x = asca$svd$scores[[term1]][,eigen1],
       y = asca$svd$scores[[term2]][,eigen2], 
       xlab = paste0("F ",term1, ", E", eigen1), 
       ylab = paste0("F ",term2, ", E", eigen2),
       main = main,
       type = "p",
       ...)
}
  

