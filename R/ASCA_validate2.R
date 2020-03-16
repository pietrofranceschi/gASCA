

ASCA_validate2<-function(asca,
                         npermutations=100,
                         upquantile=0.95){
  
  shf_varimp <- function(asca){  #internal function to create a shuffle-based varimp
    
    rnd_d <- as.data.frame(lapply(asca$d,sample))
    rnd_asca <- ASCA_decompose(rnd_d,asca$x,asca$f)
    rnd_varimp<-rnd_asca$svds$varimp # the randomised varimp
    rnd_variation<-rnd_asca$exp_variation 
    return(list(rnd_variation=rnd_variation, rnd_varimp=rnd_varimp))
  }
  
  rnd_asca.list<- replicate(npermutations, shf_varimp(asca)) # replicate the shuffle ntimes -> array of random variances + random varimp
  
  rnd_varimp.array<- # isolate the array of random varimp and assign names
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



###########

plot.validate.varimp<-function(asca, asca_validate, var, fact){
  hist(asca_validate$random_varimp.array[var,fact,], breaks = 30, xlim=c(0,1),
       main=substitute(paste("Histogram of ", asca_validate, sep=", ", "variable ", var, sep="_", "factor ", fact)))
  abline(v=asca$svds$varimp[var,fact], col="red")
}


