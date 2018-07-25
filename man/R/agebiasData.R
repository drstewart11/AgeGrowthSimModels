#' Simulates age:bias data
#'
#' Simulates biased age data from the simulated age sample
#' @param ages = simulated age data, bias = type of bias ("under","over","unbiased")
#' @return Age data
#' @export
agebiasData<-function(ages,bias=c("under","over","unbiased"),seed){
  if (missing(bias)) stop("'bias'must be specified as 'under','over',or 'unbiased'",call.=FALSE)
  if (length(bias)>1|missing(bias)) stop("'bias'must contain only one value",call.=FALSE)

  if(missing(seed)){
    seed=123
  }
  set.seed(seed)

  if(bias=="under"){
    if(max(ages)<=14){
      age2=log(ages)*1/1.06
      age2=abs(round(rnorm(length(ages),mean=exp(age2),sd=1)))
    }else if(max(ages)>=15 && max(ages)<=30){
      age2=log(ages)*1/1.08
      age2=abs(round(rnorm(length(ages),mean=exp(age2),sd=1)))
    }else if(max(ages)>=31){
      age2=log(ages)*1/1.10
      age2=abs(round(rnorm(length(ages),mean=exp(age2),sd=1)))
    }
  }else if(bias=="over"){
    if(max(ages)<=14){
      age2=log(ages)*1.20
      age2=abs(round(rnorm(length(ages),mean=exp(age2),sd=1)))
    }else if(max(ages)>=15 && max(ages)<=30){
      age2=log(ages)*1.07
      age2=abs(round(rnorm(length(ages),mean=exp(age2),sd=1)))
    }else if(max(ages)>=31){
      age2=log(ages)*1.05
      age2=abs(round(rnorm(length(ages),mean=exp(age2),sd=1)))
    }
  }else if(bias=="unbiased"){
    if(max(ages)<=14){
      age2=log(ages)*1
      age2=abs(round(rnorm(length(ages),mean=exp(age2),sd=2)))
    }else if(max(ages)>=15 && max(ages)<=30){
      age2=log(ages)*1
      age2=abs(round(rnorm(length(ages),mean=exp(age2),sd=2)))
    }else if(max(ages)>=31){
      age2=log(ages)*1
      age2=abs(round(rnorm(length(ages),mean=exp(age2),sd=2)))
    }
  }

  return(age2)

}
