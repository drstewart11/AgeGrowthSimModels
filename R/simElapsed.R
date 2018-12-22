#' Individual-based model to simulate mean length at annulus for the ith fish of age a
#'
#' Simulates data to estimate back-calculated lengths of age a for the ith fish
#' @param ages = ages simulated from simageData function
#' @param Linf.mu = Linf mean
#' @param Linf.cv = CV of Linf
#' @param k.mu = k mean
#' @param k.cv = CV of k
#' @param t0.mu = t0 mean
#' @param t0.cv = CV of t0
#' @param seed = seed
#' @references
#' @return Simulated length data
#' @export
simElapsed<-function(ages,Linf.mu,Linf.cv,k.mu,k.cv,t0.mu,t0.cv,seed){
  #Error bounds
  if (length(Linf.mu)>1|missing(Linf.mu)) stop("'Linf.mu'must contain only one value",call.=FALSE)
  if (Linf.mu<1) stop(" 'Linf.mu' must be a positive number",call.=FALSE)
  if (length(Linf.cv)>1|missing(Linf.cv)) stop("'Linf.cv'must contain only one value",call.=FALSE)
  if (Linf.cv<0) stop(" 'Linf.cv' must be a positive number",call.=FALSE)
  if (length(k.mu)>1|missing(k.mu)) stop(" 'k.mu' must contain only one value", call.=FALSE)
  if (k.mu<0) stop(" 'k.mu' must be a positive number",call.=FALSE)
  if (length(k.cv)>1|missing(k.cv)) stop("'k.cv' must contain only one value",call.=FALSE)
  if (k.cv<0) stop(" 'k.cv' must be a positive number",call.=FALSE)
  if (length(t0.mu)>1|missing(t0.mu)) stop(" 't0.mu' must contain only one value",call.=FALSE)
  if (length(t0.cv)>1|missing(t0.cv)) stop("'t0.cv' must contain only one value",call.=FALSE)
  if (t0.cv<0) stop(" 't0.cv' must be a positive number",call.=FALSE)
  if(missing(seed)){
    seed=123
  }
  set.seed(seed)

  Ninds=length(ages)

  #Create dummy params for each individual of fish age a
  make.dummy.inds<-data.frame(
    Linf=Linf.mu*rlnorm(Ninds,0,Linf.cv),
    k=k.mu*rlnorm(Ninds,0,k.cv),
    t0=t0.mu*rlnorm(Ninds,0,t0.cv),
    #randomly extract individual ages from the age structured population model
    age.inds=ages,
    deltaT=rgamma(Ninds,shape=0.7,scale=2)
  )

  #Create empty vector to store initial predicted length for the ith individual of age a using von Bertalanffy function
  lt=rep(0,Ninds)

  #Individual-based von Bertalanffy growth model to predicted initial lengths for individuals of age a
  for(i in 1:Ninds){
    lt[i]<-make.dummy.inds$Linf[i]*(1-exp(-make.dummy.inds$k[i]*(make.dummy.inds$age.inds[i]-make.dummy.inds$t0[i])))
  }

  #Empty arrray's to store seasonal changes in length
  lengthT=array(0,dim=c(Ninds,nperiod=2))

  #Initial lengths
  lengthT[,1]<-lt

  #Difference
  deltaLt=(make.dummy.inds$Linf-lengthT[,1])*(1-exp(-make.dummy.inds$k*make.dummy.inds$deltaT))
  lengthT[,2]<-lengthT[,1]+deltaLt

  lengthT=cbind(lengthT,make.dummy.inds$deltaT)
  return(lengthT)
}
