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
#' Cardinale, M., and F. Arrhenius. 2001. Is the decrease in growth rate of Atlantic Herring in the Baltic Sea Density-Dependent? A Geostatistical application. Herring: Expectations for aNew Millenium. Alaska Sea Grant College Program AK-SG-01-04:153-197.
#' Kurbanov, A.R., and B. Kamilov. 2015. Age and growth of bighead carp (Hypophthalmichthys nobilis R.) in Tudakul reservoir, Uzbekistan. Turkish Journal of Fisheries and Aquatic Sciences 3:229-232.
#' Michaletz, P.H., D.M. Nicks, and E.W. Buckner Jr. 2009. Accuracy and precision of estimates of back-calculated channel catfish lengths and growth increments using pectoral spines and otoliths. North American Journal of Fisheries Management 29:1664-1675.
#' Munro, A.R. 2004. Identifciation of life history variation in salmonids using otolith microchemistry and scale patterns: implications for illegal introductions and for whirling diesease in Missouri River Rainbow Trout. Dissertation. Montana State University. pp219.
#' Zalachowski, W., and K. Wieski. 1998. Growth rate of Bream in Lake Dabie. Journal of Polish Agricultural Universities: http://www.ejpau.media.pl/volume1/issue1/fisheries/art-03.html
#' @return Simulated length data
#' @export
simannuliR<-function(ages,Linf.mu,Linf.cv,k.mu,k.cv,t0.mu,t0.cv,seed){
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

  #Age data from simageData function
  ageC=ages
  yearC<-rep(format(as.Date(Sys.Date(),"%Y/%m/%d"),"%Y"),length(ageC))
  Ninds<-length(ageC)
  Nage<-max(ageC)

  #Matrix to store radius measurements of ages a of the ith fish
  Lprev=matrix(NA,nrow=Ninds,ncol=Nage)

  #von B parameter estimates to simulate annual growth rate of the ith fish
  make.dummy.inds<-data.frame(
    Linf=Linf.mu*rlnorm(Ninds,0,Linf.cv),
    k=k.mu*rlnorm(Ninds,0,k.cv),
    t0=t0.mu*rlnorm(Ninds,0,t0.cv)
  )

  #Time step
  t=1:Nage

  #Individual-based von B model
  for(i in 1:Ninds){
    for(j in 1:Nage){
      Lprev[i,j]<-make.dummy.inds$Linf[i]*(1-exp(-make.dummy.inds$k[i]*(t[j]-make.dummy.inds$t0[i])))
    }
  }

  #Assign NA to all measurements greater than age of capture
  Lprev[col(Lprev) > ageC] <- NA
  namescol=c(paste0("age", 1:Nage))
  colnames(Lprev)<-namescol

  #Identify length at capture of the ith individual
  lengthC<-apply(Lprev,1,max,na.rm=TRUE)
  lengthC<-round(lengthC)
  
  #length-otolith regression parameters attained from:
  #Zalachowski et al. 1998; Cardinale et al. 2001; Munro 2004; Michaletz et al. 2009; Kurbanov et al. 2015

  beta0.obs<-c(-0.15,-0.146,-1.480,-0.19,-1.48,-0.26,-1.620,0.00,0.40,-0.022,-0.017,-0.006)
  beta1.obs<-c(0.0054,0.06,0.0230,0.0140,0.0230,0.0190,0.0290,0.09,0.006,0.006,0.008,0.004)
  eta.obs<-rnorm(lengthC,0,0.2)

  beta0.mu<-mean(rnorm(Ninds,mean(beta0.obs),0.1))
  beta1.mu<-mean(rnorm(Ninds,mean(beta1.obs),0.1))

  #Generate asymptotic length of the otolith of the ith fish of age a at capture
  radC<-beta0.mu+beta1.mu*lengthC+eta.obs

  #Store data in data.frame
  dat.new<-data.frame(id=1:Ninds,yearC,lengthC,ageC,radC)

  #Create matrix to store measurement of radius of the ith fish of age a
  annu.mat<-matrix(NA,nrow=Ninds,ncol=Nage)

  #Dahl-Lea model to estimate radius measurements based on the length of the ith fish of age a
  for(j in 1:Nage){
    annu.mat[,j]<-radC*(Lprev[,j]/lengthC)
  }

  #Assign NA to all measurements greater than age of capture
  annu.mat[col(annu.mat) > dat.new$ageC] <- NA
  namescol=c(paste0("annu", 1:Nage))
  colnames(annu.mat)<-namescol

  datR<-cbind(dat.new,annu.mat)
  
  dat<-gather(datR,ageR,annuR,annu1:colnames(datR[ncol(datR)]),factor_key=TRUE)%>%arrange(id,ageR)
  str_sub(dat$ageR,start=1,end=4)<-""
  dat%<>%mutate(ageR=as.numeric(ageR))%>%
    filter(!is.na(annuR))%>%
    filter(ageR<=ageC)
 
  #Identify year in which a fish in year class will be of age a
  dat$year=as.numeric(as.character(dat$yearC))-(dat$ageC-dat$ageR)
  return(dat)
}
