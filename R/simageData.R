#' Simulates age data
#'
#' User-specified age-structured population model to generate ages from fishery-dependent or independent data such that samples would vary in age composition and abundance
#' @param A = max number of ages classes, M = Average Natural Mortality, FM = Average Fishing Mortality, f.vul = Age first vulnerable to fishery, sd.fvul = sd of age vulnerable to fishery, Fa = Fishing mortality rate for fish of age a, FM = Fishing mortality rate, M = Average natural mortality rate (per year), Na = Total numbers of fish of age a, Ca = The number of fish of age a captured, Popsize = Population size, fishery = The number of fish of age a retained for age sample from fishery (TRUE) and survey (FALSE)
#' @return Simulated age data
#' @export
simageData<-function(A,M,FM,f.vul,sd.fvul,s.vul,sd.svul,Popsize,seed,fishery=TRUE){
  #Error bounds
  if (length(A)>1|missing(A)) stop("'A'must contain only one value",call.=FALSE)
  if (A<1) stop(" 'A' must be a positive number",call.=FALSE)
  if (length(M)>1|missing(M)) stop(" 'M' must contain only one value", call.=FALSE)
  if (M<0) stop(" 'M' must be a positive number",call.=FALSE)
  if (length(FM)>1|missing(FM)) stop(" 'FM' must contain only one value",call.=FALSE)
  if (FM<0) stop(" 'FM' must be a positive number",call.=FALSE)
  if (length(f.vul)>1|missing(f.vul)) stop(" 'f.vul' must contain only one value",call.=FALSE)
  if (f.vul<1) stop(" 'f.vul' must be a positive number and greater than 1",call.=FALSE)
  if (length(sd.fvul)>1|missing(sd.fvul)) stop(" 'sd.fvul' must contain only one value",call.=FALSE)
  if (sd.fvul<0) stop(" 'sd.fvul' must be a positive number", call.=FALSE)
  if (length(Popsize)>1|missing(Popsize)) stop(" 'Popsize' must contain only one value",call.=FALSE)
  if (Popsize<0) stop(" 'Popsize' must be a positive number",call.=FALSE)
  if (A<f.vul) stop(" 'Number of age classes must be greater than F.vul'",call.=FALSE)


  if(missing(seed)){
    seed=123
  }
  set.seed(seed)


  #Number of age classes from 1 to max(age)
  age=1:A

  #Fishery selectivity: Selectivity for the fishery is a logistic function
  va = (1+exp(-(age-f.vul)/sd.fvul))^-1

  #Fishing mortality rate for fish of age a
  Fa = FM*va

  #Mortality (both fishing and natural mortality) rate for fish of age a
  Za = Fa + M

  #Age structured population model; The total number of fish of age a
  Na=Popsize*exp(-Za*age)

  #Number of fish of age a captured using fishing gear
  Ca = Na*(1-exp(-Za*age))*Fa/Za

  #Survey selectivity
  sa = (1+exp(-(age-s.vul)/sd.svul))^-1

  #Survey catchability
  q.surv = 0.05*exp(-0.5*0.3)

  #Survey abundance
  la = q.surv*sa*Na

  if(fishery){
    age_sample = sample(age,sum(Ca),replace=TRUE,prob=Ca/sum(Ca))
  }else{
    age_sample = sample(age,sum(la),replace=TRUE,prob=la/sum(la))
  }

  par(mfrow=c(3,1),mar=c(3,5,0.5,0.5),oma=c(.2,.2,.2,.2),mgp=c(2,0.5,0),cex=1,tck=-0.02)
  plot(age,Na,"h",col="red",lwd=10,ylab="",xlab="",las=1,bty="n",xlim=c(0,A),cex.axis=1.5,cex.lab=1.5)
  mtext(side=2,"No.in the pop.",cex=1.5,line=3.5)
  plot(age,Ca,"h",col="blue",lwd=10,ylab="",xlab="",las=1,bty="n",xlim=c(0,A),cex.axis=1.5,cex.lab=1.5)
  mtext(side=2,"No. in fishery",cex=1.5,line=3.5)
  plot(age,la,"h",col="black",lwd=10,ylab="",xlab="",las=1,bty="n",xlim=c(0,A),cex.axis=1.5,cex.lab=1.5)
  mtext(side=2,"No. in survey",cex=1.5,line=3.5)
  mtext(side=1,"Age (yr)",cex=1.5,line=1.5)
  ages=age_sample
  return(ages)
}
