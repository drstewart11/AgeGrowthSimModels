#' Individual-based model to simulate age:length and length:frequency data from a von Bertalanffy growth function
#'
#' Simulates length frequency data
#' @param A,Linf.mu,Linf.cv,k,k.cv,t0.mu,t0.cv,Popsize,FM,M,ts,C,bin.num = vonB growth params, Popsize=Population size,FM and M = Fishing and Nat Mortality, ts and C = seasonalized growth parameters, bin.num= The number of bins for histogram, seed=seed
#' @return Simulated length frequency data
#' @references
#' Somers, I.F. 1988. On a seasonally oscillating growth function. Fishbyte 6:8-11.
#' @export
simlengthFreq=function(A,Linf.mu,Linf.cv,k.mu,k.cv,t0.mu,t0.cv,Popsize,FM,M,ts.mu,ts.cv,C.mu,C.cv,bin.num,type=c("monthly","annual"),seed){
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
  if (length(FM)>1|missing(FM)) stop(" 'FM' must contain only one value",call.=FALSE)
  if (FM<0) stop(" 'FM' must be a positive number",call.=FALSE)
  if (length(M)>1|missing(M)) stop(" 'M' must contain only one value", call.=FALSE)
  if (Popsize<0) stop(" 'Popsize' must be a positive number",call.=FALSE)
  if (M<0) stop(" 'M' must be a positive number",call.=FALSE)
  if (length(ts.mu)>1|missing(ts.mu)) stop(" 'ts.mu' must contain only one value", call.=FALSE)
  if (length(ts.cv)>1|missing(ts.cv)) stop(" 'ts.cv' must contain only one value", call.=FALSE)
  if (length(C.mu)>1|missing(C.mu)) stop(" 'C.mu' must contain only one value", call.=FALSE)
  if (length(C.cv)>1|missing(C.cv)) stop(" 'C.cv' must contain only one value", call.=FALSE)
  if (length(bin.num)>1|missing(bin.num)) stop(" 'bin.num' must contain only one value, call.=FALSE")
  if(missing(seed)){
    seed=123
  }

set.seed(seed)

#Number of age classes from 1 to max(age)
ages=1:A

#Mortality (both fishing and natural mortality) rate for fish of age a
Za=FM+M

#Age structured population model; The total number of fish of age a
Na = round(Popsize * exp(-Za * ages))
  
# Normalize Na so the sum equals Popsize
Na_normalized = round(Na * Popsize / sum(Na))

# Now the sum should be equal to Popsize
Ninds = sum(Na_normalized)
Na = Na_normalized



#Create dummy params for each individual of fish age a
make.dummy.inds<-data.frame(
  Linf=Linf.mu*rlnorm(Ninds,0,Linf.cv),
  k=k.mu*rlnorm(Ninds,0,k.cv),
  t0=t0.mu*rlnorm(Ninds,0,t0.cv),
  #randomly extract individual ages from the age structured population model
  age.inds=sample(ages,Ninds,replace=TRUE,prob=Na/Ninds)
)


#Create empty vector to store initial predicted length for the ith individual of age a using von Bertalanffy function
lt=rep(0,Ninds)

#Individual-based von Bertalanffy growth model to predicted initial lengths for individuals of age a
for(i in 1:Ninds){
  lt[i]<-make.dummy.inds$Linf[i]*(1-exp(-make.dummy.inds$k[i]*(make.dummy.inds$age.inds[i]-make.dummy.inds$t0[i])))
}

#Parameters for seasonal growth model
lfq.inds<-data.frame(
  ts=ts.mu*rlnorm(Ninds,0,ts.cv),
  C=C.mu*rlnorm(Ninds,0,C.cv)
)

if(type=="monthly"){
#Empty arrray's to store seasonal changes in length
deltaT=array(0,dim=c(Ninds,13))
lengthT=array(0,dim=c(Ninds,13))

for(t in 1:1){
  #monthly time increment;
  t1 = (t - 1) * t.incr
  t2 = t * t.incr

  lengthT[,1] <- lt

  deltaT[,1] <- (make.dummy.inds$Linf - lengthT[,1]) * 
                (1 - exp(-(
                  make.dummy.inds$k * t.incr
                  - ((lfq.inds$C * make.dummy.inds$k) / (2 * pi)) * sin(2 * pi * (t1 - lfq.inds$ts))
                  + ((lfq.inds$C * make.dummy.inds$k) / (2 * pi)) * sin(2 * pi * (t2 - lfq.inds$ts))
                )))
}
#Second loop to model delta length across months (13)
for(t in 2:13){
  t1 = (t - 1) * t.incr
  t2 = t * t.incr
  
  deltaT[,t] <- (make.dummy.inds$Linf - lengthT[,t-1]) * 
                (1 - exp(-(
                  make.dummy.inds$k * t.incr
                  - ((lfq.inds$C * make.dummy.inds$k) / (2 * pi)) * sin(2 * pi * (t1 - lfq.inds$ts))
                  + ((lfq.inds$C * make.dummy.inds$k) / (2 * pi)) * sin(2 * pi * (t2 - lfq.inds$ts))
                )))
  lengthT[,t] <- lengthT[,t-1] + deltaT[,t]
}
#Format data to construct histogram through time

length.range<-range(unlist(lengthT))
lengthClass<-seq(floor(length.range[1]),ceiling(length.range[2])+1)
lengthMids <- lengthClass[-length(lengthClass)] + bin.num/2

#Histogram of counts for each age class through time t=1:13
freq=apply(lengthT,2,FUN=function(x){
  hist(x,breaks=lengthClass,plot=FALSE,include.lowest=TRUE)$counts})


timemin=as.Date("2018-01-01");timemax=as.Date("2019-01-01");dates=seq(timemin,timemax,'month')
dates<-as.Date(dates,"%Y-%d-%m")

#Create list object to store dates, midLengths, and histogram frequencies to use in TropFishR
lfq.dat<-list(dates=dates,
              midLengths=lengthMids,
              catch=as.matrix(freq))
return(lfq.dat)
}else if(type=="annual"){
  lengthT<-matrix(NA,nrow=length(lt),ncol=10)
  for(i in 1:length(lt)){
    lengthT[,]<-lt[]
  }
  length.range<-range(unlist(lengthT))
  lengthClass<-seq(floor(length.range[1]),ceiling(length.range[2])+1)
  lengthMids <- lengthClass[-length(lengthClass)] + bin.num/2

  #Histogram of counts for each age class through time t=1:13
  freq=apply(lengthT,2,FUN=function(x){
    hist(x,breaks=lengthClass,plot=FALSE,include.lowest=TRUE)$counts})

  timemin=as.Date("2018-01-01");timemax=as.Date("2027-01-01");dates=seq(timemin,timemax,'year')
  dates<-as.Date(dates,"%Y-%d-%m")

  #Create list object to store dates, midLengths, and histogram frequencies to use in TropFishR
  lfq.dat<-list(dates=dates,
                midLengths=lengthMids,
                catch=as.matrix(freq))
  return(lfq.dat)
}

#Format data to construct histogram through time
#length.range<-range(unlist(lengthT))
#lengthClass<-seq(floor(length.range[1]),ceiling(length.range[2])+1)
#lengthMids <- lengthClass[-length(lengthClass)] + bin.num/2

#Histogram of counts for each age class through time t=1:13
#freq=apply(lengthT,2,FUN=function(x){
#  hist(x,breaks=lengthClass,plot=FALSE,include.lowest=TRUE)$counts})


#timemin=as.Date("2018-01-01");timemax=as.Date("2019-01-01");dates=seq(timemin,timemax,'month')
#dates<-as.Date(dates,"%Y-%d-%m")

#Create list object to store dates, midLengths, and histogram frequencies to use in TropFishR
#lfq.dat<-list(dates=dates,
#              midLengths=lengthMids,
#              catch=as.matrix(freq))


}

