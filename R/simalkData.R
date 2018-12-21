#' Individual-based model to simulate age:length data for aged and unaged fish using a von Bertalanffy growth function
#'
#' Simulates length and age data for individual fish
#' @param ages = ages simulated from simageData function
#' @param Linf.mu = Linf mean
#' @param Linf.cv = CV of Linf
#' @param k.mu = k mean
#' @param k.cv = CV of k
#' @param t0.mu = t0 mean
#' @param t0.cv = CV of t0
#' @param seed = seed
#' @return Simulated length data
#' @export
simalkData<-function(ages,Linf.mu,Linf.cv,k.mu,k.cv,t0.mu,t0.cv,seed){
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


    N=length(ages)
    
    #Create empty vector to store Linf, k, and t0 parameter estimates for each individual
    #Linf=rep(0,N)
    #k=rep(0,N)
    #t0=rep(0,N)
    
    #Create dummy params for each individual of fish age a
    make.dummy.inds<-data.frame(
      Linf=Linf.mu*rlnorm(N,0,Linf.cv),
      k=k.mu*rlnorm(N,0,k.cv),
      t0=t0.mu*rlnorm(N,0,t0.cv),
      ages=ages
    )


    #Create empty vector to store initial predicted length for the ith individual of age a using von Bertalanffy function
    lt=rep(0,N)

    #Individual-based von Bertalanffy growth model to predicted initial lengths for individuals of age a
    for(i in 1:N){
      lt[i]<-make.dummy.inds$Linf[i]*(1-exp(-make.dummy.inds$k[i]*(make.dummy.inds$ages[i]-make.dummy.inds$t0[i])))
    }

    lt.na=rep(0,N)

    for(i in 1:N){
      lt.na[i]<-make.dummy.inds$Linf[i]*(1-exp(-make.dummy.inds$k[i]*(make.dummy.inds$ages[i]-make.dummy.inds$t0[i])))
    }


    #plot length estimates using ggplot2
    df.age=data.frame(Age=ages,Length=round(lt))

    df.na=data.frame(Age=rep(NA,length(lt.na)),Length=round(lt.na))

    df=rbind(df.age,df.na)
    #p<-ggplot(data=df,aes(x=Age,y=Length))+
    # geom_point(alpha=0.05,size=4,col="black")+
    # scale_x_continuous(breaks=1:max(ages))+
    # ylim(0,max(df$Length)+50)+
    # ylab("Length (mm)")+
    # xlab("Age (yr)")+
    # theme_bw()+
    # theme(panel.grid.major = element_blank(),
     #     panel.grid.minor = element_blank(),
     #     axis.text = element_text(size=20),
     #     axis.title = element_text(size=20))

    return(df)
  }
