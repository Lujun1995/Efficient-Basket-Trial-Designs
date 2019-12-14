#----------------------------------------------------------------#
# This file includes R code for implementing the proposed CBHM   #
# design for basket trials.                                      #
# JAGS is needed to run MCMC part.                               #
# Two functions are included:                                    #
# 1) decidePar() for determining tuning parameter a and b;       #
# 2) CBHMdesign() for generating operating characteristics;      #
# Interim analysis and early stopping are implemented.           #
# Note: this code is reorganized in a way that can be performed  #
# on a single-core machine. The simulations were orginially run  #
# on parallel on clusters.                                       #
#----------------------------------------------------------------#



#-----Function to decide tuning parameter a and b-------------------------------#
# Inputs:                                                                       #
# cohortsize: vector format, sample size for each cohort per disease type       #
# ntype: number of disease types/subgroups                                      #
# ntrial: the number of simulated trials                                        #
# p0: the value of null hypothesis rate                                         #
# p1: the value of alternative hypothesis rate                                  #
# var.small: the prespecified small value of shrinkage parameter sigma^2 for    #
# calibrating a and b                                                           #
# var.big: the prespecified big value of shrinkage parameter sigma^2 for        #
# calibrating a and b                                                           #
#                                                                               #
# Outputs:                                                                      #
# a: the tuning parameter that characterize the relationship between test       #
# statistic T and the shrinkage parameter sigma^2                               #
# b: the tuning parameter that characterize the relationship between test       #
# statistic T and the shrinkage parameter sigma^2                               #
#-------------------------------------------------------------------------------#


decidePar<-function(cohortsize,ntype,ntrial, p0, p1, var.small, var.big){
  
  
  
  presponse <- matrix(NA, nrow= ntype, ncol = ntype)
  for (i in 1:(ntype-1)){
    presponse[i,]<- c(rep(p0, i),rep(p1, ntype-i))
  }
  presponse[ntype,] <-  rep(p1, ntype)
  ncohort <- length(cohortsize)
  medianT<-NULL
  
  #----------------------------#
  for (j in 1:ntype){
    test.stat<-matrix(NA,nrow=ntrial,ncol=ncohort)
    for (sim in 1:ntrial){
      
      set.seed(100+sim)
      n<-numeric(ntype); y<-numeric(ntype)
      for (i in 1:ncohort){
        
        y = y + rbinom(rep(1,ntype), cohortsize[i], presponse[j,]);
        n = n + cohortsize[i];
        
        p<-sum(y)/sum(n)
        x<-cbind(y, n-y)
        E <- cbind(n * p, n * (1 - p))
        T <- sum((abs(x - E))^2/E)
        if (is.nan(T)){T<-0}
        test.stat[sim,i]<-T       # store the test statistic  
        
      }
    }
    medianT<-c(medianT, apply(test.stat, 2, median)[ncohort])
    
    
  }
  heteroT <- min(medianT[1:(ntype-1)]); homoT <- medianT[ntype]
  b<-(log(var.big)-log(var.small))/(log(heteroT)-log(homoT)) 
  a<-log(var.small)-b*log(homoT)
  
  results<-list(a = a, b=b)
  results
}


#---------------------------------------------------------------------------------------#
#-----------CBHMdesign() Function Arguments-------------------------------------#
# cohortsize: vector format, sample size for each cohort per disease type       #
# ntype: number of disease types/subgroups                                      #
# p.true: the vector of true response rate in each disease type                 #
# p.null: the vector of null hypothesis rate in each type                       #
# p.target: the vector of alternative hypothesis rate in each type              #
# ntrial: the number of simulated trials                                        #
# mu.par: the mean parameter of mu's prior distribution                         #
# eff.ref: the efficacy evaluation criteria, calibrated through simulations     #
# a: the tuning parameter that characterize the relationship between test       #
# statistic T and the shrinkage parameter sigma^2                               #
# b: the tuning parameter that characterize the relationship between test       #
# statistic T and the shrinkage parameter sigma^2                               #
#-------------------------------------------------------------------------------#

CBHMdesign<-function(cohortsize,ntype=4,p.true, p.null, p.target, ntrial=5000, mu.par, eff.ref, a, b){
  
  library(R2jags)
  # effect size: log odds of response
  effectsize<-function(pa){
    
    log(pa/(1-pa))
  }
  
  # JAGS main model for MCMC 
  cbhm<-function(){
    
    for (j in 1:ntype){                #J=4, the number of arms
      y[j]~ dbin(p[j],n[j])
      p[j]<-exp(theta[j])/(1+exp(theta[j]))
      
    }
    for (j in 1:ntype){
      theta[j] ~ dnorm(mu,tau)      #hierarchical model for theta
    }
    mu ~ dnorm(mu.par, 0.01)         #prior on mu
  }
  
  
  # parameters of interest
  jags.params<-c("mu","theta","p")
  # initial value 
  jags.inits<-function(){
    list("theta"=effectsize(p.true), "mu"=mu.par)
  }
  
  
  
  p.est<-matrix(0, nrow=ntrial, ncol=ntype)                   # store the estimated response rate for each disease type
  sample.size<-matrix(0, nrow=ntrial,ncol=ntype)              # store the maximum sample size used for each disease type
  ncohort<-length(cohortsize)
  arm.count<-matrix(0, nrow=ntrial, ncol=ncohort)             # store the number of disease groups that did not stop early at each interim analysis
  test.stat<-matrix(0, nrow=ntrial, ncol=ncohort)             # store the test statistic value at each interim analysis
  
  efficacy<-NULL
  eff.prob.store<-matrix(0,nrow=ntrial,ncol=ntype)
  futstop<-0.05                                              # futulity stopping cutoff
  
  
  for (trial in 1:ntrial){
    
    set.seed(100+trial)
    n<-numeric(ntype)
    y<-numeric(ntype)
    stopping<-numeric(ntype)
    presponse<-p.true
    csize<-matrix(cohortsize, nrow=ncohort, ncol=ntype)
    
    for (i in 1:ncohort){
    
    # generate data for a new cohort of patients
      y = y + rbinom(rep(1,ntype), cohortsize[i], presponse);
      n = n + csize[i,];
    
      if (i != ncohort){ # interim analysis
        
        phat<-sum(y)/sum(n)
        obs<-cbind(y, n-y); E <- cbind(n * phat, n * (1 - phat))
        T <- sum((abs(obs - E))^2/E)
        if (is.nan(T) | (T<1)) {T<-1}
        test.stat[trial, i]<-T       # store the test statistic
        sigma2<-exp(a+b*log(T))
        if (sigma2==Inf) {sigma2<-10^4}
        tau<-1/sigma2
        arm.count[trial, i]<-length(y[which(stopping==0)])
        jags.data<-list("y","n","ntype","tau","mu.par")
        jagsfit<-jags(data=jags.data,inits=jags.inits, jags.params, n.iter=20000,model.file=cbhm)
        jagsfit.upd<-autojags(jagsfit, n.update=100, n.iter=10000)
        pres.est<-jagsfit.upd[[2]]$sims.list$p                                                    # extract the mcmc samples of p
        fut.prob<-sapply(seq(1,ntype,1), function(x) mean(pres.est[,x]>(p.null[x]+p.target[x])/2)) # calculate the futility probability
        stopping[which(fut.prob<futstop)]<-1                                                      # stop if the probability <0.05
        sample.size[trial, which(stopping==1)]<-n[which(stopping==1)]
        
        
        if (1 %in% stopping){ # update the response rate and sample size for those early stopped arms to ensure that no more patients will be treated
          presponse[which(stopping==1)]<-rep(0,length(which(stopping==1)))
          csize[(i+1),which(stopping==1)]<-rep(0,length(which(stopping==1)))
        }
        if (!(0 %in% stopping)){ # when all groups have stopped for futility
          arm.count[trial, i+1]<-0
          eff.prob<-sapply(seq(1,ntype,1), function(x) mean(pres.est[,x]>p.null[x]))
          eff.prob.store[trial,]<-eff.prob
          eff.tmp<-(eff.prob>eff.ref)*1                                                # summarize the efficacy for each disease type
          efficacy<-rbind(efficacy,eff.tmp)                                            
          p.est[trial,]<-jagsfit.upd[[2]]$mean$p                                       # store the rate estimate for each disease type
          sample.size[trial,]<-n                                                       # store the maximum sample size used for each disease subgroup
          break                                                                        # stop simulating patients if all subgroups have stopped
        }  
      } else{ # final analysis
          phat<-sum(y)/sum(n)
          obs<-cbind(y, n-y); E <- cbind(n * phat, n * (1 - phat))
          T <- sum((abs(obs - E))^2/E)
          if (is.nan(T) | (T<1)) {T<-1}
          test.stat[trial, i]<-T       # store the test statistic
          sigma2<-exp(a+b*log(T))
          if (sigma2==Inf) {sigma2<-10^4}
          tau<-1/sigma2
          arm.count[trial, i]<-length(y[which(stopping==0)])
          jags.data<-list("y","n","ntype","tau","mu.par")
          jagsfit<-jags(data=jags.data,inits=jags.inits, jags.params, n.iter=20000,model.file=cbhm)
          jagsfit.upd<-autojags(jagsfit,n.update=100, n.iter=10000)
          pres.est<-jagsfit.upd[[2]]$sims.list$p
          eff.prob<-sapply(seq(1,ntype,1), function(x) mean(pres.est[,x]>p.null[x]))
          eff.prob.store[trial,]<-eff.prob
          eff.tmp<-(eff.prob>eff.ref)*1                                                # summarize the efficacy for each disease type
          efficacy<-rbind(efficacy,eff.tmp)                                            
          p.est[trial,]<-jagsfit.upd[[2]]$mean$p                                       # store the rate estimate for each disease type      
          sample.size[trial,which(stopping==0)]<-n[which(stopping==0)]                 # store the maximum sample size used for each disease type      
          
      }  
    }
    
    
    
    
    
  }
  
  # summarize results 
  print("probability of claiming efficacy by subgroup")
  cat(formatC(c(colMeans(efficacy)), digits=2, format="f"), sep ="  ", "\n") 
  print("average number of patients used") 
  cat(formatC(c(mean(rowSums(sample.size))), digits=1, format="f"), sep ="  ", "\n") 
  print("estimated response rate by subgroup")
  cat(formatC(c(colMeans(p.est)), digits=2, format="f"), sep ="  ", "\n") 
  
  results<-cbind(test.stat, sample.size, arm.count,efficacy,eff.prob.store,p.est)
  results
  
}


# One example: scenario 3 from Table 1 in the manuscript: 

cohortsize<-c(10,10,10)
p.null<-c(0.2, 0.2, 0.2, 0.2) 
p.target<-c(0.35, 0.35, 0.35, 0.35)

p.true<-c(0.45, 0.2, 0.2, 0.45)          # scenario 3 rate
eff.ref<-c(0.891,0.8915,0.891,0.8905)    # the calibrated efficacy cutoff for CBHM design based on null scenario (i.e., scenario 1)               
ntrial <- 10
mu.par <- -1.39
a.par <- decidePar(cohortsize, ntype=4, ntrial= 5000, p0 = 0.2, p1 = 0.35, var.small = 1, var.big = 80)$a
b.par <- decidePar(cohortsize, ntype=4, ntrial= 5000, p0 = 0.2, p1 = 0.35, var.small = 1, var.big = 80)$b
#------final setting in the revision-----------#
# p.null = 0.2, p.target =0.35                 #
# log(0.2/0.8) = -1.39                         #
#----------------------------------------------#




basket3 <- CBHMdesign(cohortsize=cohortsize,ntype=4,p.true=p.true, p.null=p.null, p.target=p.target, ntrial=ntrial, mu.par=mu.par, eff.ref=eff.ref, a = a.par, b= b.par )
#save(basket3,file="cHBM30_scen3org_rev1.RData")
