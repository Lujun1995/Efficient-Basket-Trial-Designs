#----------------------------------------------------------------#
# This file includes R code for implementing Berry et al (2013)  #
# "Bayesian hierarchical modeling of patient subpopulations:     #
# Efficient designs of phase II oncology trials".                #
#                                                                #
# JAGS is needed to run MCMC part.                               #
# Interim analysis and early stopping are implemented.           #
# Note: this code is reorganized in a way that can be performed  #
# on a single-core machine. The simulations were orginially run  #
# on parallel on clusters.                                       #
#----------------------------------------------------------------#


#-----------------Function Arguments--------------------------------------------#
# cohortsize: vector format, sample size for each cohort per disease type       #
# ntype: number of disease types/subgroups                                      #
# p.true: the vector of true response rate in each disease type                 #
# p.null: the vector of null hypothesis rate in each type                       #
# p.target: the vector of alternative hypothesis rate in each type              #
# ntrial: the number of simulated trials                                        #
# mu.par: the mean parameter of mu's prior distribution                         #
# eff.ref: the efficacy evaluation criteria, calibrated through simulations     #
#-------------------------------------------------------------------------------#






berrymethod<-function(cohortsize,ntype=4,p.true, p.null, p.target, ntrial=10000, mu.par, eff.ref){
  
  library(R2jags)
  # effect size: log odds of response
  effectsize<-function(pa,pb){
    
    log(pa/(1-pa))-log(pb/(1-pb))
  }
  
  # JAGS main model for MCMC 
  berrymodel<-function(){
    for (j in 1:ntype){                
      y[j]~ dbin(p[j],n[j])
      p[j]<-exp(theta[j]+log(p.target[j]/(1-p.target[j])))/(1+exp(theta[j]+log(p.target[j]/(1-p.target[j]))))
      
    }
    for (j in 1:ntype){
      theta[j]~ dnorm(mu,tau)      
    }
    sigma2<-1/tau                    #variance=sigma2
    mu ~ dnorm(mu.par, 0.01)         #prior on mu
    tau ~ dgamma(0.0005,0.000005)    #prior on tau
  }
  
  # parameters of interest
  jags.params<-c("mu","sigma2","theta","p")
  # initial value 
  jags.inits<-function(){
    list("theta"=effectsize(p.true,p.target),"mu"=mu.par,"tau"=1/15)
  }
  
  
  
  p.est<-matrix(0, nrow=ntrial, ncol=ntype)                   # store the estimated response rate for each disease type
  sample.size<-matrix(0, nrow=ntrial,ncol=ntype)              # store the maximum sample size used for each disease type
  ncohort<-length(cohortsize)
  var.est<-matrix(0, nrow=ntrial, ncol=ncohort)               # store the posterior mean estimate of shrinkage parameter at each interim analysis
  arm.count<-matrix(0, nrow=ntrial, ncol=ncohort)             # store the number of disease groups that did not stop early at each interim analysis
  
  
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
        arm.count[trial, i]<-length(y[which(stopping==0)])
        jags.data<-list("y","n","ntype","p.target","mu.par")
        jagsfit<-jags(data=jags.data,inits=jags.inits, jags.params, n.iter=20000,model.file=berrymodel)
        jagsfit.upd<-autojags(jagsfit, n.update=100, n.iter=10000)
        pres.est<-jagsfit.upd[[2]]$sims.list$p                                                    # extract the mcmc samples of p
        fut.prob<-sapply(seq(1,ntype,1), function(x) mean(pres.est[,x]>(p.null[x]+p.target[x])/2)) # calculate the futility probability
        stopping[which(fut.prob<futstop)]<-1                                                      # stop if the probability <0.05
        sample.size[trial, which(stopping==1)]<-n[which(stopping==1)]
        var.est[trial, i]<-jagsfit.upd[[2]]$mean$sigma2                                           # store the sigma2 mean estimate
        
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
        
          arm.count[trial, i]<-length(y[which(stopping==0)])
          jags.data<-list("y","n","ntype","p.target","mu.par")
          jagsfit<-jags(data=jags.data,inits=jags.inits, jags.params, n.iter=20000,model.file=berrymodel)
          jagsfit.upd<-autojags(jagsfit,n.update=100, n.iter=10000)
          pres.est<-jagsfit.upd[[2]]$sims.list$p
          eff.prob<-sapply(seq(1,ntype,1), function(x) mean(pres.est[,x]>p.null[x]))
          eff.prob.store[trial,]<-eff.prob
          eff.tmp<-(eff.prob>eff.ref)*1                                                # summarize the efficacy for each disease type
          efficacy<-rbind(efficacy,eff.tmp)                                            
          p.est[trial,]<-jagsfit.upd[[2]]$mean$p                                       # store the rate estimate for each disease type      
          sample.size[trial,which(stopping==0)]<-n[which(stopping==0)]                 # store the maximum sample size used for each disease type      
          var.est[i]<-jagsfit.upd[[2]]$mean$sigma2                                     # store the sigma2 mean estimate
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
  
  results<-cbind(var.est, sample.size, arm.count,efficacy,eff.prob.store,p.est)
  results
  
}


# One example: scenario 3 from Table 1 in the manuscript: 

cohortsize<-c(10,10,10)
p.null<-c(0.2, 0.2, 0.2, 0.2) 
p.target<-c(0.35, 0.35, 0.35, 0.35)

p.true<-c(0.45, 0.2, 0.2, 0.45)          # scenario 3 rate
eff.ref<-c(0.853,0.853,0.850,0.850)      # the calibrated efficacy cutoff based on null scenario (i.e., scenario 1)               
ntrial <- 5000
mu.par <- -0.767
#------final setting in the revision-----------#
# p.null = 0.2, p.target =0.35                 #
# log(0.2/0.8)-log(0.35/0.65) = -0.767         #
#----------------------------------------------#




berry_scen3 <- berrymethod(cohortsize=cohortsize,ntype=4,p.true=p.true, p.null=p.null, p.target=p.target, ntrial=ntrial, mu.par=-0.767, eff.ref=eff.ref )
save(berry_scen3,file="berry30_scen3org_rev1.RData")
