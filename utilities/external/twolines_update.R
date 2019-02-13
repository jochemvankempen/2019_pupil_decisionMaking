
#ADDED TO DATACOLADA.ORG/62
#
#Yair Heller contacted me, Uri, via Andrew Gelman, with an example where the two-lines test had an elevated 
#False-Positive rate. I determined the source of the problem was hetersokedasticity and found that estimating
#the interrupted regression with robust standard errors fixed the problem. Then I generated different forms of 
#heteroskedasticity and found that again robust standard errors worked well.

#The code below has Yair's example, estimated with and without robust SE, and a few additional examples.
#
#This version: 2017 11 03
# Written by: Uri Simonsohn (urisohn@gmail.com)
##Please email me directly if you spot any errors or have questions or comments.
########################################################################################

#OUTLINE
#0) Load libraries 
#1) Function 1.  Simplified two-line regression function, with interruption at point xc
#2) Function 2.  Simplified Robin Hood function to find the optimal cutoff and then run  two-lines regression
#3) Original example by Yair Heller (his code was different, but achieved the same simulation)
#4) Function 3.  Simulation for monte-carlo with flexible forms of heteroskedasticity
#5) Run monte-carlo simulation of four extreme forms of heteroskedasticity relying on Functions 1-3


########################################################################################
#0) Load Libraries and set working directory
  library(sandwich)
  library(lmtest)
  library(mgcv)
  

  
    
#1) Simplified two-line regression function, with interruption at point xc
    reg2=function(x,y,xc,robust)
      { 
      #robust=1 estimate robust standard errors
      #robust=0 estimate traditional (homoskedastic-assuming) standard errors
      
      #Create new variables
        xlow1=ifelse(x<=xc,x-xc,0)     #xlow=x-xc when x<xc, 0 otherwise
        xhigh1=ifelse(x>xc,x-xc,0)     #xhigh=x when x<xmax, 0 otherwise
        high1=ifelse(x>xc,1,0)         #high dummy, allows interruption
      #Run the regressions 
        lm1=lm(y~xlow1+xhigh1+high1)     #estimate regression
        lmc1=summary(lm1)$coefficients   #Get b,t,se,p
      #Robust standard errors
        rob1=coeftest(lm1, vcov=vcovHC(lm1,"HC3"))
        
      #Report robust or homoskedastic-assuming standard errors? 
        if (robust==0) c1=lmc1   #if not robust, just use lm()
        if (robust==1) c1=rob1   #if yes robust, select results from the coeftest() output
        
      #Turn into individual results (scalars)
        b1=c1[2,1]
        t1=c1[2,3]
        p1=c1[2,4]
        
        b2=c1[3,1]
        t2=c1[3,3]
        p2=c1[3,4]
      
      #Is the u-shape significant?
        u.sig =ifelse(b1*b2<0 & p1<.05 & p2<.05,1,0)                     
      #All results
        res=list(b1=b1,b2=b2,t1=t1,t2=t2,u.sig=u.sig)
        res
    }


#2) Simplified Robin Hood function to find the optimal cutoff and then run  two-lines regression
  reg2hood=function(x,y,robust)
  #Robust =1 estimates robust standard errors (HC3), robust=0 assumes homoskedasticity
  {
  #Syntax:
  #1 Run gam()
    g=gam(y~s(x,bs="cr"))  
  #2 Get fitted values
    g.fit=predict.gam(g,se.fit=T)
    y.hat=g.fit$fit
    y.se =g.fit$se.fit
  #3 Focus on the middle 80% of the x-values (smoothed values are not reliable near the end)
    x10=quantile(x,.1)
    x90=quantile(x,.9)
    middle=(x>x10 & x<x90)       #Don't consider extreme values for cutoff
    y.ub=y.hat+1*y.se            #+1 SE is the flat region
    x.middle=x[middle]           #X values associated with flat region
    xc.max=x.middle[match(max(y.hat[middle]),y.hat[middle])]       #find value of x associated with highest predicted value
  #4 Find flat maximum  
    flat=(y.ub>max(y.hat) & middle)
    xflat=x[flat] 
  #5 If empty (e.g., if best fitting line is monotonic) use median x
    if (length(xflat)==0) xflat=median(x) 
  
  #6 Regression split based on predicted maximum
     rmax=reg2(x,y,xc=xc.max,robust=robust) 
  
  #7  Regression split based median of xflat
     rmed=reg2(x,y,xc=median(xflat),robust=robust)  #At the median of  xflat
  #8 Adjust split point based on ratio of t1,t2, move split away from less significant slope, to give more of a chance
      t1=abs(rmed$t1)             
      t2=abs(rmed$t2)             
      xc.prop=quantile(xflat,t2/(t1+t2))  #Choose the percentile value in flat region proportional to t1,t2
  #For example, if t2=t1 it stays at median flat, if t2>t1, it moves lower
  #9 Regression split based on adjusted based on t1,t2    
      rprop=reg2(x,y,xc=xc.prop,robust=robust)
      res=c(rmax$u.sig,rmed$u.sig, rprop$u.sig)
   res
}

##########################################  
#3) Original example by Yair Heller (his code was different, but achieved the same simulation)
  simtot=5000
  set.seed(123)
  res=matrix(nrow=simtot,ncol=2)
  for (simk in 1:simtot) {
    x=sort(runif(100))
    y.raw=ifelse(x<.7,x,.7)
    e1=rnorm(sum(x<.7),sd=.1) #in segment 1 sd(y)=.1
    e2=rnorm(sum(x>.7),sd=.5) #in segment 2 sd(y)=.5
    e=c(e1,e2)
    y=y.raw+e
    traditional  =reg2hood(x,y,robust=0)[3]   #reg2hood() reports 3 results, the 3rd is for Robin Hood. 1/0 is u-shape significant 
    robust       =reg2hood(x,y,robust=1)[3]   #reg2hood() reports 3 results, the 3rd is for Robin Hood. 1/0 is u-shape significant 
    res[simk,]=c(traditional,robust)          #Store both results for this simulation as a row in a matrix with 2 columns
    if (simk%%100==0) cat("...",simk)         #counter
    }
  
  #Name both columns
    names=c("traditional","robust")
    colnames(res)=names
  #The means of 1/0 variables on significant u-shape are the false-positive rates
    colMeans(res)             

    
#If curious this is what Yair's (extreme, in my view) example looks like, i use larger n to make it visually easy to spo
    x=sort(runif(1000))  #So, n=1000 for the graph to be clearer
    y.raw=ifelse(x<.7,x,.7)
    e1=rnorm(sum(x<.7),sd=.1) #in segment 1 sd(y)=.1
    e2=rnorm(sum(x>.7),sd=.5) #in segment 2 sd(y)=.5
    e=c(e1,e2)
    y=y.raw+e
    plot(x,y)
    
################<<<<--->>>>###############################################    
#let's now consider a broader set of possible heteroskedastic errors
  
#4) Function 3-  Simulation for monte-carlo with flexible forms of heteroskedasticity  
    sim=function(seed,n,noise,simtot,xc)
    {
    set.seed(seed)
    res=matrix(nrow=simtot,ncol=2)
        for (simk in 1:simtot) {
      x=sort(runif(n))
      y.raw=ifelse(x<xc,x,xc)
      e=eval(parse(text=noise)) 
      y=y.raw+e
      traditional  =reg2hood(x,y,robust=0)[3]   #reg2hood() reports 3 results, the 3rd is for Robin Hood. 1/0 is u-shape significant 
      robust       =reg2hood(x,y,robust=1)[3]   #reg2hood() reports 3 results, the 3rd is for Robin Hood. 1/0 is u-shape significant 
      res[simk,]=c(traditional,robust)          #Store both results for this simulation as a row in a matrix with 2 columns
      if (simk%%100==0) cat("...",simk)         #counter
      
    #Plot the first simulation for user to see what's happening
      if (simk==1) {
        plot(x,y,col="gray77",cex=.5,pch=16,main=paste0("e=",noise,"\nn=",n,"|   xc=",xc))
        points(x,y.raw,type='l')
      }
    }
      output=colMeans(res)    
      cat("\n False-positive rates for traditional and robust are ", output, "respectively.")
      output
    
  }
  
#5) Run the 4 scenarios, with n=100 and with n=10000
    s=matrix(nrow=8,ncol=2)
  
  #n=100
    s[1,]=sim(seed=100,simtot=500,n=100,xc=.7,noise="rnorm(n,sd=x)")            #increasing on x
    s[2,]=sim(seed=102,simtot=500,n=100,xc=.7,noise="rnorm(n,sd=.1+sqrt(x))")   #increasing, sqrt(x)
    s[3,]=sim(seed=103,simtot=500,n=100,xc=.7,noise="rnorm(n,sd=1.1-x)")        #decreasing linearly
    s[4,]=sim(seed=104,simtot=500,n=100,xc=.7,noise="rnorm(n,sd=.1+abs(x-.5))") #extremes are more error prone
    
  #Do again with n=1000
    s[5,]=sim(seed=105,simtot=500,n=1000,xc=.7,noise="rnorm(n,sd=x)")              
    s[6,]=sim(seed=107,simtot=500,n=1000,xc=.7,noise="rnorm(n,sd=.1+sqrt(x))")    
    s[7,]=sim(seed=108,simtot=500,n=1000,xc=.7,noise="rnorm(n,sd=1.1-x)")         
    s[8,]=sim(seed=109,simtot=500,n=1000,xc=.7,noise="rnorm(n,sd=.1+abs(x-.5))")  

  #Name the columns
    names=c("traditional","robust")
    colnames(s)=names
  #Show results
    s
    
  #Compare performance of both
    plot(s[,2],ylim=c(0,1),pch=16,col='blue',
         ylab='False-positive rate',
         xlab="Scenarios",
         main="Robust SE-->Robust two-lines test")
    points(s[,1],pch=16,col='red')
    abline(h=.05,col='green')
    legend(2,.8,c("Two-lines","Two-lines, with Robust SE"),col=c("red","blue"),pch=c(16,16))
    
    
    
#Now just do the plots for the data that were simulated, with n=1000
  plotit=function(seed,n,noise)
    {
      set.seed(seed)
      x=sort(runif(n))
      y.raw=ifelse(x<.5,x,.5)
      e=eval(parse(text=noise)) 
      y=y.raw+e
      plot(x,y,col="gray77",cex=.5,pch=16,main=paste0("e=",noise,"\nn=",n))
      points(x,y.raw,type='l')
    }
    
  #2x2 graph
  dev.off()
  par(mfrow=c(2,2))
    plotit(seed=105,n=1000,noise="rnorm(n,sd=x)")             #increasing on x
    plotit(seed=107,n=1000,noise="rnorm(n,sd=.1+sqrt(x))")    #increasing, sqrt(x)
    plotit(seed=108,n=1000,noise="rnorm(n,sd=1.1-x)")         #decreasing linearly
    plotit(seed=109,n=1000,noise="rnorm(n,sd=abs(x-.5))")     #more at the ends
  

  
  
    