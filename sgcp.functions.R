# sgcp.functions.R
# JRCP 15.04.08

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~ COPYRIGHT NOTICE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#						                           #
# ¬Â© 2008, Juliet R.C. Pulliam and Jonathan Dushoff. Some Rights Reserved.  #
#						                           #
# This work is Copyright Juliet R.C. Pulliam and Jonathan Dushoff 2008     #
# and is made available under a Creative Commons Attribution Noncommercial #
# Share Alike United States 3.0 License. You are free to copy, distribute, #
# and modify this work, provided that you do not use the work for          #
# commercial purposes, you credit the authors, and you include a copy of   #
# this notice.  If this work is used for academic purposes, the associated #
# paper should be cited:                                                   #
#						                           #
# Pulliam, JRC and Dushoff, J (2009) Ability to replicate in the cytoplasm #
# predicts zoonotic transmission of livestock viruses. JID                 #
# 199(4): 565-568. DOI: 10.1086/596510                                     #
#                                                                          #
# If you modify this work, you may distribute the resulting derivative     #
# work only under the same or a similar license to this one and with       #
# proper attribution of the original work, as stated above. To see further #
# details of this license, including a link to the legal code, visit       #
# http://creativecommons.org/licenses/by-nc-sa/3.0/.                       #
#						                           #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


#################
##- Functions -##
#################

#-----------------------------#
#-- Model-related functions --#
#-----------------------------#

## PROB() -- Calculate the expected probability of infecting humans for a given set of parameter values
prob<-function(betaI,betaSeg,betaGM,betaSR){
  x<-betaI+data$seg*betaSeg+data$gm*betaGM+data$sr*betaSR
  p<-1/(1+exp(-x))
  (log(p)^data$human)*(log(1-p)^abs(data$human-1))
}

## MLOGLIKE() -- Calculate minus log likelihood for logistic regression parameter estimation using MLE()
mloglike<-function(bI=0,bSeg=0,bGM=0,bSR=0){
  ll<--sum(prob(bI,bSeg,bGM,bSR))
  if(is.finite(ll)){ll}else{2*10^9}
}

## AIC.C() -- Calculate AIC with small sample size correction (AICc)
aic.c<-function(fit,k=sum(coef(fit)!=0),n=dim(data)[1]){
  -2*logLik(fit)[1]+2*k*(n/(n-k-1))
}

#-------------------------#
#-- Graphical functions --#
#-------------------------#

## RECTANGLE() -- Draw a rectangle
rectangle<-function(bottom,top,left,right,color="black",...){
  segments(left,bottom,left,top,col=color,...)
  segments(left,top,right,top,col=color,...)
  segments(right,top,right,bottom,col=color,...)
  segments(right,bottom,left,bottom,col=color,...)
}

#-----------------------#
#-- Utility functions --#
#-----------------------#

## SE.PROP() -- Calculate standard errors around an obsevered proportion
se.prop<-function(samp.size,p.hat){
	se<-sqrt(p.hat*(1-p.hat)/samp.size)
  return(data.frame(minus.se=p.hat-se,plus.se=p.hat+se))
}

## SHUFFLE() -- Randomly reorder a vector 
shuffle<-function(v){
  return(v[order(runif(length(v)))])
}


