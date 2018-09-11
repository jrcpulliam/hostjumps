# sgcp.model.R
# JRCP (Last modified: 28.04.08)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~ COPYRIGHT NOTICE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#									   #
# © 2008, Juliet R.C. Pulliam and Jonathan Dushoff. Some Rights Reserved.  #
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
# 199(4): 565-568. DOI: 10.1086/596510	                                   #
#                                                                          #
# If you modify this work, you may distribute the resulting derivative     #
# work only under the same or a similar license to this one and with       #
# proper attribution of the original work, as stated above. To see further #
# details of this license, including a link to the legal code, visit       #
# http://creativecommons.org/licenses/by-nc-sa/3.0/.                       #
#						                           #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


rm(list=ls())
library("utils")
library("grDevices")
library("stats4")

### CONTROL PARAMETERS ###

CREATE.FIGS<-TRUE
# Note: Change value to FALSE if you do not want to create and save figures. 
FIG.TYPE<-"pdf"
# Note: Default is to save figures as PDF files. Change value to "eps" if you would prefer to save them as EPS files.

SAVE.VALS<-TRUE
# Note: Change value to FALSE if you do not want to save the logistic regression table. 

INCLUDE.FLU<-FALSE
# Note: By default, all analyses are done with Influenza A virus excluded (see text). To include Influenza A virus, change value of INCLUDE.FLU to TRUE.


### PREPARATION FOR ANALYSIS ###

source("sgcp.functions.R"); source("sgcp.data.prep.R")


### BUILD MODELS ###

#-- Determine maximum likelihood estimates of the best fit parameters for each of the 7 models. --#
# Note: Using mle() for estimation of regression parameters gives the same result as using glm() for estimation of regression parameters (see example below); the advantage of using mle() is that each of the models specifies all three coefficients (including zeros) in the mle-class object produced, facilitating concatenation of output.
fit.full<-suppressWarnings(mle(mloglike))    # Equivalently: fit.full<-glm(human~seg+gm+sr,data=data,family="binomial")
fit.seg.gm<-suppressWarnings(mle(mloglike,fixed=list(bSR=0)))
fit.seg.sr<-suppressWarnings(mle(mloglike,fixed=list(bGM=0)))
fit.gm.sr<-suppressWarnings(mle(mloglike,fixed=list(bSeg=0)))
fit.seg<-suppressWarnings(mle(mloglike,fixed=list(bGM=0,bSR=0)))
fit.gm<-suppressWarnings(mle(mloglike,fixed=list(bSeg=0,bSR=0)))
fit.sr<-suppressWarnings(mle(mloglike,fixed=list(bSeg=0,bGM=0)))

#-- Create table of parameter values and aic.c values for each model. --#
models<-NULL
for(i in 1:length(ls(pattern="fit"))){
  temp<-get(ls(pattern="fit")[i])
  models<-data.frame(rbind(models,c(coef(temp),logLik(temp),sum(coef(temp)!=0),aic.c(temp))))
}
rm(temp); names(models)[5:7]<-c("log.l","k","aic.c")

#-- Reorder table from best to worst model (by aic.c). --#
models<-models[match(sort(models$aic.c),models$aic.c),]

#-- Add delta.i and w.i (Akaike weight) values to table. --#
models<-data.frame(models,models$aic.c-min(models$aic.c),exp(-0.5*(models$aic.c-min(models$aic.c)))/sum(exp(-0.5*(models$aic.c-min(models$aic.c)))))
names(models)[8:9]<-c("delta.i","w.i")

#-- Calculate importance weights (in the following order: seg, gm, sr). --#
importance<-c(sum(models$w.i[models$bSeg!=0]),sum(models$w.i[models$bGM!=0]),sum(models$w.i[models$bSR!=0]))

#-- Define functions for model-based prediction. --#
exponent<-function(seg,gm,sr,model){x<-model[1]+seg*model[2]+gm*model[3]+sr*model[4]}
prediction<-function(exponent){1/(1+exp(-exponent))}

#-- Calculate the average model based on Akaike weights. --#
ave.model<-NULL
ave.model<-sapply(1:4,function(i){b<-sum(models[1:7,i]*models$w.i[1:7])/sum(models$w.i[1:7]); c(ave.model,b)})

#-- Determine which of the three traits have the most/least impact based on parameter estimates for the averaged model. --#
most<-which(ave.model[2:4]==max(ave.model[2:4]))
least<-which(ave.model[2:4]==min(ave.model[2:4]))

if(SAVE.VALS){ #-- Save control parameters, data, and output to an R data file. --#
  if(!INCLUDE.FLU){		# If Influenza A is not included...
    save(file="sgcp.models.Rdata",list=c("sgcp","fam","gen","INCLUDE.FLU","use.stat","use.human","models","ave.model"))
  }else{				# If Influenza A is included...
	save(file="sgcp.models.flu.Rdata",list=c("sgcp","fam","gen","INCLUDE.FLU","use.stat","use.human","models","ave.model"))
  }
}


### MODEL PREDICTIONS ###

#-- Calculate model-based probabilities of human infection for each possible combination of traits. --#
model.pred<-function(model){
  model<-as.numeric(model)

  x000<-exponent(0,0,0,model)	# exponent: Non-segmented, DNA genome, Nuclear entry required
  x001<-exponent(0,0,1,model)	# exponent: Non-segmented, DNA genome, Replication completed in cytoplasm
  x010<-exponent(0,1,0,model)	# exponent: Non-segmentd, RNA genome, Nuclear entry required
  x011<-exponent(0,1,1,model)	# exponent: Non-segmented, RNA genome, Replication completed in cytoplasm
  x100<-exponent(1,0,0,model)	# exponent: Segmented, DNA genome, Nuclear entry required
  x101<-exponent(1,0,1,model)	# exponent: Segmented, DNA genome, Replication completed in cytoplasm
  x110<-exponent(1,1,0,model)	# exponent: Segmented, RNA genome, Nuclear entry required
  x111<-exponent(1,1,1,model)	# exponent: Segmented, RNA genome, Replication completed in cytoplasm

  p000<-prediction(x000)		# probability: Non-segmented, DNA genome, Nuclear entry required
  p001<-prediction(x001)		# probability: Non-segmented, DNA genome, Replication completed in cytoplasm
  p010<-prediction(x010)		# probability: Non-segmentd, RNA genome, Nuclear entry require
  p011<-prediction(x011)		# probability: Non-segmented, RNA genome, Replication completed in cytoplasm
  p100<-prediction(x100)		# probability: Segmented, DNA genome, Nuclear entry required
  p101<-prediction(x101)		# probability: Segmented, DNA genome, Replication completed in cytoplasm
  p110<-prediction(x110)		# probability: Segmented, RNA genome, Nuclear entry required
  p111<-prediction(x111)		# probability: Segmented, RNA genome, Replication completed in cytoplasm

  return(c(p000,p001,p010,p011,p100,p101,p110,p111))
}
temp<-model.pred(ave.model)

#-- Create standard order for output. --#
ind<-match(sort(temp),temp)
#-- Put predictions from averaged model into standard order. --#
pred<-temp[ind]


# ATTACH DATA SET
attach(data)
is000<-!seg&!gm&!sr		# T/F: Non-segmented, DNA genome, Nuclear entry required
is001<-!seg&!gm&sr 		# T/F: Non-segmented, DNA genome, Replication completed in cytoplasm
is010<-!seg&gm&!sr		# T/F: Non-segmentd, RNA genome, Nuclear entry required
is011<-!seg&gm&sr		# T/F: Non-segmented, RNA genome, Replication completed in cytoplasm
is100<-seg&!gm&!sr		# T/F: Segmented, DNA genome, Nuclear entry required
is101<-seg&!gm&sr		# T/F: Segmented, DNA genome, Replication completed in cytoplasm
is110<-seg&gm&!sr		# T/F: Segmented, RNA genome, Nuclear entry required
is111<-seg&gm&sr		# T/F: Segmented, RNA genome, Replication completed in cytoplasm


if(CREATE.FIGS){ #-- Create Figure 1. --#
slicesDNA<-c(sum(is000),sum(is001),sum(is100),sum(is101))
slicesRNA<-c(sum(is010),sum(is011),sum(is110),sum(is111))
l.txt<-c("Non-segmented, Requires nuclear entry",
	"Non-segmented, Replicates in cytoplasm",
	"Segmented, Requires nuclear entry",
	"Segmented, Replicates in cytoplasm")

if(FIG.TYPE%in%c("pdf","PDF",".pdf",".PDF")){
  pdf(file="Figure1.pdf",width=10.75,height=6.6)
}else{if(FIG.TYPE%in%c("eps","EPS",".eps",".EPS")){
  postscript(file="Figure1.eps",width=10.75,height=6.6)
}else{
  warning("Unknown figure type. Figure will not be saved.")
}}
  par(mfcol=c(1,2),cex.main=2,mar=c(5,4,3,3))
  temp<-which(slicesRNA>0)
  rad<-0.9
  pie(slicesRNA[temp],labels=NA,col=c("black","dark grey","black","dark grey"),main="RNA viruses",radius=rad,density=c(-1,-1,20,20),lwd=1.5)
  temp<-which(slicesDNA>0)
  pie(slicesDNA[temp],labels=NA,col=c("black","dark grey"),main="DNA viruses",radius=rad*sqrt(sum(!data$gm)/sum(data$gm)))
  par(new=T,fig=c(0,1,0,1),mar=c(0,0,0,0),lwd=1.5)
  legend("bottom",legend=l.txt[which(slicesRNA>0|slicesDNA>0)],fill=c("black","dark grey","black","dark grey"),
    bty="n",y.intersp=1.4,density=c(-1,-1,20,20),pt.cex=2,cex=1.3,inset=-.18)
  mtext("Copyright Pulliam and Dushoff 2008",2,col="dark grey",cex=.9,adj=0.5,line=-1)
  
dev.off()
}

#-- Caluclate observed proportion able to infect humans for each set of viral traits. --#
obs<-c(sum(is000&human)/sum(is000),sum(is001&human)/sum(is001),sum(is010&human)/sum(is010),sum(is011&human)/sum(is011),
  sum(is100&human)/sum(is100),sum(is101&human)/sum(is101),sum(is110&human)/sum(is110),sum(is111&human)/sum(is111))
#-- Put observed proportions into standard order. --#
obs<-obs[ind]

#-- Calculate sample sizes for each set of viral traits. --#
n<-c(sum(is000),sum(is001),sum(is010),sum(is011),sum(is100),sum(is101),sum(is110),sum(is111))
#-- Put sample sizes into standard order. --#
n<-n[ind]

# DETACH DATA SET
detach(data)


#-- Determine labels for Figure 2. --#
labels<-NULL
if(ave.model[2]>0){
  labels<-c("Segmented")
}else{
  labels<-c("Non-segmented")
}
if(ave.model[3]>0){
  labels<-c(labels,"RNA")
}else{
  labels<-c(labels,"DNA")
}
if(ave.model[4]>0){
  labels<-c(labels,"Cytoplasm")
}else{
  labels<-c(labels,"Nucleus")
}
color<-gray(c(.4,.8,0))


#-- Calculate standard errors around the proportions able to infect humans for each set of viral traits. --#
ses<-se.prop(n,obs)
use<-which(!(is.na(ses[,1])| (ses[,1]-ses[,2])==0))


if(CREATE.FIGS){ #-- Create Figure 2. --#

if(FIG.TYPE%in%c("pdf","PDF",".pdf",".PDF")){
  pdf(file="Figure2.pdf",width=15.489583,height=6.864583)
}else{if(FIG.TYPE%in%c("eps","EPS",".eps",".EPS")){
  postscript(file="Figure2.eps",width=15.489583,height=6.864583)
}else{
  warning("Unknown figure type. Figure will not be saved.")
}}

  par(xaxt="n",omd=c(0,1,.28,1),xpd=NA,cex.axis=1.6,las=1,cex.lab=1.3,mar=c(5,3.1,2,.9))
  mid<-barplot(obs,ylim=c(0,1),width=1,yaxt="n",border="dark grey")
  if(1){  # Add standard errors to plot.
  	arrows(mid[use],ses[use,1],mid[use],ses[use,2],col="black",lwd=1,code=3,angle=90,length=.05)
  }
  axis(2,line=-2,lwd=2)
  mtext("Proportion able to infect humans",2,las=0,line=1.8,cex=1.5)
  points(mid,pred,pch=19,col="black",cex=1.2) # Predictions from the averaged model.
  text(x=mid,y=0,paste("(",n,")",sep=""),xpd=T,pos=1,cex=1.3)

  polygon(x=c(mid[5]-.6,mid[5]-.6,mid[8]+.6,mid[8]+.6),y=c(-.35,-.15,-.15,-.35),col=color[most],lty=NULL,density=-1)

  polygon(x=c(mid[3]-.6,mid[3]-.6,mid[4]+.6,mid[4]+.6),y=c(-.55,-.35,-.35,-.55),col=color[c(-most,-least)],lty=NULL,density=-1)
  polygon(x=c(mid[7]-.6,mid[7]-.6,mid[8]+.6,mid[8]+.6),y=c(-.55,-.35,-.35,-.55),col=color[c(-most,-least)],lty=NULL,density=-1)

  polygon(x=c(mid[2]-.6,mid[2]-.6,mid[2]+.6,mid[2]+.6),y=c(-.75,-.55,-.55,-.75),col=color[least],lty=NULL,density=-1)
  polygon(x=c(mid[4]-.6,mid[4]-.6,mid[4]+.6,mid[4]+.6),y=c(-.75,-.55,-.55,-.75),col=color[least],lty=NULL,density=-1)
  polygon(x=c(mid[6]-.6,mid[6]-.6,mid[6]+.6,mid[6]+.6),y=c(-.75,-.55,-.55,-.75),col=color[least],lty=NULL,density=-1)
  polygon(x=c(mid[8]-.6,mid[8]-.6,mid[8]+.6,mid[8]+.6),y=c(-.75,-.55,-.55,-.75),col=color[least],lty=NULL,density=-1)

  rectangle(-.35,-.15,mid[1]-.6,mid[4]+.6,lwd=2)
  rectangle(-.35,-.15,mid[5]-.6,mid[8]+.6,lwd=2)
  text(x=mid[7]-.6,y=-.25,labels[most],cex=2,adj=0.5,col="white")

  rectangle(-.55,-.35,mid[1]-.6,mid[2]+.6,lwd=2)
  rectangle(-.55,-.35,mid[3]-.6,mid[4]+.6,lwd=2)
  rectangle(-.55,-.35,mid[5]-.6,mid[6]+.6,lwd=2)
  rectangle(-.55,-.35,mid[7]-.6,mid[8]+.6,lwd=2)
  text(x=mid[8]-.6,y=-.45,labels[c(-most,-least)],cex=2,adj=0.5,col="white")
  text(x=mid[4]-.6,y=-.45,labels[c(-most,-least)],cex=2,adj=0.5,col="white")

  rectangle(-.75,-.55,mid[1]-.6,mid[1]+.6,lwd=2)
  rectangle(-.75,-.55,mid[2]-.6,mid[2]+.6,lwd=2)
  rectangle(-.75,-.55,mid[3]-.6,mid[3]+.6,lwd=2)
  rectangle(-.75,-.55,mid[4]-.6,mid[4]+.6,lwd=2)
  rectangle(-.75,-.55,mid[5]-.6,mid[5]+.6,lwd=2)
  rectangle(-.75,-.55,mid[6]-.6,mid[6]+.6,lwd=2)
  rectangle(-.75,-.55,mid[7]-.6,mid[7]+.6,lwd=2)
  rectangle(-.75,-.55,mid[8]-.6,mid[8]+.6,lwd=2)
  text(x=mid[8],y=-.65,labels[least],cex=2,adj=0.5)
  text(x=mid[4],y=-.65,labels[least],cex=2,adj=0.5)
  text(x=mid[6],y=-.65,labels[least],cex=2,adj=0.5)
  text(x=mid[2],y=-.65,labels[least],cex=2,adj=0.5)

  segments(mid[5]-.6,as.numeric(prediction(exponent(1,1,1,model=models[1,1:4]))),mid[8]+.6,as.numeric(prediction(exponent(1,1,1,model=models[1,1:4]))),lty=2,lwd=1) # Prediction from best model (SR only) for viruses replicating in cytoplasm.
  segments(mid[1]-.6,as.numeric(prediction(exponent(1,1,0,model=models[1,1:4]))),mid[4]+.6,as.numeric(prediction(exponent(1,1,0,model=models[1,1:4]))),lty=2,lwd=1) # Prediction from best model (SR only) for viruses replicating in nucleus.
  
  mtext("Copyright Pulliam and Dushoff 2008",3,col="dark grey",cex=.9,adj=1,line=0)
dev.off()
}

