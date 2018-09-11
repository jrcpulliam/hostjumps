# sgcp.hypothesis.R
# JRCP (Last modified: 10.08.08)


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
# 199(4): 565-568. DOI: 10.1086/596510                                     #
#                                                                          #
# If you modify this work, you may distribute the resulting derivative     #
# work only under the same or a similar license to this one and with       #
# proper attribution of the original work, as stated above. To see further #
# details of this license, including a link to the legal code, visit       #
# http://creativecommons.org/licenses/by-nc-sa/3.0/.                       #
#						                           #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


rm(list=ls())
library(utils)
library(grDevices)
library(stats)


### CONTROL PARAMETERS ###

INCLUDE.FLU<-FALSE	
# Note: By default, all analyses are done with Influenza A virus excluded (see text). To include Influenza A virus, change value of INCLUDE.FLU to TRUE.

SAVE.VALS<-TRUE
# Note: Change value to FALSE if you do not want to save the table of calculated p-values. 

REPS<-1000#00		# Number of permutations to perform for p-value calculation. 
# Note: Increase REPS to increase precision; decrease REPS to reduce run time.

#-- Set seed for random number generator --#
SEED<-484425
set.seed(SEED)
# Note: Leaving this line uncommented and setting REPS to 100000 will produce the exact results presented 
# in the paper (Pulliam and Dushoff 2009).


### PREPARATION FOR ANALYSIS ###

source("sgcp.functions.R"); source("sgcp.data.prep.R")

#-- Determine regression parameters for full model --#
fit.full<-glm(human~seg+gm+sr,data=data,family="binomial")
# Note that the regression parameters for the full model from glm() are the same as the maximum likelihood estimates of the best fit parameters for the full model from mle(), which was used in sgcp.model.R

#-- Define observed AICc statistic for comparison --#
AICc<-aic.c(fit.full)
# Note that this is exactly analogous to using -log likelihood as the test statistic since the number of parameter (k) and the number of data points (n) remain constant across permutations.

### FAMILY-LEVEL ANALYSIS ###

#-- Define function that calculates AICc for regression with a random permutation of the given trait at the family level --#
find.stat.fam<-function(i){ # 1=Seg, 2=GM, 3=SR
  switch(i,	
    { # Seg
	  temp<-rep(0,length(fam[,i+1]))
      temp[sample(c(f.011,f.111),length(f.111))]<-1 # random RNA/cyto -> seg
      temp[sample(c(f.010,f.110),length(f.110))]<-1 # random RNA/nuc -> seg
    },
    { # GM
	  temp<-rep(1,length(fam[,i+1]))
      temp[sample(c(f.001,f.011),length(f.001))]<-0 # random non-seg/cyto -> DNA 
      temp[sample(c(f.000,f.010),length(f.000))]<-0 # random non-seg/nuc  -> DNA
    },
    { # SR
      temp<-rep(1,length(fam[,i+1]))
      temp[sample(c(f.110,f.111),length(f.110))]<-0 # random seg/RNA -> nuc 
      temp[sample(c(f.010,f.011),length(f.010))]<-0 # random non-seg/RNA -> nuc 
      temp[sample(c(f.000,f.001),length(f.000))]<-0 # random non-seg/DNA -> nuc
    }
  )
  data[,i]<-temp[index]			### NOTE: only changes values in this call of this function ### 
  stat<-aic.c(glm(human~seg+gm+sr,data=data,family="binomial"))
  return(stat)
}

attach(fam)
  f.000<-which(!seg&!gm&!sr)	# Indices: Non-segmented, DNA genome, Nuclear entry required
  f.001<-which(!seg&!gm&sr) 	# Indices: Non-segmented, DNA genome, Replication completed in cytoplasm
  f.010<-which(!seg&gm&!sr)     # Indices: Non-segmentd, RNA genome, Nuclear entry required
  f.011<-which(!seg&gm&sr)      # Indices: Non-segmented, RNA genome, Replication completed in cytoplasm
  f.100<-which(seg&!gm&!sr)     # Indices: Segmented, DNA genome, Nuclear entry required
  f.101<-which(seg&!gm&sr)      # Indices: Segmented, DNA genome, Replication completed in cytoplasm
  f.110<-which(seg&gm&!sr)      # Indices: Segmented, RNA genome, Nuclear entry required
  f.111<-which(seg&gm&sr)       # Indices: Segmented, RNA genome, Replication completed in cytoplasm
detach(fam)

#-- Check that algorithm is ok for data set being used --#
if(length(f.100)!=0||length(f.100)!=0){
  warning("Segmented DNA viruses in data set. All 3 permutation algorithms must be modified to account for these.")
}

#-- Calculate p-value for Segmentation with permutation at the family level --#
Seg.stat.fam<-{invisible(sapply(seq(1:REPS),function(i){find.stat.fam(1)}))}
pvalSeg.fam<-length(which(Seg.stat.fam<=AICc))/REPS
# p-value: 0.803
# p-value with flu A added to database: 0.635

#-- Calculate p-value for genomic material with permutation at the family level --#
GM.stat.fam<-{invisible(sapply(seq(1:REPS),function(i){find.stat.fam(2)}))}
pvalGM.fam<-length(which(GM.stat.fam<=AICc))/REPS
# p-value: 0.950
# p-value with flu A added to database: 0.747

#-- Calculate p-value for site of replication with permutation at the family level --#
SR.stat.fam<-{invisible(sapply(seq(1:REPS),function(i){find.stat.fam(3)}))}
pvalSR.fam<-length(which(SR.stat.fam<=AICc))/REPS
# p-value: 0.015
# p-value with flu A added to database: 0.036


### GENUS-LEVEL ANALYSIS ###

#-- Add trait data to genus table by matching --#
gen<-data.frame(cbind(gen,fam$seg[gen.index],fam$gm[gen.index],fam$sr[gen.index]))
names(gen)[3:5]<-names(fam)[2:4]

#-- Create an index that matches each virus species to its genus (and therefore its genus-level trait) --#
new.gen.index<-match(paste(sgcp$family,sgcp$genus),paste(gen$family,gen$genus))

#-- Define function that calculates AICc for regression with a random permutation of the given trait at the genus level --#
find.stat.gen<-function(i){ # 1=Seg, 2=GM, 3=SR
  switch(i,	
    { # Seg
	  temp<-rep(0,length(gen[,i+2]))
      temp[sample(c(g.011,g.111),length(g.111))]<-1 # random RNA/cyto -> seg
      temp[sample(c(g.010,g.110),length(g.110))]<-1 # random RNA/nuc -> seg
    },
    { # GM
	  temp<-rep(1,length(gen[,i+2]))
      temp[sample(c(g.001,g.011),length(g.001))]<-0 # random non-seg/cyto -> DNA 
      temp[sample(c(g.000,g.010),length(g.000))]<-0 # random non-seg/nuc  -> DNA
    },
    { # SR
      temp<-rep(1,length(gen[,i+2]))
      temp[sample(c(g.110,g.111),length(g.110))]<-0 # random seg/RNA -> nuc 
      temp[sample(c(g.010,g.011),length(g.010))]<-0 # random non-seg/RNA -> nuc 
      temp[sample(c(g.000,g.001),length(g.000))]<-0 # random non-seg/DNA -> nuc
    }
  )
  data[,i]<-temp[new.gen.index]			### NOTE: only changes values in this call of this function ### 
  stat<-aic.c(glm(human~seg+gm+sr,data=data,family="binomial"))
  return(stat)
}

attach(gen)
  g.000<-which(!seg&!gm&!sr)	# Indices: Non-segmented, DNA genome, Nuclear entry required
  g.001<-which(!seg&!gm&sr) 	# Indices: Non-segmented, DNA genome, Replication completed in cytoplasm
  g.010<-which(!seg&gm&!sr)     # Indices: Non-segmentd, RNA genome, Nuclear entry required
  g.011<-which(!seg&gm&sr)      # Indices: Non-segmented, RNA genome, Replication completed in cytoplasm
  g.100<-which(seg&!gm&!sr)     # Indices: Segmented, DNA genome, Nuclear entry required
  g.101<-which(seg&!gm&sr)      # Indices: Segmented, DNA genome, Replication completed in cytoplasm
  g.110<-which(seg&gm&!sr)      # Indices: Segmented, RNA genome, Nuclear entry required
  g.111<-which(seg&gm&sr)       # Indices: Segmented, RNA genome, Replication completed in cytoplasm
detach(gen)

#-- Calculate p-value for Segmentation with permutation at the genus level --#
Seg.stat.gen<-{invisible(sapply(seq(1:REPS),function(i){find.stat.gen(1)}))}
pvalSeg.gen<-length(which(Seg.stat.gen<=AICc))/REPS
# p-value: 0.820
# p-value with flu A added to database: 0.672

#-- Calculate p-value for genomic material with permutation at the genus level --#
GM.stat.gen<-{invisible(sapply(seq(1:REPS),function(i){find.stat.gen(2)}))}
pvalGM.gen<-length(which(GM.stat.gen<=AICc))/REPS
# p-value: 0.958
# p-value with flu A added to database: 0.786

#-- Calculate p-value for site of replication with permutation at the genus level --#
SR.stat.gen<-{invisible(sapply(seq(1:REPS),function(i){find.stat.gen(3)}))}
pvalSR.gen<-length(which(SR.stat.gen<=AICc))/REPS
# p-value: 0.018
# p-value with flu A added to database: 0.042


### SPECIES-LEVEL ANALYSIS ###

#-- Define function that calculates AICc for regression with a random permutation of the given trait at the species level --#
find.stat.sp<-function(i){ # 1=Seg, 2=GM, 3=SR
  temp<-rep(1,length(data[,i]))
  switch(i,	
    { # Seg
	  temp<-rep(0,length(data[,i]))
      temp[sample(c(s.011,s.111),length(s.111))]<-1 # random RNA/cyto -> seg
      temp[sample(c(s.010,s.110),length(s.110))]<-1 # random RNA/nuc -> seg
    },
    { # GM
	  temp<-rep(1,length(data[,i]))
      temp[sample(c(s.001,s.011),length(s.001))]<-0 # random non-seg/cyto -> DNA 
      temp[sample(c(s.000,s.010),length(s.000))]<-0 # random non-seg/nuc  -> DNA
    },
    { # SR
      temp<-rep(1,length(data[,i]))
      temp[sample(c(s.110,s.111),length(s.110))]<-0 # random seg/RNA -> nuc 
      temp[sample(c(s.010,s.011),length(s.010))]<-0 # random non-seg/RNA -> nuc 
      temp[sample(c(s.000,s.001),length(s.000))]<-0 # random non-seg/DNA -> nuc
    }
  )
  data[,i]<-temp			### NOTE: only changes values in this call of this function ### 
  stat<-aic.c(glm(human~seg+gm+sr,data=data,family="binomial"))
  return(stat)
}

attach(data)
  s.000<-which(!seg&!gm&!sr)	# Indices: Non-segmented, DNA genome, Nuclear entry required
  s.001<-which(!seg&!gm&sr) 	# Indices: Non-segmented, DNA genome, Replication completed in cytoplasm
  s.010<-which(!seg&gm&!sr)     # Indices: Non-segmented, RNA genome, Nuclear entry required
  s.011<-which(!seg&gm&sr)      # Indices: Non-segmented, RNA genome, Replication completed in cytoplasm
  s.100<-which(seg&!gm&!sr)     # Indices: Segmented, DNA genome, Nuclear entry required
  s.101<-which(seg&!gm&sr)      # Indices: Segmented, DNA genome, Replication completed in cytoplasm
  s.110<-which(seg&gm&!sr)      # Indices: Segmented, RNA genome, Nuclear entry required
  s.111<-which(seg&gm&sr)       # Indices: Segmented, RNA genome, Replication completed in cytoplasm
detach(data)

#-- Calculate p-value for Segmentation with permutation at the species level --#
Seg.stat.sp<-{invisible(sapply(seq(1:REPS),function(i){find.stat.sp(1)}))}
pvalSeg.sp<-length(which(Seg.stat.sp<=AICc))/REPS
# p-value: 0.731
# p-value with flu A added to database: 0.500

#-- Calculate p-value for genomic material with permutation at the species level --#
GM.stat.sp<-{invisible(sapply(seq(1:REPS),function(i){find.stat.sp(2)}))}
pvalGM.sp<-length(which(GM.stat.sp<=AICc))/REPS
# p-value: 0.930
# p-value with flu A added to database: 0.648

#-- Calculate p-value for site of replication with permutation at the species level --#
SR.stat.sp<-{invisible(sapply(seq(1:REPS),function(i){find.stat.sp(3)}))}
pvalSR.sp<-length(which(SR.stat.sp<=AICc))/REPS
# p-value: 0.0002
# p-value with flu A added to database: 0.0009


### RESULTS ###

#-- Create table analogous to Table 1 --#
pval.table<-data.frame(cbind(rbind(pvalSeg.sp,pvalSeg.gen,pvalSeg.fam),rbind(pvalGM.sp,pvalGM.gen,pvalGM.fam),rbind(pvalSR.sp,pvalSR.gen,pvalSR.fam)))
names(pval.table)<-c("Segmentation","Genomic material","Site of Replication")
rownames(pval.table)<-c("Species","Genus","Family")
pval.table
# Note: Because p-values are calculated via randomized permutation tests, p-values may vary slightly but should be consistent to at least two decimal places (increasing the value of REPS will increase precision but may increase the run time considerably).

if(SAVE.VALS){ #-- Save control parameters, data, and output to an R data file. --#
  if(!INCLUDE.FLU){		# If Influenza A is not included...
    save(file="sgcp.pval.Rdata",list=c("sgcp","fam","gen","INCLUDE.FLU","REPS","use.stat","use.human","pval.table","SEED"))
  }else{				# If Influenza A is included...
	save(file="sgcp.pval.flu.Rdata",list=c("sgcp","fam","gen","INCLUDE.FLU","REPS","use.stat","use.human","pval.table","SEED"))
  }
}
