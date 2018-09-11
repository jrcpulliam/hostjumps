# sgcp.data.prep.R
# JRCP (Last modified: 29.04.08)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~ COPYRIGHT NOTICE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#						                           #
# © 2008, Juliet R.C. Pulliam and Jonathan Dushoff. Some Rights Reserved.  #
#									   #
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


##########
##-Data-##
##########

if(!exists("INCLUDE.FLU")){
  warning("Control parameter INCLUDE.FLU does not exist. Creating with value FALSE.")
  INCLUDE.FLU<-FALSE
}

if(INCLUDE.FLU){		
  #-- Import data on viral species and human infection --#
  sgcp<-read.delim(file="database.with.flu.txt",comment.char="#")
}else{
  #-- Import data on viral species and human infection --#
  sgcp<-read.delim(file="database.txt",comment.char="#")
}


#-- Import data on characteristics at the family level --#
fam<-read.delim(file="vir.fam.txt",header=T,comment.char="#")

if(0){ #-- Check that coding matches text descriptions of traits; all commands below should return "FALSE" --#
  any((fam$segments=="multiple")!=fam$seg)
  any((fam$nucleic.acid=="RNA")!=fam$gm)
  any((fam$site=="cytoplasm")!=fam$sr)
}
#-- Remove text descriptions of traits, leaving coded values for analysis --#
fam<-subset(fam,select=-c(segments,nucleic.acid,site))

#-- Import list of families with corresponding genera --#
gen<-read.delim(file="vir.gen.txt",header=T,comment.char="#")
  
#-- Check that all families and genera included in database are also listed in family/genus tables. --#
if(any(!sgcp$family%in%fam$family)){stop("Not all families in database have an entry in <fam>.")}
if(any(!sgcp$genus%in%gen$genus)){stop("Not all genera in database have an entry in <gen>.")}
if(any(!paste(sgcp$family,sgcp$genus)%in%paste(gen$family,gen$genus))){stop("Not all species in database match a genus represented in <gen>.")}


#-- Restrict dataset to viral species with no indication of human infection or confirmed human infection (serology considered insufficient for inclusion in analysis). --#
use.human<-c("N","Y","Y1","Y2","V")    
# Note: Adding "S" to this vector will include viral species with serological evidence of human infections in analysis; if they are included, these species will be coded by default as being able to infect humans.
sgcp<-subset(sgcp,sgcp$human%in%use.human)

use.stat<-c("tentative","approved","unassigned")
sgcp<-subset(sgcp,sgcp$status%in%use.stat)
# Note: The 2 lines above are unnecessary but show how inclusion criteria can be modified based on taxonomic status.

#-- Check that all families and genera included in family/genus tables are also listed in database. --#
if(any(!fam$family%in%sgcp$family)){fam<-subset(fam,fam$family%in%sgcp$family)}
if(any(!gen$genus%in%sgcp$genus)){gen<-subset(gen,gen$genus%in%sgcp$genus)}
if(any(!paste(gen$family,gen$genus)%in%paste(sgcp$family,sgcp$genus))){stop("Not all family/genus combinations included in database.")}


if(0){
  #-- Determine number of families represented in database and included in analysis. --#
  length(unique(sgcp$family)) 			
  #-- Determine number of genera represented in database and included in analysis. --#
  length(unique(paste(sgcp$family,sgcp$genus))) 
}

#-- Create an index that matches genera to family-level traits.
gen.index<-match(gen$family,fam$family) 
# Note: This is used in sgcp.hypothesis.R in preparation of genus table for genus-level permutation analysis.

index<-match(sgcp$family,fam$family)
# Note: This is the same index as used in sgcp.hypothesis.R for trait matching in the family-level analysis.

#-- Create data table for regression parameter estimation --#
data<-data.frame(cbind(fam$seg[index],fam$gm[index],fam$sr[index],!sgcp$human%in%c("N")))
names(data)<-c("seg","gm","sr","human")
