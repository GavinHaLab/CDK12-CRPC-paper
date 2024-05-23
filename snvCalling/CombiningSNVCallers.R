# Author: Manasvita Vashisth
# The script processes data from three different variant calling tools (Varscan, Mutect2, and Strelka), cleans and merges the data, and extracts various pieces of information for further analysis. 

rm(list=ls())
setwd("/set/working/directory/here")
library(tidyverse)
library(vcfR)
library(UpSetR)
library(ggplot2)
library(data.table)
library(dplyr)

######## CHANGE PATHS TO YOUR OWN ########

mutect_path <- '/path/to/mutect2/results/'
varscan_path <- '/path/to/varscan/results/'
strelka_path <- '/path/to/strelka/results/variants/'

######## CHANGE PATH TO YOUR OWN, SEE EXAMPLE IN DIR ########
## read in sample names from file, store in vector named sample_names
sample_names <- read.table("/path/to/sample/list.txt", header=FALSE) # read in sample names from file, store in vector named sample_names. For every line, the first column is the sample name 
sample_names <- as.vector(sample_names$V1) # convert to vector

# Loop over each sample in list of sample names
for(i in 1:length(sample_names))
{
  # Sets working directory for the project
  # setwd("/fh/fast/ha_g/projects/ProstateTAN/analysis_mutational")
  
  # Initializes variables, setting each to NULL
  mutect2=NULL
  strelka=NULL
  varscan=NULL
  comb1=NULL
  comb2=NULL
  comb3=NULL
  print(sample_names[i])
  
  # Reads the data from Varscan, Mutect2, and Strelka for each sample.
  varscan=as.data.table(fread(paste(varscan_path,sample_names[i],"/",sample_names[i],".snp.hg38_multianno.txt",sep=""),header=TRUE,sep = "\t",stringsAsFactors = FALSE,na.strings=c(".", "NA")))
  mutect2=as.data.table(fread(paste(mutect_path,sample_names[i],"/pass_snvs.hg38_multianno.txt",sep=""),header=TRUE,sep = "\t",stringsAsFactors = FALSE,na.strings=c(".", "NA")))
  strelka=as.data.table(fread(paste(strelka_path,sample_names[i],"/results/variants/somatic.snvs.hg38_multianno.txt", sep=""),header=TRUE,sep = "\t",stringsAsFactors = FALSE,na.strings=c(".", "NA")))
  
  # Renames columns of each dataset to make them consistent across all three variant callers.
  colnames(mutect2)[which(colnames(mutect2)=="Otherinfo1"):ncol(mutect2)]=c("Entry1","Entry2","Entry3","Chromosome_Number","Position","ID","Reference","Alternate","Quality","Mutect_Filter","Mutect_Info","Mutect_Format","Mutect_Normal","Mutect_Tumor")
  colnames(strelka)[which(colnames(strelka)=="Otherinfo1"):ncol(strelka)]=c("Entry1","Entry2","Entry3","Chromosome_Number","Position","ID","Reference","Alternate","Quality","Strelka_Filter","Strelka_Info","Strelka_Format","Strelka_Normal","Strelka_Tumor")
  colnames(varscan)[which(colnames(varscan)=="Otherinfo1"):ncol(varscan)]=c("Entry1","Entry2","Entry3","Chromosome_Number","Position","ID","Reference","Alternate","Quality","Varscan_Filter","Varscan_Info","Varscan_Format","Varscan_Normal","Varscan_Tumor")
  
  # Filters out variants that did not pass the specific variant caller's quality control (i.e. FILTER != PASS).
  mutect2=mutect2[mutect2$Mutect_Filter=="PASS",]
  strelka=strelka[strelka$Strelka_Filter=="PASS",]
  varscan=varscan[varscan$Varscan_Filter=="PASS",]
  
  # Creates a unique identifier for each variant, concatenating Chromosome, Start Position, Reference Allele, and Alternate Allele with '_'.
  mutect2$varID=paste(mutect2$Chr,"_",mutect2$Start,"_",mutect2$Ref ,"_",mutect2$Alt, sep="")
  strelka$varID=paste(strelka$Chr,"_",strelka$Start,"_",strelka$Ref,"_",strelka$Alt, sep="")
  varscan$varID=paste(varscan$Chr,"_",varscan$Start,"_",varscan$Ref,"_",varscan$Alt, sep="")
  
  # Initialize columns in each dataframe to NA in place for the upcoming merge.
  mutect2[,c("Strelka_Filter","Strelka_Info","Strelka_Format","Strelka_Normal","Strelka_Tumor","Varscan_Filter","Varscan_Info","Varscan_Format","Varscan_Normal","Varscan_Tumor"):=NA]
  strelka[,c("Mutect_Filter","Mutect_Info","Mutect_Format","Mutect_Normal","Mutect_Tumor","Varscan_Filter","Varscan_Info","Varscan_Format","Varscan_Normal","Varscan_Tumor"):=NA]
  varscan[,c("Strelka_Filter","Strelka_Info","Strelka_Format","Strelka_Normal","Strelka_Tumor","Mutect_Filter","Mutect_Info","Mutect_Format","Mutect_Normal","Mutect_Tumor"):=NA]
  
  
  # Merge all datasets based on newly created unique variant identifier ('VarID')
  comb1=bind_rows(mutect2[,c(1:129,143)],strelka[,c(1:129,143)],varscan[,c(1:129,143)])
  comb2=distinct(comb1,varID,.keep_all=TRUE)
  comb3=merge(comb2,mutect2[,138:143],by="varID",all=TRUE)
  comb3=merge(comb3,varscan[,138:143],by="varID",all=TRUE)
  comb3=merge(comb3,strelka[,138:143],by="varID",all=TRUE)
  
  # Adds a flag/indicator for each variant caller if the variant has PASSed their respective FILTER
  comb3$mutect_flag=as.integer(comb3$varID %in% mutect2$varID)
  comb3$strelka_flag=as.integer(comb3$varID %in% strelka$varID)
  comb3$varscan_flag=as.integer(comb3$varID %in% varscan$varID)
  
  setwd("/path/to/desired/output/dir") # CHANGE THIS TO YOUR WORKING DIRECTORY
  fwrite(comb3,file=paste(sample_names[i],"_Combined_MutationCalls.txt",sep=""), sep="\t",col.names = TRUE)
  
}

