library(stringr)
library(data.table)
setDTthreads(threads = 2)
library(GenomicRanges)
library(dplyr)
source("/fh/fast/ha_g/user/tpersse/newSrc/utils.R")

args <- commandArgs(TRUE)

options(stringsAsFactors=FALSE, width=180)

inDir <- args[1] #list of files for APOLLOH LOH
geneFile <- args[2] ##/cga/meyerson/home/gavinha/references/GRCh38/gencode/GRCh38.p12.ensembl.gene.annotations.sorted.txt
sampleFile <- args[3]
method <- args[4] #{'common','severity','complete'}; default is 'severity'
headerType <- args[5] #{'Gene','chrPosn'}
filterLen <- as.numeric(args[6])
outPrefix <- args[7]
# inDir <- '/Volumes/ha_g-1/projects/NelsonLab/MedGenome_WGS_LNCaP/TitanCNA_tumorOnly_FGC_as_hetSite/results/titan/hmm/optimalClusterSolution' #list of files for APOLLOH LOH
# geneFile <- '/Volumes/ha_g-1/projects/NelsonLab/MedGenome_WGS_LNCaP/TitanCNA_tumorOnly_FGC_as_hetSite/results/titan/hmm/GRCh38.p12.ensembl.gene.annotations.sorted.include.AREnhancer.txt' ##/cga/meyerson/home/gavinha/references/GRCh38/gencode/GRCh38.p12.ensembl.gene.annotations.sorted.txt
# sampleFile <- '/Volumes/ha_g-1/projects/NelsonLab/MedGenome_WGS_LNCaP/TitanCNA_tumorOnly_FGC_as_hetSite/results/titan/hmm/sampleList_allCases.txt'
# method <- 'common' #{'common','severity','complete'}; default is 'severity'
# headerType <- 'Gene' #{'Gene','chrPosn'}
# filterLen <- 1000
# outPrefix <- 'LNCaP_allCases'
outImage <- paste(outPrefix,"_geneMats.RData",sep="")

cnExt <- ".titan.ichor.seg.txt"
overlapType <- "any"
cnCol <- "Corrected_Copy_Number"
callCol <- "Corrected_Call"
filterSnp <- 0
genomeStyle <- "UCSC"
genomeBuild <- "hg38"
seqinfo <- Seqinfo(genome=genomeBuild)
chrs <- c(1:22,"X")
# cat(chrs)
chrs <- paste0("chr",chrs, sep="")
# seqlevelsStyle(chrs) <- genomeStyle

samples <- fread(sampleFile)
# keep only column 1 of the sample file
# samples <- samples[[1]]
# print(samples)
# print(samples)

genes <- fread(geneFile)
# print(genes)
if (headerType == "Gene"){
	setnames(genes, c("cdsStart", "cdsEnd"), c("Start", "End"))
	genes <- genes[, .(Chr = Chr[1], Start = min(Start), End = max(End), Karyotype_band=Karyotype_band[1], strand=strand[1]), by=Gene]
	genes[, Length := End - Start]
	genes <- genes[Length > 1000]
	genes[, chrPosn := paste0(Chr,":", Start, "-", End)]
}else{ #if (!headerType %in% colnames(genes)){
   colnames(genes)[1:3] <- c("Chr", "Start", "End")
}
# seqlevelsStyle(genes$Chr) <- genomeStyle
genes$Chr <- paste0("chr",genes$Chr, sep="")
# cat('did we get here?')
genes <- genes[Chr %in% chrs]
genes$Chr <- factor(genes$Chr, levels = chrs)
genes <- genes[order(Chr, Start)]
numGenes <- nrow(genes)
#setkey(genes, chrPosn)

# stateKey <- array(0:9); names(stateKey) <- c('HOMD','DLOH','HET','NLOH','GAIN','ALOH','BCNA','UBCNA','ASCNA','HLAMP')
# severity <- array(c(6,5,4,5,1,2,3,3,4,5)); names(severity) <- c('HOMD','DLOH','NLOH','ALOH','HET','GAIN','BCNA','UBCNA','ASCNA','HLAMP')

## want to capture this severity: 
# Most severe: high level amp/homdel
# HLAMP: 10 < CN
# Homdel: CN = 0
# Next most: Deletion
# CN = 1
# Next most: Amp
# AMP: 5 < CN < 10
# Next: Gain
# GAIN  2 < CN < 5

# so, set 

save.image(outImage)
#########################################################################################
####################### FUNCTION: FIND OVERLAP OF SEGMENT AND GENE ######################
#########################################################################################

## input: lohHits = loh rows that overlap region of interests
# start = start coordinate of region of interest
# end = end coordinate of region of interest
## NOT USED
getOverlapLength <- function(lohHits, start, end){
	coords <- cbind(lohHits[, c("Start","Stop")], as.numeric(start), as.numeric(end))
	coordsSort <- t(apply(coords, 1, sort))
	dist <- coordsSort[, 3] - coordsSort[, 2] + 1
	return(dist)
}

#Input: states is an array of state names (e.g. (DLOH,NLOH,...,))
#Uses global variable "severity"
## NOT USED
getMostSevereState <- function(states){	
    severityValue <- 0 # initialize to 0, lowest severity
    severeState <- states[1] # take in the first state
    for (i in states){
        if (!is.na(severity[i]) && length(severity[i]) == 1){
            if (severity[i] > severityValue){ # if the severity of the current state is greater than the severity of the previous state
                severeState <- i # set the severeState to the current state
                severityValue <- severity[i] # set the severityValue to the severity of the current state
            }
		# else
        }
    }
    return(severeState)
}

# output the matrix to file
#Input: Matrix to output; output file name
writeMatrixToFile <- function(mat,outfile){
	outMat <- cbind(rownames(mat),mat)
	if (!is.null(colnames(outMat))){
		colnames(outMat)[1] <- "Sample"
	}
	write.table(outMat,file=outfile,row.names=F,col.names=T,quote=F,na="NaN",sep="\t")
}

## finds the bin-level overlap
# bins from CN data that overlap a gene 
getBinCNOverlap <- function(x, y, func = "min", type = "any", colToReturn = "Corrected_Copy_Number"){
	# print(y)
	cn <- rep(NA, nrow(x))
	x.gr <- as(x, "GRanges") # x = gene
	y.gr <- as(y, "GRanges") # y = bin-level CN
	# print(y.gr)
	# # print column ColToReturn
	# print(y$colToReturn) # this is returning NULL
	# print(y[, "Corrected_Copy_Number"]) # this is returning NULL
	# print('did we get here?')
	hits <- findOverlaps(query = y.gr, subject = x.gr, type = type)
	hits.dt <- as.data.table(hits)
	
	avgBinCN <- hits.dt[, match.fun(func)(y[queryHits, round(get(colToReturn))], na.rm=T), by = subjectHits] # match.fun(func) is a function that takes the mean of the bin-level CN
	# minBinCN <- hits.dt[, as.numeric(match.fun(min)(y[queryHits, get(colToReturn)], na.rm=T)), by = subjectHits] # get the minimum value from coltoReturn
	# minBinCN$V1[minBinCN$V1 < 0.25] <- 0
	# avgBinCN$V1[minBinCN$V1 == 0] <- 0
	cn[avgBinCN$subjectHits] <- avgBinCN$V1
	return(cn)
} 
# print(list.files(inDir))
files <- list.files(inDir, pattern = cnExt, full.names = TRUE)
# print(files)
ids <- str_extract(basename(files), "[0-9A-Za-z-_]+")
# print(ids)
names(files) <- ids
print(samples[[1]])
files <- files[samples[[1]]] # this is now returning NAs... why?
#
print(files)
numSamples <- length(files)
numGenes <- nrow(genes)

geneCallmat <- matrix(NaN,nrow=numSamples,ncol=numGenes, dimnames=list(names(files), genes[, get(headerType)]))
geneCNmat <- matrix(NaN,nrow=numSamples,ncol=numGenes, dimnames=list(names(files), genes[, get(headerType)]))
# geneCNmatSev <- matrix(NaN,nrow=numSamples,ncol=numGenes, dimnames=list(names(files), genes[, get(headerType)]))
geneLOHmat <- matrix(NaN,nrow=numSamples,ncol=numGenes, dimnames=list(names(files), genes[, get(headerType)]))
geneARmat <- matrix(NaN,nrow=numSamples,ncol=numGenes, dimnames=list(names(files), genes[, get(headerType)]))
geneHRmat <- matrix(NaN,nrow=numSamples,ncol=numGenes, dimnames=list(names(files), genes[, get(headerType)])) 
geneLogRmat <- matrix(NaN,nrow=numSamples,ncol=numGenes, dimnames=list(names(files), genes[, get(headerType)])) 
geneCFmat <- matrix(NaN,nrow=numSamples,ncol=numGenes, dimnames=list(names(files), genes[, get(headerType)]))
geneCCmat <- matrix(NaN,nrow=numSamples,ncol=numGenes, dimnames=list(names(files), genes[, get(headerType)])) 
geneLenmat <- matrix(NaN,nrow=numSamples,ncol=numGenes, dimnames=list(names(files), genes[, get(headerType)])) 
geneBinCNmat <- matrix(NaN,nrow=numSamples,ncol=numGenes, dimnames=list(names(files), genes[, get(headerType)]))

for (i in 1:numSamples){
	caseId <- names(files[i])
	cat("Analyzing sample: ",caseId," for file: ",files[i],"...(",i,"file(s) in)\n")

	#LOH
	loh <- fread(files[i])
	colnames(loh)[c(2,3,4)] <- c("Chr","Start","End")
	loh[, Length.bp := End - Start + 1]
	## filter by length threshold ##
	# print(loh)
	loh <- loh[Length.bp>=filterLen & Length.snp.>=filterSnp, ]
	loh[, LOH := as.numeric(MinorCN == 0)] # 1 if LOH, 0 if not 
	
	# CNA bin file to use for chrX CN
	bin <- fread(gsub("seg.txt", "cna.txt", files[i]))
	#binChrXind <- bin[Chr == "chrX", which = TRUE]
	#geneChrXind <- genes[Chr == "chrX", which = TRUE]
	# return and fix...
	# geneBinCN <- getBinCNOverlap(x=genes, y=bin, func = "mean", type = "any", colToReturn="logR_Copy_Number")
	geneBinCN <- getBinCNOverlap(x=genes, y=bin, func = "mean", type = "any", colToReturn=cnCol)
	# print('geneBinCN')
	# print(geneBinCN)
	#convert to gene scaffold
	geneCN <- getOverlap(x=genes, y=loh, type=overlapType, colToReturn=cnCol, method='common')
	# geneSevCN <- getOverlap(x=genes, y=loh, type=overlapType, colToReturn=cnCol, method='severity')
	cat('Done with CN calls for sample: ', caseId, '...\n')
	# print(geneCN)
	geneLOH <- getOverlap(x=genes, y=loh, type=overlapType, colToReturn="LOH")
	geneLogR <- getOverlap(x=genes, y=loh, type=overlapType, colToReturn="Median_logR")
	geneAR <- getOverlap(x=genes, y=loh, type=overlapType, colToReturn="Median_Ratio")
	if ("Median_HaplotypeRatio" %in% colnames(loh)){
    	geneHR <- getOverlap(x=genes, y=loh, type=overlapType, colToReturn="Median_HaplotypeRatio")
  	}
	cat('Done with LOH,LOGR,AR calls for sample: ', caseId, '...\n')
	geneCF <- getOverlap(x=genes, y=loh, type=overlapType, colToReturn="Cellular_Prevalence")
	geneCC <- getOverlap(x=genes, y=loh, type=overlapType, colToReturn="Clonal_Cluster")
	geneState <- getOverlap(x=genes, y=loh, type=overlapType, colToReturn="TITAN_state")
	geneLen <- getOverlap(x=genes, y=loh, type=overlapType, colToReturn="Length.bp")
	#build matrices
	geneCallmat[caseId,] <- geneState  #call matrix
	geneCNmat[caseId,] <- geneCN #call matrix
	geneLOHmat[caseId,] <- geneLOH #LOH matrix
	geneARmat[caseId,] <- geneAR   #AR matrix
	geneLogRmat[caseId,] <- geneLogR   #logR matrix
	geneCFmat[caseId,] <- geneCF  #cellular prevalence matrix
	geneCCmat[caseId,] <- geneCC  #clonal cluster matrix
	geneBinCNmat[caseId, ] <- geneBinCN # average bin-level CN 
	# geneCNmatSev[caseId,] <- geneSevCN #severity CN matrix
	if ("Median_HaplotypeRatio" %in% colnames(loh)){
		geneHRmat[caseId,] <- geneHR   #haplotypeRatio matrix
	}
	geneLenmat[caseId, ] <- geneLen
	save.image(outImage)
}
save.image(outImage)


outCall <- paste(outPrefix,"_geneCalls.txt",sep="")
writeMatrixToFile(geneCallmat,outCall)
outCN <- paste(outPrefix,"_geneCN.txt",sep="")
writeMatrixToFile(geneCNmat,outCN)
outLOH <- paste(outPrefix,"_geneLOH.txt",sep="")
writeMatrixToFile(geneLOHmat,outLOH)
outAR <- paste(outPrefix,"_geneAR.txt",sep="")
writeMatrixToFile(geneARmat,outAR)
outLogR <- paste(outPrefix,"_geneLogR.txt",sep="")
writeMatrixToFile(geneLogRmat,outLogR)
outCF <- paste(outPrefix,"_geneCF.txt",sep="")
writeMatrixToFile(geneCFmat,outCF)
outCC <- paste(outPrefix,"_geneCC.txt",sep="")
writeMatrixToFile(geneCCmat,outCC)
outHR <- paste(outPrefix,"_geneHR.txt",sep="")
writeMatrixToFile(geneHRmat,outHR)
outLen <- paste(outPrefix,"_geneLength.txt",sep="")
writeMatrixToFile(geneLenmat,outLen)
outBinCN <- paste(outPrefix,"_geneBinCN.txt",sep="")
writeMatrixToFile(geneBinCNmat,outBinCN)
# outSevCN <- paste(outPrefix,"_geneSevCN.txt",sep="")
# writeMatrixToFile(geneCNmatSev,outSevCN)


# save(geneCallmat, geneCNmat, geneLOHmat, geneARmat, geneLogRmat, geneCFmat, geneCCmat, geneHRmat, geneLenmat, geneBinCNmat, geneCNmatSev, file=outImage)
#save.image(outImage)
save(geneCallmat, geneCNmat, geneLOHmat, geneARmat, geneLogRmat, geneCFmat, geneCCmat, geneHRmat, geneLenmat, geneBinCNmat, file=outImage)
# '''
# Rscript makeMatrixFromTITAN-ICHOR_segBased_hg38.R ../../CN-SV_calling/TitanCNA/NNI/titanSegsv2/ ../../CN-SV_calling/TitanCNA/GeneCnStatus/GRCh38.p12.ensembl.gene.annotations.sorted.include.AREnhancer.txt ../../CN-SV_calling/TitanCNA/GeneCnStatus/purity_ploidy.txt severity Gene 1000 binTest 
# '''