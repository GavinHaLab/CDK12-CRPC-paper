#The following is a modified version of /fh/fast/ha_g/user/gha/software/git/scripts/titan/TDPfromTITAN.R
#Changes made can be found documented next to "########" --Anna Hoge 11/12/19

## requires R-3.3 #
library(VariantAnnotation)
library(data.table)
library(dplyr)
library(ggplot2)
library(reshape2)
library(stringr)
library(ggrepel)
library(GenomicRanges)
 
source("/fh/fast/ha_g/user/gha/software/git/scripts/util/Rscripts/ggMultiPlot.R") # should we include this one...?

##########commented out
#seqinfo <- readRDS("~/software/code/git/scripts/references/Seqinfo_hg19.rds")
#seqinfo <- readRDS("/fh/fast/ha_g/user/gha/software/git/scripts/references/Seqinfo_hg19.rds")

#source("/fh/fast/ha_g/user/gha/software/git/scripts/snowman/utils.R")
args <- commandArgs(TRUE)
options(stringsAsFactors=F, width=175, scipen=999, bitmapType="cairo")

cnDir <- args[1]
maxDupLength <- as.numeric(args[2]) # 1e7
scoresFile <- args[3]
geneFile <- args[4]
keyFile <- args[5]
outPrefix <- args[6]

#########added
genomeBuild <- args[7]
seqinfo <- Seqinfo(genome="GRCh38.p13")[c(1:22,"X","Y")]
# add chr prefix to seqnames
seqinfo <- renameSeqlevels(seqinfo, paste0("chr", seqlevels(seqinfo)))
bsg <- paste0("BSgenome.Hsapiens.UCSC.", genomeBuild)
if (!require(bsg, character.only=TRUE, quietly=TRUE, warn.conflicts=FALSE)) {
    seqinfo <- Seqinfo(genome=genomeBuild)
} else {
    seqinfo <- seqinfo(get(bsg))
}


outImage <- paste0(outPrefix, ".RData")
scoreType <- "gene"
scoreValues <- c("CDK12", "TP53")
excludeValues <- c("Intron")
custom_palette = c("Multiple" = '#8c47d1', "Nonsense_Mutation"="#d668b9", "Frame_Shift_Ins"='#d3ab4e',
                   "Frame_Shift_Del"='#d28e4e', "In_Frame_Ins"='#b943cc',
                   "In_Frame_Del"='#b943cc', "Nonstop_Mutation"='#d668b9',
                   "Translation_Start_Site"='#9159CF', "Splice_Site"='#f4fc14',
                   "Missense_Mutation"='#94d153', "5\'Flank"='#00A8A8',
                   "3\'Flank"='#79F2F2', "5\'UTR"='#006666',
                   "3\'UTR"='#002AA8', "RNA"='#5977CF', "Intron"='#F37812',
                   "Start_Codon_SNP" = '#9159CF', "Stop_Codon_Del" = '#9159CF',
                   "IGR"='#F2B079', "Silent"='#888811',
                   "Targeted_Region"='#FDF31C', "Rearrangement"="grey", "2 Copy Loss"='black', "Homozygous Deletion" ='black', "Amplification" = '#e60000', "Superenhancer" = '#9b0000',
                   "Copy neutral LOH" = '#80b3ff', "Gain" = "#ed867d", "1 Copy Loss - LOH" = "#80b3ff", "Deletion - LOH" = "#80b3ff", "Germline frameshift" = "#39406d", "Germline stopgain" = "#3e4db2",
                   "Chrplxy - DEL" = '#e5d487', "Chrplxy - Break" = '#eda710', "None" = "white")

#minTF <- 0.05
#maxDupLength <- 10e7
#minNeighborDupDist <- 2e6
TPD.positive.threshold <- 0.5
plotNNIthreshold <- 0.3
plotNumGainsThreshold <- 30
minLength <- 1e5
flankCNdiff <- 0
gainFrac <- 0.75
maxMad <- 0.25
minPurity <- 0.20
ncol <- 5
nrow <- 5
plotHeight <- 30
plotWidth <- 30

chrs <- c(1:22, "X")

segFiles <- list.files(cnDir, pattern = "titan.ichor.seg.txt", full.names=T)
numSamples <- length(segFiles)
## Change TP, 2023.09.18: changed regex to reflect addition of read length in sample names (-50, -100) ##
segIds <- str_match(basename(segFiles), "([0-9]+-[A-Za-z0-9_-]+)_cluster")[, 2]
key <- data.table(read.delim(keyFile, header=T, as.is=T)[,1:2])
setkey(key, ID)
ids <- sapply(segIds, function(x){ key[x, SubjectID][1] })
names(segFiles) <- ids


getGenomeWideCoords <- function(segs, chrLen){
  gwSegs <- copy(segs)
  gwSegs[, Start := as.numeric(Start)]
  gwSegs[, End := as.numeric(End)]

  chrs <- as.numeric(unique(segs$Chromosome))

  chrLen <- as.numeric(chrLen)

  gwSegs[Chromosome == 1, Start.GW := Start]
  gwSegs[Chromosome == 1, End.GW := End]

  for (i in 2:length(chrs)){
    gwSegs[Chromosome == chrs[i], Start.GW := Start + sum(chrLen[1:(chrs[i]-1)])]
    gwSegs[Chromosome == chrs[i], End.GW := End + sum(chrLen[1:(chrs[i]-1)])]
  }
  return(gwSegs)
}
save.image(outImage)

segsAll <- NULL
for (i in 1:numSamples){
  id <- names(segFiles)[i]
  message("Loading data for ", id, appendLF = F)
  ## load copy number data from TITAN ##
  cna <- fread(gsub("titan.ichor.seg.txt", "titan.txt", segFiles[i]))
  # print(cna)
  params <- read.delim(gsub("titan.ichor.seg.txt", "params.txt", segFiles[i]), header=F, as.is=T, nrow=5)
  purity <- 1 - as.numeric(params[1,2])
  diffLR <- diff(cna$LogRatio); diffLR <- diffLR[diffLR > 0]
  # print(diffLR)
  # print(mad(diffLR))
  diffMAD <- mad(diffLR)
  if (diffMAD > maxMad){
    message("... high mad")
  #  next
  }
  segs <- fread(segFiles[i])
  colnames(segs)[c(1:4)] <- c("ID", "Chromosome", "Start", "End")

  #######added this so that "chr1" naming in seg file is converted to "1" naming,
  #######for compatibility with getGenomeWideCoords function
  #if "chr" is in the first chromosome name in segs
  if(grepl("chr", segs$Chromosome[1])){
    #remove the first 3 letters ("chr") of every chromosome name
    segs$Chromosome <- gsub("^.{0,3}", "", segs$Chromosome)
  }
# note: := is a data.table operator that modifies by reference, i.e. it modifies the original data.table
  segs <- segs[Chromosome %in% chrs]
  segs <- getGenomeWideCoords(segs, seqlengths(seqinfo)) # get genome wide coordinates
  cat(colnames(segs), "\n")
  segs[, ID2 := id]
  cat(colnames(segs), "\n")
  segs[, Length.snp. := Length.snp.] 
  segs[, length := End - Start + 1] # length of segment in question
  segs[, Segment_Mean := Median_logR]
  segs[, Purity := purity]
  segs[, snpMAD := diffMAD]
  segs <- segs[, .(ID, ID2, Chromosome, Start, End, Start.GW, End.GW, Length.snp., length,
                    Segment_Mean, Copy_Number, MinorCN, Clonal_Cluster, Purity, snpMAD)] # this is where some the output columns for _TDPsegs.txt are defined... note that startGW and endGW are later removed

  segsAll <- rbind(segsAll, segs)
  message("")
}
segsAll[, ID := ID2]
segs <- copy(segsAll)


## filter segments ##
segs <- segs[Length.snp. >= 10]# & length >= minLength & length <= 20e6] # keep only those segments whose snp length is >= 10

## assign distance to beginning or end of chromosomes ##
segs[, minDistToChrEnd.GW := {chrStart = Start.GW[1]
            chrEnd = End.GW[.N]
            pmin(Start.GW - chrStart, chrEnd - End.GW)
            }, by=c("ID")]
segs[, minDistToChrEnd := {chrStart = Start[1]
            chrEnd = End[.N]
            pmin(Start - chrStart, chrEnd - End)
            }, by=c("ID", "Chromosome")]
segs[, MAD := mad(abs(diff(Segment_Mean)), na.rm=T), by = c("Chromosome", "ID")]
##################################################
########## DEFINE TANDEM DUPLICATIONS ############
##################################################
#####commented this out
#segs[, flankInd:=NULL]
segs[, flankInd := {
                    lengthInd = length < maxDupLength
                    #leftDiff = Segment_Mean - c(NA, Segment_Mean[-.N])
                    leftDiff = Copy_Number - c(NA, Copy_Number[-.N])
                    #leftLen = c(NA, length[-.N]) > length
                    leftInd = leftDiff > 0 #& leftLen
                    #rightDiff = Segment_Mean - c(Segment_Mean[-1], NA)
                    rightDiff = Copy_Number - c(Copy_Number[-1], NA)
                    #rightLen = c(length[-1], NA) > length
                    rightInd = rightDiff > 0 #& rightLen
                    #flankDiff = abs(c(NA, Segment_Mean[-.N]) - c(Segment_Mean[-1], NA))
                    flankDiff = abs(c(NA, Copy_Number[-.N]) - c(Copy_Number[-1], NA)) # difference between copy number of flanking segments
                    flankInd = flankDiff <= 1#MAD/2 # keep only those segments whose copy number difference between flanking segments is less than half the MAD
                    leftInd & rightInd & flankInd & lengthInd # keep only those segments that pass all of the above filters
                    }, by=c("Chromosome", "ID")]
# hist(log10(segs$length), 1000, xlab="Length (log10)")
segs[, segID := .I]
#setkey(segs, segID)
segs.gain <- segs[flankInd == TRUE] # gains in new table are those with flanking segments

## analyze dispersion ##
# assumes sorted by coordinate within chromosomes
# get minimum distance to left and right event (passing size threshold before)
segs.gain[, interDupDist.GW := {
            next.start = c(Start.GW[-1], NA) - End.GW
            #min(next.start, minDistToChrEnd.GW)
            }, by=c("ID")]
segs.gain[, interDupDist := {
            next.start = c(Start[-1], NA) - End
            }, by=c("ID", "Chromosome")]
# assign first and last event of each chr with dist to end of chr
segs.gain[is.na(interDupDist.GW), interDupDist.GW := minDistToChrEnd.GW]
segs.gain[is.na(interDupDist), interDupDist := minDistToChrEnd]
segs[segs.gain$segID, interDupDist.GW := segs.gain$interDupDist.GW]
segs[segs.gain$segID, interDupDist := segs.gain$interDupDist]
#segs[flankInd==T, interDupDist := segs.gain$interDupDist]
#ulp.gain <- ulp.gain[interDupDist > minNeighborDupDist]

## nearest neighbor index (NNI) ##
# need to compute from original seg dt because it has full chromosome coordinates
segs[, NNI := {
          nni <- NA
          N <- sum(!is.na(interDupDist)) # sum inter duplication distances
          if (N > 1){
            nni <- mean(interDupDist, na.rm=T) / (0.5*((End[.N]-Start[1])/sum(!is.na(interDupDist))))
          #}else if (N == 1){
          #  nni <- -1
          }
          as.numeric(nni)
          }, by=c("Chromosome","ID")]
segs[is.na(NNI), NNI := 0]
segs[NNI == -1, NNI := NA]
#segs[, NNI.GW := mean(interDupDist.GW, na.rm=T) / (0.5*((End.GW[.N]-Start.GW[1])/(sum(!is.na(interDupDist.GW))))), by=c("ID")]
#segs[is.na(NNI.GW), NNI.GW := 0]
# reassign NNI to segs.gain
segs.gain[, NNI := segs[segs.gain$segID, NNI]]
segs.gain[, Start.GW:=NULL]
segs.gain[, End.GW:=NULL]
segs.gain[, minDistToChrEnd.GW:=NULL]
segs.gain[, interDupDist.GW:=NULL]

#setkey(ulp.gain.all, ID, chrom)
## summary by sample ##
#segs.gain.summary <- segs.gain[, .(Num.Gains=.N, Median.Length=median(length), Median.InterDupDist=median(interDupDist,na.rm=T)), by=ID]
segs[flankInd==FALSE | is.na(flankInd), length := NA] # set length to NA for non-TDPs
segs.gain.summary <- segs[, .(Num.Gains=sum(flankInd, na.rm=T), Median.Length=median(length,na.rm=T), Median.InterDupDist=median(interDupDist,na.rm=T)), by=ID]
#segs.gain.summary[, patient := str_match(ID, "(.+).ctDNA")[, 2]]
## summary by sample ##

## summarize by chromosome and sample
# group by ID and chrom
segs.gain.per.chr <- segs[, .(Num.Gains.per.chr=sum(flankInd, na.rm=T), Purity=mean(Purity,na.rm=T), snpMAD=mean(snpMAD,na.rm=T), Median.Length.per.chr=median(length,na.rm=T), Median.InterDupDist.per.chr=median(interDupDist,na.rm=T), NNI=mean(NNI,na.rm=T)), by=c("Chromosome","ID")]
# per chr column names
perChrCols <- c("Num.Gains.per.chr", "Median.Length.per.chr", "Median.InterDupDist.per.chr", "NNI", "Purity", "snpMAD")
# group by ID to get median across all chromosomes
#ulp.gain.all.per.chr <- ulp.gain.all.per.chr[, lapply(.SD, mean, na.rm=T), by=ID, .SDcols=perChrCols]
segs.gain.per.chr.summary <- segs.gain.per.chr[, .(Num.Gains.per.chr=sum(Num.Gains.per.chr) / length(chrs), Purity=mean(Purity), snpMAD=mean(snpMAD),
    Median.Length.per.chr=median(Median.Length.per.chr,na.rm=T),
    Median.InterDupDist.per.chr=median(Median.InterDupDist.per.chr,na.rm=T),
    NNI=sum(NNI, na.rm=T)/length(chrs)), by=ID]
segs.gain.summary[, (perChrCols) := segs.gain.per.chr.summary[, perChrCols, with=F]]

# add back to original summary which is grouped by ID
#segs.gain.summary[order(NNI, decreasing=T)][Median.Length <= 3.5e6 & Num.Gains.per.chr >= 2 & Median.InterDupDist < 20e6,]

outFile <- paste0(outPrefix, "_summary.txt")
write.table(segs.gain.summary, file=outFile, col.names=T, row.names=F, quote=F, sep="\t")

outFile <- paste0(outPrefix, "_TDPsegs.txt")
mat <- copy(segs[Purity > minPurity & snpMAD <= maxMad]) # filter by purity and snpMAD (median absolute deviation)
mat[flankInd==FALSE | is.na(flankInd), Copy_Number := NA] # set copy number to NA for non-TDPs
mat[, Start.GW := NULL]; mat[, End.GW := NULL]; mat[, interDupDist.GW := NULL] 
write.table(mat, file=outFile, col.names=T, row.names=F, quote=F, sep="\t")

save.image(file = outImage)

###############################################################
###### Plot NNI vs length ###########
###############################################################


#mat <- segs.gain.summary[segs.gain.summary[, .I[tumor.fraction == max(tumor.fraction)], by = patient]$V1]
mat <- copy(segs.gain.summary)
mat <- mat[Purity > minPurity & snpMAD <= maxMad]
mat[NNI > plotNNIthreshold, ID.label := ID]
gp <- ggplot(mat, aes(x=Median.Length, y=NNI, size=Num.Gains.per.chr)) +
      geom_point(alpha = 0.5) +  #scale_y_continuous(trans="log10") +
      #geom_text(aes(label=ID), size=2, hjust=0) +
      geom_text_repel(aes(label=ID.label), size=3) +
      #ylab("Median Inter-duplication Distance (Mb)") +
      #ylab("Mean TD Segments Per Chromosome") +
      #ylab("Number of TD Segments") +
      #ylab("Nearest Neighbor Index (NNI)") + #ylim(0,1) +
      xlab("Median TD Segment Length (Mb)") +
      ylab("TDP Dispersion Score (NNI)") +
      scale_x_continuous(trans="log10", limits=c(1e4,10e6),
         breaks=c(1e4,1e5,1e6,2e6,3e6,4e6,5e6,10e6), labels=c(0.01,0.1,1,2,3,4,5,10)) +
      #scale_y_continuous(trans="log10") +
      #scale_y_continuous(trans="log10", limits=c(5,200),
      #   breaks=c(1,5,10,25,50,75,100,200), labels=c(1,5,10,25,50,75,100,200)) +
      theme_bw() +
      theme(axis.text.y=element_text(size=20), axis.text.x=element_text(size=16),
          axis.title=element_text(size=20),
          plot.title=element_text(size=24, face="bold"))
outPlot <- paste0(outPrefix, "_nni-vs-length.pdf")
ggsave(gp, file=outPlot, width=8, height=7)

mat <- copy(segs.gain.summary)
mat <- mat[Purity > minPurity & snpMAD <= maxMad]
mat[Num.Gains > plotNumGainsThreshold, ID.label := ID]
#mat[NNI > 1.5, ID.label := ID]
gp <- ggplot(mat, aes(x=Median.Length, y=Num.Gains, size=Num.Gains.per.chr)) +
      geom_point(alpha = 0.5) +  #scale_y_continuous(trans="log10") +
      #geom_text(aes(label=ID), size=2, hjust=0) +
      geom_text_repel(aes(label=ID.label), size=3) +
      #ylab("Median Inter-duplication Distance (Mb)") +
      #ylab("Mean TD Segments Per Chromosome") +
      ylab("Number of TD Segments") +
      xlab("Median TD Segment Length (Mb)") +
      scale_x_continuous(trans="log10", limits=c(1e4,10e6),
         breaks=c(1e4,1e5,1e6,2e6,3e6,4e6,5e6,10e6), labels=c(0.01,0.1,1,2,3,4,5,10)) +
      scale_y_continuous(trans="log10") +
      #scale_y_continuous(trans="log10", limits=c(5,200),
      #   breaks=c(1,5,10,25,50,75,100,200), labels=c(1,5,10,25,50,75,100,200)) +
      theme_bw() +
      theme(axis.text.y=element_text(size=20), axis.text.x=element_text(size=16),
          axis.title=element_text(size=20),
          plot.title=element_text(size=24, face="bold"))
outPlot <- paste0(outPrefix, "_numTD-vs-length.pdf")
ggsave(gp, file=outPlot, width=9, height=7)

mat <- copy(segs.gain.summary)
mat[NNI < 0.20, ID := NA] 
mat[, TDP := "N"]
mat[NNI >= 0.20, TDP := "Y"]
gp <- ggplot(mat, aes_string(x="Median.Length", y="NNI", color="TDP", fill="TDP")) +
      #geom_point(aes(alpha = tumor.fraction)) +  #scale_y_continuous(trans="log10") +
      geom_point(size=4, shape=21, color="black") +
      #geom_text(aes(label=ID), size=2, hjust=0) +
      #scale_color_manual(breaks=c("N", "Y"), values=c("grey", "red")) +
      scale_fill_manual(values=c("red", "grey")) +
      geom_text_repel(aes(label=ID), size=2, color = "black") +
      #ylab("Median Inter-duplication Distance (Mb)") + 
      #ylab("Mean TD Segments Per Chromosome") +
      #ylab("Number of TD Segments") +
      ylab("TD Dispersion Score (NNI)") + ylim(0, 1.8) +
      xlab("Median TD Segment Length (Mb)") +
      scale_x_continuous(trans="log10", limits=c(5e5,5e6),
         breaks=c(5e5,1e6,2e6,3e6,4e6,5e6), labels=c(0.5, 1,2,3,4,5)) +
      #scale_y_continuous(trans="log10") +
      #scale_y_continuous(trans="log10", limits=c(5,200),
      #   breaks=c(1,5,10,25,50,75,100,200), labels=c(1,5,10,25,50,75,100,200)) +
      theme_bw() + theme(legend.position = "none") +
      theme(axis.text.y=element_text(size=20), axis.text.x=element_text(size=16),
          axis.title=element_text(size=20), 
          plot.title=element_text(size=24, face="bold")) 
outPlot <- paste0(outPrefix, "_NNI-vs-length_labels.pdf")
ggsave(gp, file=outPlot, width=5, height=4)

###############################################################
###### Plot provided score (e.g. expr) vs NNI ###########
###############################################################
if (scoresFile != "0"){
  scoresOrig <- fread(scoresFile)
  scoreValues <- read.delim(geneFile, header=T, as.is=T)[,1]
  for (i in 1:length(scoreValues)){
    scoreVal <- scoreValues[i]
    scores <- scoresOrig[get(scoreType) %in% scoreVal & !variant_class %in% excludeValues, ]
    if (nrow(scores) == 0){ next }
    message("Plotting for ", scoreVal)
    key <- data.table(read.delim(keyFile, header=T, as.is=T)[,1:2])
    setkey(key, otherID)
    ids <- sapply(scores$sample, function(x){ key[x, sample][1] })
    scores[, ID := sample]

    mat <- copy(segs.gain.summary)
    mat[is.na(Num.Gains), Num.Gains := 0]; mat[is.na(NNI), NNI := 0];
    mat[is.na(Num.Gains.per.chr), Num.Gains.per.chr := 0]
    #mat[is.na(Median.Length), Median.Length := 0]
    setkey(mat, ID)
    mat[[scoreVal]] <- "None"
    mat[scores$sample, eval(scoreVal) := scores$variant_class]
    varOrder <- unique(mat[[scoreVal]]); varInd <- which(varOrder=="None" | varOrder=="Multiple")
    varOrder <- c("None", as.vector(varOrder[-c(varInd)]), "Multiple")
    mat[[scoreVal]] <- factor(mat[[scoreVal]], levels = varOrder)
    mat <- mat[order(get(scoreVal))]
    #mat[Num.Gains > plotNumGainsThreshold, ID.label := ID]
    mat[, ID.label := NULL]
    mat[NNI > plotNNIthreshold, ID.label := ID]
    mat.unFilt <- mat[ID %in% key$sample];
    mat.plot <- mat[Purity > minPurity & snpMAD <= maxMad]
    gp <- ggplot(mat.plot, aes_string(x="Median.Length", y="NNI", fill=scoreVal, color=scoreVal)) +
#gp <- ggplot(mat, aes(x=Median.Length, y=NNI, size=Num.Gains.per.chr)) +
      geom_point(size=4, alpha=0.75, shape=21) +
      #geom_text_repel(aes(label=ID.label), color="black", size=3) +
      scale_fill_manual(values=custom_palette) +
      scale_colour_manual(values=rep("black",length(mat[,unique(get(scoreVal))]))) +
       xlab("Median TD Segment Length (Mb)") +
      ylab("TDP Dispersion Score (NNI)") +
      scale_x_continuous(trans="log10", limits=c(1e4,10e6),
         breaks=c(1e4,1e5,1e6,2e6,3e6,4e6,5e6,10e6), labels=c(0.01,0.1,1,2,3,4,5,10)) +
      #scale_y_continuous(trans="log10") +
      scale_y_continuous(limits=c(0,1.6),
         breaks=c(0,0.5,1,1.5), labels=c(0,0.5,1,1.5)) +
      theme_bw() +
      theme(axis.text.y=element_text(size=20), axis.text.x=element_text(size=16),
          axis.title=element_text(size=20),
          plot.title=element_text(size=24, face="bold"))
    outPlot <- paste0(outPrefix, "_nni-vs-",scoreVal,".pdf")
    ggsave(gp, file=outPlot, width=6, height=4)

  outFile <- paste0(outPrefix, "_summary_",scoreVal,".unfiltered.txt")
  write.table(mat.unFilt, file=outFile, col.names=T, row.names=F, quote=F, sep="\t")
  outFile <- paste0(outPrefix, "_summary_",scoreVal,".txt")
  write.table(mat, file=outFile, col.names=T, row.names=F, quote=F, sep="\t")
  outFile <- paste0(outPrefix, "_summary_",scoreVal,"_plot.txt")
  write.table(mat.plot, file=outFile, col.names=T, row.names=F, quote=F, sep="\t")

    mat[, TDP.status := NNI > TPD.positive.threshold]
    mat <- mat[Purity > minPurity & snpMAD <= maxMad]
    gp <- ggplot(mat, aes_string(x=scoreVal, y="NNI")) +
          geom_boxplot() + geom_jitter() +
          ylab("TDP Dispersion Score (NNI)") +
          stat_summary(fun.data = function(x){
              c(y = median(x)*1.05, label = length(x))},
              geom = "text", fun.y = median) +
          #geom_text(corCoeff, aes(label=V1, group=NULL), x=0.0 , y=4, hjust="left", size=3) +
          theme_bw() +
          theme(axis.text.y=element_text(size=20), axis.text.x=element_text(size=16),
              axis.title=element_text(size=20),
              plot.title=element_text(size=24, face="bold"))
    outPlot <- paste0(outPrefix, "_nni-vs-",scoreVal,"_boxplot.pdf")
    ggsave(gp, file=outPlot, width=7, height=7)
    }
}

###############################################################
###### Plot length distributions ###########
###############################################################
gp.length <- NULL
gp.interDist <- NULL
numSamples <- length(unique(segs$ID))
for (i in 1:numSamples){
  plotInd <- i
  id <- unique(segs.gain$ID)[i]
  ind <- which(segs.gain$ID == id)
  #message(id)
  segs.gain.sample <- segs.gain[ind]
    ## median dups per chromosome ##
  numSegs.per.chr <- segs.gain.sample[, .N, by=Chromosome][, median(N)]
  text.y <- 30 #max(build$data[[1]]$ymax)
  bin <- 0.05
  gp.length[[plotInd]] <- ggplot(segs.gain.sample, aes(x=length)) + theme_bw() +
          #facet_wrap(~ ID, ncol = ncol, scales = "free_y") +
          geom_histogram(aes(y=..count..), binwidth=bin, color="black", fill="grey") +
          #geom_density(aes(y=..density..), alpha=0.2, fill="red", ) +
          stat_density(aes(bin=bin, y=bin*..count..), alpha=0.2, color="red", fill="red") +
          xlab("Length (Mb)") + ylab("Number of TD Segments") +
          scale_x_continuous(trans="log10", limits=c(1e6,50e6),
              breaks=c(1e6,2e6,3e6,4e6,5e6,10e6,50e6), labels=c(1,2,3,4,5,10,50)) +
          ylim(0, text.y) +
          ggtitle(id) +
          theme(axis.text.y=element_text(size=20), axis.text.x=element_text(size=16),
                    axis.title=element_text(size=20),
                    plot.title=element_text(size=20, face="bold"))

  bin <- 0.05
  text.y <- 10
  gp.interDist[[plotInd]] <- ggplot(segs.gain.sample, aes(x=interDupDist)) + theme_bw() +
          geom_histogram(aes(y=..count..), binwidth=bin, color="black", fill="grey") +
          ggtitle(id) +
          stat_density(aes(bin=bin, y=bin*..count..), alpha=0.2, color="red", fill="red") +
          xlab("Inter-duplication Distance (Mb)") + ylab("Number of TD Segments") +
          scale_x_continuous(trans="log10", limits=c(1e6,100e6),
              breaks=c(1e6,2e6,3e6,4e6,5e6,10e6,50e6,100e6), labels=c(1,2,3,4,5,10,50,100)) +
          ylim(0, text.y) +
          theme(axis.text.y=element_text(size=20), axis.text.x=element_text(size=16),
                    axis.title=element_text(size=20),
                    plot.title=element_text(size=24, face="bold"))
}

## plot TD lengths in multipanel for all samples ##
#par(mfrow=c(nrow,ncol))
## for each page of col x row number of plots ##
#ceiling(length(gp.length) / (ncol*nrow))
ind <- c(seq(1, length(gp.length), ncol*nrow), length(gp.length))
for (i in 2:length(ind)){
  outPlot <- paste0(outPrefix, "_lengthDist_pg", i-1, ".pdf")
#  pdf(outPlot, width=plotWidth, height=plotHeight)
#  ggMultiPlot(plotlist=gp.length[(ind[i-1]):(ind[i]-1)], ncol=ncol)
#  dev.off()
  #print(gp.length[[i]])
}

## plot TD interdup distance in multipanel for all samples ##
ind <- c(seq(1, length(gp.length), ncol*nrow), length(gp.length))
for (i in 2:length(ind)){
  outPlot <- paste0(outPrefix, "_interDupDist_pg", i-1,".pdf")
#  pdf(outPlot, width=30, height=25)
#  ggMultiPlot(plotlist=gp.interDist[(ind[i-1]):(ind[i]-1)], ncol=ncol)
#  dev.off()
}
