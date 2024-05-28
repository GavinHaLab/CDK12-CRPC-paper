library(data.table)
library(ggplot2)
library(mclust)
library(dplyr)
library(stringr)
library(plyranges)
library(ape)

inFiles <- c("../TDPanalysis/wgs/UW_WGS_TDPsegs.txt")
outPrefix <- "TAN_TDPheterogeneity"

#GENOME_SIZE <- 3e9
RECIPROCAL_THRESHOLD <- 0.9
MIN_TD_LENGTH <- 1000
DUP_TD_BUFFER <- 2000

##################################################
# load all segs from all cohorts
##################################################
ptIDs <- c("05-217")
segsAll <- NULL
for (i in 1:length(inFiles)){
  cohortID <- gsub("_TDPsegs.txt", "", basename(inFiles[i]))
  segs <- fread(inFiles[i])
  if (is.null(segs$length) | sum(is.na(segs$length)) == nrow(segs)){
    segs[, length := NULL]
    segs[Chromosome == chromosome_2, length := (End - Start + 1)]
  }
  if (is.null(segs$Copy_Number)){
    segs[, Copy_Number := Copy_Number_1_2_mean]
  }
  if (is.null(segs$TD)){
    segs[!is.na(interDupDist) | interDupDist > 0, TD := TRUE]
  }
  segs <- segs[, .(ID, Chromosome, Start, End, length, interDupDist, TD, Copy_Number, NNI)]
  segs <- cbind(cohortID, segs)
  segsAll <- rbind(segsAll, segs)
}

## process segments
# filter to only samples in ptIDs
segsAll <- segsAll[ID %like% paste(ptIDs, collapse="|")]
# remove non TD segs
#segs <- segsAll[!is.na(length)]
segs <- segsAll[!is.na(length) & length >= MIN_TD_LENGTH &
               TD == TRUE &
               #!is.na(interDupDist) & interDupDist > 0 & # exclude rows that are not TDs (no length)
               !is.na(Start) & !is.na(End)
               ] 


for (id in ptIDs){
  segs.pt <- segs[ID %like% id]
  sampleIDs <- segs.pt[, unique(ID)]
  
  numSamples <- length(sampleIDs)
  
  ## get union of all TD events
  tdsAll.gr <- as(segs.pt[ID == sampleIDs[1]], "GRanges")
  for (i in 2:numSamples){
    tds.gr <- as(segs.pt[ID == sampleIDs[i]], "GRanges")
    #tdsAll.gr <- union(tdsAll.gr, tds.gr)
    tdsAll.gr <- c(tdsAll.gr, tds.gr)
    mcols(tdsAll.gr) <- NULL
    tdsAll.gr <- sort(unique(tdsAll.gr))
  }
  # remove unique TDs but allow buffer for both start and ends (maxgap for equal type overlap)
  tdsAll.dt <- as.data.table(tdsAll.gr)
  #mcols(tdsAll.gr) <- data.frame(width = width(tdsAll.gr))
  hits <- findOverlaps(tdsAll.gr, tdsAll.gr, type="equal", maxgap=DUP_TD_BUFFER)
  dupIndToRemove <- subjectHits(hits)[duplicated(queryHits(hits), fromLast = FALSE)]
  tdsAll.gr <- tdsAll.gr[-dupIndToRemove]
  
  # Build TD matrix (samples x TD events)
  tdsAll.dt <- as.data.table(tdsAll.gr)
  segTotals <- matrix(NA, nrow=1, ncol=numSamples, dimnames = list("Total_TDs", sampleIDs))
  for (i in 1:numSamples){
    tds.gr <- as(segs.pt[ID == sampleIDs[i]], "GRanges")
    ### find event-based overlap between pair ###
    # https://support.bioconductor.org/p/72656/
    hits <- findOverlaps(tdsAll.gr, tds.gr)
    overlaps <- GenomicRanges::pintersect(tdsAll.gr[queryHits(hits)], tds.gr[subjectHits(hits)], drop.nohit.range=TRUE)
    overlapPerc1 <- width(overlaps) / width(tdsAll.gr[queryHits(hits)])
    #overlapPerc2 <- width(overlaps) / width(tds.gr[subjectHits(hits)])
    hitsReciprocal <- hits[overlapPerc1 >= RECIPROCAL_THRESHOLD]# & overlapPerc2 >= RECIPROCAL_THRESHOLD]
    tdHits <- data.frame(logical(length(tdsAll.gr)))
    tdHits[queryHits(hitsReciprocal), ] <- TRUE
    colnames(tdHits) <- sampleIDs[i]
    tdsAll.dt <- cbind(tdsAll.dt, tdHits)
    segTotals[, sampleIDs[i]] <- length(tds.gr)
  }
  

  ## add column for all overlaps
  ind <- rowSums(tdsAll.dt[, 6:(numSamples+5)]) == numSamples
  tdsAll.dt[, AllTDs := FALSE]
  tdsAll.dt[ind, AllTDs := TRUE]
  
  tdsAll.dt[!ind][`05-217_LN_F_WGS` & !`05-217_LUNG_L_WGS` & !`05-217_LN_I_WGS` & !`05-217_LN_N_WGS`, .N]
  
  
  filename <- paste0(outPrefix, "_", id, "_tree.pdf")
  pdf(filename, width = 10, height = 10)
  matDist <- dist(t(tdsAll.dt[, 6:(length(sampleIDs)+5+1)]), method = "manhattan")
  matNJ <- nj(matDist)
  #matNJ$edge.length <- round(matNJ$edge.length)
  matNJ.phy <- root(matNJ, numSamples+1)
  plot(matNJ.phy, "p", use.edge.length = TRUE)
  edgelabels(sprintf("%0.2f", round(matNJ.phy$edge.length, digits=2)), bg="white", col="black", font=1)
  dev.off()
  
  
} # end for (id in ptIDs)
