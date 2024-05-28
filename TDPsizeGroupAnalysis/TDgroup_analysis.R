library(data.table)
library(ggplot2)
library(mclust)
library(dplyr)

inFiles <- list.files("../TDPanalysis/", pattern = "TDPsegs.txt", full.names=TRUE, recursive=TRUE)
outFile <- "TD_groupClass.txt"

##################################################
# load all segs from all cohorts
##################################################
mat <- NULL
for (i in 1:length(inFiles)){
  cohortID <- gsub("_TDPsegs.txt", "", basename(inFiles[i]))
  segs <- fread(inFiles[i])
  if (is.null(segs$length) | sum(is.na(segs$length)) == nrow(segs)){
    segs[, length := NULL]
    segs[Chromosome == chromosome_2, length := (End - Start + 1)]
  }
  segs <- segs[, .(ID, Chromosome, Start, End, length, interDupDist, NNI)]
  segs <- cbind(cohortID, segs)
  mat <- rbind(mat, segs)
}

##################################################
## process segments
##################################################
mat <- mat[!is.na(length) & length >= 1 & !is.na(interDupDist)] # exclude rows that are not TDs (no length)
mat[, length := length / 1000]

## get sample-cohort key
sampleCohortKey <- unique(mat[, .(cohortID, ID)])
setkey(sampleCohortKey, "ID")
## hardcoded
sampleCohortKey[cohortID == "UW_WES" | cohortID == "ECDT", seqType := "WES"]
sampleCohortKey[cohortID == "hmf" | cohortID == "WCDT" | cohortID == "UW_WGS", seqType := "WGS"]

##################################################
## function for classifying TD group classes
##################################################
## define length range for each TD group/class - in log10 scale
grp1 <- log10(c(0.001, 100))
grp2 <- log10(c(100.001, 1000))
grp3 <- log10(c(1000.001, 10000))
# grp4 = grp1 + grp2
# grp5 = grp1 + grp3
# grp6 = grp2 + grp3

classifyTDgroup <- function(x){
  tdClass <- vector("integer", length(x))
  # check if mean length falls in each class
  tdClass[dplyr::between(x, grp1[1], grp1[2])] <- 1
  tdClass[dplyr::between(x, grp2[1], grp2[2])] <- 2
  tdClass[dplyr::between(x, grp3[1], grp3[2])] <- 3
  
  # check for mixture groups (grp4-6)
  groupClass <- NULL
  if (sum(tdClass %in% c(2,3)) == 2){
    groupClass <- 6
  }else if (sum(tdClass %in% c(1,3)) == 2){
    groupClass <- 5
  }else if (sum(tdClass %in% c(1,2)) == 2){
    groupClass <- 4
  }else{
    groupClass <- tdClass
  }
  
  return(list(groupClass = groupClass, tdClass = tdClass))
}


##################################################
## perform GMM using densityMclust
##################################################
maxVar <- 0.2
minPro <- 0.15 # minimum mixed weights
tdMeans <- list()
tdWeights <- list()
for (i in unique(mat$ID)){
  clust1 <- densityMclust(log10(mat[ID==i, length]), G=1:4, plot = FALSE)
  # start with components having minimum mixed weights 
  ind <- clust1$parameters$pro > minPro
  if (length(clust1$parameters$mean) == 1){ # if only one component
    ind <- 1
  # otherwise check and max variance 
  }else if (length(clust1$parameters$variance$sigmasq) == length(ind)){ 
    ind <- ind & clust1$parameters$variance$sigmasq < maxVar
  }
  tdMeans <- c(tdMeans, list(clust1$parameters$mean[ind]))
  tdWeights <- c(tdWeights, list(clust1$parameters$pro[ind]))
}
names(tdMeans) <- unique(mat$ID)
names(tdWeights) <- unique(mat$ID)
#tdMeans <- mat[, .(list(densityMclust(log10(length), G=1:4, plot=FALSE)$parameters$mean)), by=ID]

##################################################
## classify TD calls
##################################################
groupCalls <- lapply(tdMeans, classifyTDgroup)
groupCalls <- t(do.call('cbind', groupCalls))
groupCalls <- data.table(ID = rownames(groupCalls), 
                         groupClass = as.character(unlist(groupCalls[,"groupClass"])),
                         tdClassMix = plyr::ldply(groupCalls[, "tdClass"], paste, collapse=",")$V1)
groupCalls[, cohortID := sampleCohortKey[groupCalls$ID, cohortID]]
groupCalls[, seqType := sampleCohortKey[groupCalls$ID, seqType]]
setkey(groupCalls, ID)
groupCalls[, tdGroupMeans := {
  lapply(ID, function(x){ paste(sprintf(tdMeans[[x]], fmt='%.3f'), collapse=",") })
}]

##################################################
## plot distribution of TD segment lengths by cohort
##################################################
mat$cohortID <- factor(mat$cohortID, c("ECDT", "TAN_WES", "TAN_WGS", "WCDT", "hmf"))
gp <- ggplot(mat, aes(log10(length), fill=cohortID)) +
  geom_density(color="black", linewidth=1, alpha=0.75) +
  scale_fill_manual(values = list(ECDT="#E69F00", TAN_WES="#E69F00", TAN_WGS="#0072B2", WCDT="#0072B2", hmf="#0072B2")) +
  xlab("TD length (log10)") + xlim(c(-1,6)) + ylab("Density") +
  facet_grid(~cohortID) +
  theme_bw() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        strip.background = element_blank())
ggplot2::ggsave("TD_cohort_segSizeDist.pdf", gp, width=8, height=3)

##################################################
## plot distribution of TD segment lengths by sequencing type
mat$cohortID <- factor(mat$cohortID, c("ECDT", "TAN_WES", "TAN_WGS", "WCDT", "hmf"))
mat[, seqType := sampleCohortKey[mat$ID, seqType]]
gp <- ggplot(mat, aes(log10(length), fill=seqType)) +
  geom_density(color="black", linewidth=1, alpha=0.75) +
  scale_fill_manual(values = list(WES="#E69F00", WGS="#0072B2")) +
  xlab("TD length (log10)") + xlim(c(-1,6)) + ylab("Density") +
  facet_wrap(~seqType, scale="free_y") +
  theme_bw() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        strip.background = element_blank())
ggplot2::ggsave(filename="TD_segSizeDist.pdf", plot=gp, width=6, height=3)

##################################################
## plot TD group classes (tdClass, e.g. 2,3)
groupCalls$groupClass <- factor(groupCalls$groupClass, c("6", "3", "2", "5", "4", "1"))
groupCalls$tdClassMix <- factor(groupCalls$tdClassMix, c("2,3", "3", "2", "1,3","1,2", "1"))
groupCalls$cohortID <- factor(groupCalls$cohortID, c("ECDT", "TAN_WES", "TAN_WGS", "WCDT", "hmf"))
gp <- ggplot(groupCalls, aes(groupClass, fill=cohortID)) +
  #geom_histogram(stat="count") +
  geom_bar(aes(y= after_stat(count/sum(count)))) + 
  scale_fill_manual(values = list(ECDT="#E69F00", TAN_WES="#E69F00", TAN_WGS="#0072B2", WCDT="#0072B2", hmf="#0072B2")) +
  xlab("TD Group") + ylab("Count") + ylim(0,0.4) + 
  facet_wrap(~cohortID, nrow=1) + 
  theme_bw() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        strip.background = element_blank())
ggplot2::ggsave("TD_groupClass.pdf", plot=gp, width=8, height=3)


##################################################
## plot TD means
wesMeans <- unlist(tdMeans[sampleCohortKey[seqType=="WES", ID]])
wgsMeans <- unlist(tdMeans[sampleCohortKey[seqType=="WGS", ID]])
meansToPlot <- data.table(c(wesMeans, wgsMeans))
meansToPlot[, seqType := c(rep("WES", length(wesMeans)), rep("WGS", length(wgsMeans)))]
gp <- ggplot(meansToPlot, aes(V1, fill=seqType)) +
  geom_density(color="black", linewidth=1, alpha=0.75) +
  scale_fill_manual(values = list(WES="#E69F00", WGS="#0072B2")) +
  xlab("Mean TD class length (log10)") + xlim(c(-1,6)) + ylab("Density") +
  facet_wrap(~seqType, scales="free_y") +
  theme_bw() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        strip.background = element_blank())
ggplot2::ggsave("TD_meanDistBySeqType.pdf", plot=gp, width=6, height=3)

##################################################
## plot TD length + predicted TD means overlaid together
gp <- ggplot() +
  geom_histogram(aes(log10(length), y=stat(density), fill=seqType), data=mat, binwidth=0.25, color="black", size=0, alpha=0.75) +
  geom_density(aes(V1, color=seqType), data=meansToPlot, size=0.75, linetype="twodash", color="red") +
  scale_fill_manual(values = list(WES="#E69F00", WGS="#0072B2")) +
  #scale_color_manual(values = list(WES="#E69F00", WGS="#0072B2")) +
  xlab("TD length (log10)") + xlim(c(-1,6)) + ylab("Density") +
  facet_wrap(~seqType, scale="free_y") +
  theme_bw() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        strip.background = element_blank())
ggplot2::ggsave("TD_segSizeAndTDmeansDist.pdf", plot=gp, width=6, height=3)

## print out TD group output file
fwrite(groupCalls, file=outFile, sep = "\t", quote=FALSE)

