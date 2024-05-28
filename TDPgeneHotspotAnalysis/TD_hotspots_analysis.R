library(data.table)
library(ggplot2)
library(mclust)
library(dplyr)
library(GenomicRanges)
library(ggbio)
library(karyoploteR)
library(stringr)
library(ggrepel)

inFiles <- list.files("../TDPanalysis", pattern = "TDPsegs_allSamples.txt", full.names = TRUE, recursive = TRUE)
tdpSampleList <- fread("TDP_sampleList.txt", header=FALSE)
geneFile <- "Census_all_hg38_01_27_19_DL06Jul2021.txt" #hg38
geneFile.hg19 <- "Census_allWedFeb700_28_332024_HG19.tsv"
outPrefix <- "TD_hotspots"
binWidth <- 50000
cohorts.hg38 <- c("ECDT", "WCDT", "UW_WES", "UW_WGS")
cohorts.hg19 <- c("hmf")
#TOTAL_COHORT_SIZE <- 983

##################################################
# load all segs from all cohorts
##################################################
mat <- NULL
for (i in 1:length(inFiles)){
  cohortID <- gsub("_TDPsegs_allSamples.txt", "", basename(inFiles[i]))
  segs <- fread(inFiles[i])
  if (is.null(segs$length) | sum(is.na(segs$length)) == nrow(segs)){
    segs[, length := NULL]
    segs[Chromosome == chromosome_2, length := (End - Start + 1)]
  }
  segs <- segs[, .(ID, Chromosome, Start, End, length, interDupDist, NNI)]
  segs <- cbind(cohortID, segs)
  mat <- rbind(mat, segs)
}
numSamples <- length(unique(mat$ID))
TOTAL_COHORT_SIZE <- numSamples
TDP_SAMPLE_SIZE <- dim(tdpSampleList)[1]
NON_TDP_SAMPLE_SIZE <- TOTAL_COHORT_SIZE - TDP_SAMPLE_SIZE
##################################################
## process segments
##################################################
#seqlevelsStyle(mat$Chromosome) <- "NCBI" #stopped working for some reason
mat.all <- copy(mat)
mat[, Chromosome := gsub("chr","",Chromosome)]
mode(mat$Start) <- "integer"
mode(mat$End) <- "integer"
mat <- mat[!is.na(length) & length >= 1 & !is.na(interDupDist) &
             !is.na(Start) & !is.na(End)] # exclude rows that are not TDs (no length)
mat[, Chromosome := paste0("chr", Chromosome)]

## indicate which sample is TDP
mat[, isTDP := FALSE]
mat[ID %in% tdpSampleList$V1, isTDP := TRUE]
## convert to GRanges object
mat.gr <- as(mat, "GRanges")
seqlevelsStyle(mat.gr) <- "UCSC"
  


##################################################
## get genome tiled bins
seqinfo <- Seqinfo(genome="GRCh38.p13")[c(1:22,"X","Y")] # use NCBI style for ease of use now.
bins.gr <- tileGenome(seqlengths = seqlengths(seqinfo), tilewidth = binWidth, cut.last.tile.in.chrom = TRUE)
seqlevelsStyle(bins.gr) <- "UCSC"
bins.nonTDP.gr <- copy(bins.gr)
bins.dt <- as.data.table(bins.gr)
bins.nonTDP.dt <- as.data.table(bins.nonTDP.gr)

##################################################
## function to load and annotation genes
loadCGCgeneList <- function(fileName){
  genes <- fread(fileName)
  genes <- genes[, .(`Gene Symbol`, `Genome Location`, `Tier`, `Role in Cancer`)]
  
  chrPosn <- genes[, do.call('rbind', ((str_split(`Genome Location`, pattern=":|-"))))]
  colnames(chrPosn) <- c("chr", "start", "stop")
  genes <- cbind(chrPosn, genes)
  genes <- genes[start != "" & stop != ""]
  genes[, chr := paste0("chr", chr)]
  genes.gr <- as(genes, "GRanges")
  seqlevelsStyle(genes.gr) <- "UCSC"
  
  #genes[, type := NULL]
  genes.oncogene <- genes[`Role in Cancer` %like% "oncogene" | `Role in Cancer` == "fusion"]
  genes.tsg <- genes[`Role in Cancer` %like% "TSG" | `Role in Cancer` == "fusion"]
  
  genes.oncogene.gr <- as(genes.oncogene, "GRanges")
  seqlevelsStyle(genes.oncogene.gr) <- "UCSC"
  genes.tsg.gr <- as(genes.tsg, "GRanges")
  seqlevelsStyle(genes.tsg.gr) <- "UCSC"
  
  return(list(genes=genes, genes.gr=genes.gr, genes.oncogene=genes.oncogene, 
              genes.oncogene.gr=genes.oncogene.gr, genes.tsg=genes.tsg, 
              genes.tsg.gr=genes.tsg.gr))
}

genes <- loadCGCgeneList(geneFile)
genes.hg19 <- loadCGCgeneList(geneFile.hg19)
genes.nonTDP <- loadCGCgeneList(geneFile)
genes.hg19.nonTDP <- loadCGCgeneList(geneFile.hg19)

####################################################
## BIN LEVEL ANALYSIS ##############################
####################################################
## assign genes to bins
hits <- findOverlaps(query = bins.gr, subject = genes$genes.gr)
bins.dt[queryHits(hits), gene := genes$genes[subjectHits(hits), `Gene Symbol`]]
bins.nonTDP.dt[queryHits(hits), gene := genes$genes[subjectHits(hits), `Gene Symbol`]]

##################################################
## count overlaps of TDs for each bin - TDP SAMPLES
tdCounts <- countOverlaps(query = bins.gr, subject = mat.gr[mat$isTDP == TRUE])
bins.dt[, tdCounts := tdCounts]

## count overlaps of TDs for each bin - NON-TDP SAMPLES
tdCounts.nonTDP <- countOverlaps(query = bins.nonTDP.gr, subject = mat.gr[mat$isTDP == FALSE])
bins.nonTDP.dt[, tdCounts := tdCounts.nonTDP]


## compute fisher's test
bins.dt[, tdCounts_nonTDP := bins.nonTDP.dt$tdCounts]
bins.dt[, TDP_notGain := TDP_SAMPLE_SIZE - tdCounts]
bins.dt[, NotTDP_notGain := NON_TDP_SAMPLE_SIZE - tdCounts_nonTDP]
bins.tdCount.mat <- cbind(TDP_Gain = bins.dt$tdCounts, 
                          TDP_notGain = bins.dt$TDP_notGain, 
                          NotTDP_Gain = bins.dt$tdCounts_nonTDP, 
                          NotTDP_notGain = bins.dt$NotTDP_notGain)
chiSqTest.p <- apply(bins.tdCount.mat, 1, function(x) { chisq.test(t(matrix(x, nrow=2)))$p.value })
chiSqTest.q <- p.adjust(chiSqTest.p, method = "bonferroni")
bins.dt[, chiSqTest.p := chiSqTest.p]
bins.dt[, chiSqTest.q := chiSqTest.q]

##################################################
## plot tdCounts by bin with ideograms from karyoplotR - TDP SAMPLES ##
#maxCount <- max(bins.dt$tdCounts)
#mcols(bins.gr)$y <- bins.dt$tdCounts / maxCount * 1.2 # offset value
ind <- bins.dt[chiSqTest.q < 0.001, which = TRUE]
yScaleFactor <- 0.5

plotFile <- "TD_hotspot_bins.pdf"
pdf(plotFile, width=14, height=2.5)
kp <- plotKaryotype(plot.type = 4, main = "") 
kpDataBackground(kp, col="white")
kpBars(kp, bins.gr, y1 = bins.dt[, tdCounts/TDP_SAMPLE_SIZE/yScaleFactor], border="grey")
kpBars(kp, bins.gr, y1 = bins.nonTDP.dt[, tdCounts/NON_TDP_SAMPLE_SIZE/yScaleFactor], border="black")
kpBars(kp, bins.gr[ind], y1 = bins.dt[ind, tdCounts/TDP_SAMPLE_SIZE/yScaleFactor], border="red")
#ind <- bins.gr$y > 15
#kpBars(kp, bins.gr, y1 = bins.gr$y[ind]/maxCount, border="black")
#kpText(kp, data = bins.gr[ind], labels = bins.dt[ind, gene], col="black", srt=45)
kpAxis(kp, data.panel = 1, ymin=0, ymax=yScaleFactor)
dev.off()

## plot tdCounts by bin with ideograms from karyoplotR - NON-TDP SAMPLES ##
maxCount <- max(bins.nonTDP.dt$tdCounts)
#mcols(bins.nonTDP.gr)$y <- bins.nonTDP.dt$tdCounts / maxCount * 1.2 # offset value
ind <- bins.nonTDP.dt[tdCounts > 0, which = TRUE]

plotFile <- "TD_hotspot_bins_nonTDP.pdf"
pdf(plotFile, width=14, height=3)
kp <- plotKaryotype(plot.type = 4, main = "") 
kpDataBackground(kp)
kpBars(kp, bins.nonTDP.gr[ind], y1 = bins.nonTDP.dt[ind, tdCounts/NON_TDP_SAMPLE_SIZE], border="grey")
kpAxis(kp, data.panel = 1)#, ymin=0, ymax=1)
dev.off()

## print out bin-level file
outFile <- "TD_hotspots_bins.txt"
fwrite(bins.dt, outFile, row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")

####################################################
## GENE LEVEL ANALYSIS #############################
####################################################


##################################################
## count oncogenes found fully within TDs -- TDP SAMPLES
## hg38 - UW TAN, SU2C, WCDT
ind.hg38.TDP <- mat[cohortID %in% cohorts.hg38 & isTDP == TRUE, which=TRUE]
## hg19 - HMF
ind.hg19.TDP <- mat[cohortID %in% cohorts.hg19 & isTDP == TRUE, which=TRUE]

tdCounts <- countOverlaps(query = genes$genes.oncogene.gr, subject = mat.gr[ind.hg38.TDP], type="within")
genes$genes.oncogene[, tdCounts := NULL]
genes$genes.oncogene[, tdCounts := tdCounts]
setkey(genes$genes.oncogene, `Gene Symbol`)

hits.hg19 <- as.data.table(findOverlaps(query = genes.hg19$genes.oncogene.gr, subject = mat.gr[ind.hg19.TDP], type="within"))
hits.hg19 <- hits.hg19[, .N, by=queryHits]
genes$genes.oncogene[genes.hg19$genes.oncogene[hits.hg19$queryHits, `Gene Symbol`], 
                     tdCounts := tdCounts + hits.hg19$N]
genes$genes.oncogene <- genes$genes.oncogene[genes$genes.oncogene.gr$`Gene Symbol`] # re-order by genes in GR

##################################################
## count oncogenes found fully within TDs -- NON-TDP SAMPLES
## hg38 - UW TAN, SU2C, WCDT
ind.hg38.nonTDP <- mat[cohortID %in% cohorts.hg38 & isTDP == FALSE, which=TRUE]
## hg19 - HMF
ind.hg19.nonTDP <- mat[cohortID %in% cohorts.hg19 & isTDP == FALSE, which=TRUE]

tdCounts <- countOverlaps(query = genes$genes.oncogene.gr, subject = mat.gr[ind.hg38.nonTDP], type="within")
genes.nonTDP$genes.oncogene[, tdCounts := NULL]
genes.nonTDP$genes.oncogene[, tdCounts := tdCounts]
setkey(genes.nonTDP$genes.oncogene, `Gene Symbol`)

hits.hg19 <- as.data.table(findOverlaps(query = genes.hg19.nonTDP$genes.oncogene.gr, subject = mat.gr[ind.hg19.nonTDP], type="within"))
hits.hg19 <- hits.hg19[, .N, by=queryHits]
genes.nonTDP$genes.oncogene[genes.hg19.nonTDP$genes.oncogene[hits.hg19$queryHits, `Gene Symbol`], 
                     tdCounts := tdCounts + hits.hg19$N]
genes.nonTDP$genes.oncogene <- genes.nonTDP$genes.oncogene[genes.nonTDP$genes.oncogene.gr$`Gene Symbol`] # re-order by genes in GR

####################################################
## ONOCOGENES ##
## Test gene enrichment (fisher's exact test) for TDP vs non-TDP
geneMat <- cbind(TDP_Gain = genes$genes.oncogene$tdCounts, 
                 TDP_notGain = TDP_SAMPLE_SIZE - genes$genes.oncogene$tdCounts, 
                 NotTDP_Gain = genes.nonTDP$genes.oncogene$tdCounts, 
                 NotTDP_notGain = NON_TDP_SAMPLE_SIZE - genes.nonTDP$genes.oncogene$tdCounts)
fisherTest.p <- apply(geneMat, 1, function(x) { fisher.test(t(matrix(x, nrow=2)), alternative="two.sided")$p.value })
fisherTest.q <- p.adjust(fisherTest.p, method = "bonferroni")
fisherTest.OR <- apply(geneMat, 1, function(x) { fisher.test(t(matrix(x, nrow=2)), alternative="two.sided")$estimate })
groupPropDiff <- apply(geneMat, 1, function(x) { (x[1]/x[2]) - (x[3]/x[4]) })
#propTest.p <- apply(geneMat, 1, function(x) { prop.test(t(matrix(x, ncol=2)))$p.value })
#propTest.q <- p.adjust(propTest.p, method = "bonferroni")

genes$genes.oncogene[, fisherTest.p := fisherTest.p]
genes$genes.oncogene[, fisherTest.q := fisherTest.q]
genes$genes.oncogene[, fisherTest.sig := FALSE]
genes$genes.oncogene[fisherTest.q < 0.01, fisherTest.sig := TRUE]
genes$genes.oncogene[, fisherTest.OR := fisherTest.OR]
genes$genes.oncogene[, fisherTest.groupPropDiff := groupPropDiff]

## VOLCANO PLOT: Log Odds Ratio
ind <- genes$genes.oncogene[-log10(fisherTest.p) > 6 & log(fisherTest.OR) > 2, which=TRUE]
genes$genes.oncogene[, genesToLabel := ""]
genes$genes.oncogene[ind, genesToLabel := `Gene Symbol`]
genes$genes.oncogene[`Gene Symbol` %in% c("AR","IDH2","CCND1","KRAS","ETV5","CDK4","GATA2","FGF4","IDH2"), genesToLabel := `Gene Symbol`]
gp <- ggplot(genes$genes.oncogene, aes(x=log(fisherTest.OR), y=-log10(fisherTest.p))) +
  geom_point(aes(color=fisherTest.sig)) +
  scale_color_manual(values=c(`TRUE`="red", `FALSE`="lightgrey")) +
  xlab("Log odds ratio") + xlim(-1,4) +
  #scale_y_continuous(name="-log_10(p-value)", breaks=c(0,2,4,6), labels=c("0.0","2.0","4.0","6.0")) + 
  ylab("-log_10(p-value)") +
  geom_text_repel(aes(label = genesToLabel), box.padding = 0.50, max.overlaps = 100, force = 15) +
  theme_bw() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        strip.background = element_blank())
ggplot2::ggsave("TD_Volcano_logOR_oncogene.pdf", gp, width=5.25, height=4)

##################################################
## plot tdCounts by oncogene ideograms from karyoplotR ##
## overlap of gene fully within TD event - TDP SAMPLES
maxCount <- max(genes$genes.oncogene$tdCounts)
mcols(genes$genes.oncogene.gr)$y <- genes$genes.oncogene$tdCounts / maxCount * 1.10 # offset value for gene labels
ind <- genes$genes.oncogene[fisherTest.sig == TRUE, which = TRUE]
ind.labels <- genes$genes.oncogene[fisherTest.sig == TRUE | `Gene Symbol` %in% c("AR"), which = TRUE]
yScaleFactor <- 0.5
plotFile <- "TD_hotspot_oncogenes.pdf"
pdf(plotFile, width=14, height=2.5)
kp <- plotKaryotype(plot.type = 4, main = "")
kpDataBackground(kp, col="white")
kpBars(kp, genes$genes.oncogene.gr, y1 = genes$genes.oncogene[, tdCounts/TDP_SAMPLE_SIZE/yScaleFactor], border="lightgrey")
kpBars(kp, genes$genes.oncogene.gr[ind], y1 = genes$genes.oncogene[ind, tdCounts/TDP_SAMPLE_SIZE/yScaleFactor], border="red")
#ind <- genes.gr$y > 15
#kpBars(kp, genes.gr, y1 = genes.gr$y[ind]/maxCount, border="black")
kpText(kp, data = genes$genes.oncogene.gr[ind.labels], labels = genes$genes.oncogene[ind.labels, `Gene Symbol`], col="black", srt=45)
kpAxis(kp, data.panel = 1, ymin=0, ymax=yScaleFactor)
dev.off()

## overlap of gene fully within TD event - NON-TDP SAMPLES
maxCount <- max(genes.nonTDP$genes.oncogene$tdCounts)
mcols(genes.nonTDP$genes.oncogene.gr)$y <- genes.nonTDP$genes.oncogene$tdCounts / maxCount * 1.10 # offset value for gene labels
ind <- genes$genes.oncogene[fisherTest.sig == TRUE, which = TRUE]
ind.labels <- genes.nonTDP$genes.oncogene[tdCounts/NON_TDP_SAMPLE_SIZE > 0.15, which=TRUE]
yScaleFactor <- 0.5
plotFile <- "TD_hotspot_oncogenes_noTDP.pdf"
pdf(plotFile, width=14, height=2.5)
kp <- plotKaryotype(plot.type = 4, main = "")
kpDataBackground(kp, col="white")
kpBars(kp, genes.nonTDP$genes.oncogene.gr, y1 = genes.nonTDP$genes.oncogene[, tdCounts/NON_TDP_SAMPLE_SIZE/yScaleFactor], border="lightgrey")
kpBars(kp, genes.nonTDP$genes.oncogene.gr[ind], y1 = genes.nonTDP$genes.oncogene[ind, tdCounts/NON_TDP_SAMPLE_SIZE/yScaleFactor], border="red")
#ind <- genes.gr$y > 15
#kpBars(kp, genes.gr, y1 = genes.gr$y[ind]/maxCount, border="black")
kpText(kp, data = genes.nonTDP$genes.oncogene.gr[ind.labels], y = genes.nonTDP$genes.oncogene[ind.labels, tdCounts/NON_TDP_SAMPLE_SIZE/yScaleFactor],
       labels = genes.nonTDP$genes.oncogene[ind.labels, `Gene Symbol`], col="black", srt=45)
kpAxis(kp, data.panel = 1, ymin=0, ymax=yScaleFactor)
dev.off()


##################################################
## count tsg found broken by TDs
mat.start <- mat
mat.start <- mat.start[, End := Start]
mat.end <- mat
mat.end <- mat.end[, Start := End]
mat.start.gr <- as(mat.start, "GRanges")
mat.end.gr <- as(mat.end, "GRanges")

######### for TDP samples
## hg38 - UW TAN, SU2C, WCDT
hits.start <- findOverlaps(query = genes$genes.tsg.gr, subject = mat.start.gr[ind.hg38.TDP], type="any")
hits.end <- findOverlaps(query = genes$genes.tsg.gr, subject = mat.end.gr[ind.hg38.TDP], type="any")
hits <- unique(rbind(as.data.table(hits.start), as.data.table(hits.end)))
## hg19 - HMF
hits.start <- findOverlaps(query = genes$genes.tsg.gr, subject = mat.start.gr[ind.hg19.TDP], type="any")
hits.end <- findOverlaps(query = genes$genes.tsg.gr, subject = mat.end.gr[ind.hg19.TDP], type="any")
hits.hg19 <- unique(rbind(as.data.table(hits.start), as.data.table(hits.end)))

tdCounts.tsg <- hits[,.N, by=queryHits] 
tdCounts.tsg.hg19 <- hits.hg19[,.N, by=queryHits] 
genes$genes.tsg[, tdCounts := 0]
genes$genes.tsg[tdCounts.tsg$queryHits, tdCounts := tdCounts.tsg$N]
genes$genes.tsg[tdCounts.tsg.hg19$queryHits, tdCounts := tdCounts + tdCounts.tsg.hg19$N]


######## for non TDP samples ##
## hg38 - UW TAN, SU2C, WCDT
hits.start <- findOverlaps(query = genes$genes.tsg.gr, subject = mat.start.gr[ind.hg38.nonTDP], type="any")
hits.end <- findOverlaps(query = genes$genes.tsg.gr, subject = mat.end.gr[ind.hg38.nonTDP], type="any")
hits <- unique(rbind(as.data.table(hits.start), as.data.table(hits.end)))
## hg19 - HMF
hits.start <- findOverlaps(query = genes$genes.tsg.gr, subject = mat.start.gr[ind.hg19.nonTDP], type="any")
hits.end <- findOverlaps(query = genes$genes.tsg.gr, subject = mat.end.gr[ind.hg19.nonTDP], type="any")
hits.hg19 <- unique(rbind(as.data.table(hits.start), as.data.table(hits.end)))

tdCounts.nonTDP.tsg <- hits[,.N, by=queryHits] 
tdCounts.nonTDP.tsg.hg19 <- hits.hg19[,.N, by=queryHits] 
genes.nonTDP$genes.tsg[, tdCounts := 0]
genes.nonTDP$genes.tsg[tdCounts.nonTDP.tsg$queryHits, tdCounts := tdCounts.nonTDP.tsg$N]
genes.nonTDP$genes.tsg[tdCounts.nonTDP.tsg$queryHits, tdCounts := tdCounts + tdCounts.nonTDP.tsg$N]

###############################
## TUMOR SUPPRESSOR GENES ##
## Test gene enrichment (fisher's exact test) for TDP vs non-TDP
geneMat <- cbind(TDP_Gain = genes$genes.tsg$tdCounts, 
                 TDP_notGain = TDP_SAMPLE_SIZE - genes$genes.tsg$tdCounts, 
                 NotTDP_Gain = genes.nonTDP$genes.tsg$tdCounts, 
                 NotTDP_notGain = NON_TDP_SAMPLE_SIZE - genes.nonTDP$genes.tsg$tdCounts)
fisherTest.p <- apply(geneMat, 1, function(x) { fisher.test(t(matrix(x, nrow=2)), alternative="two.sided")$p.value })
fisherTest.q <- p.adjust(fisherTest.p, method = "bonferroni")
fisherTest.OR <- apply(geneMat, 1, function(x) { fisher.test(t(matrix(x, nrow=2)), alternative="two.sided")$estimate })
groupPropDiff <- apply(geneMat, 1, function(x) { (x[1]/x[2]) - (x[3]/x[4]) })
#propTest.p <- apply(geneMat, 1, function(x) { prop.test(t(matrix(x, ncol=2)))$p.value })
#propTest.q <- p.adjust(propTest.p, method = "bonferroni")

genes$genes.tsg[, fisherTest.p := fisherTest.p]
genes$genes.tsg[, fisherTest.q := fisherTest.q]
genes$genes.tsg[, fisherTest.sig := FALSE]
genes$genes.tsg[fisherTest.q < 0.01, fisherTest.sig := TRUE]
genes$genes.tsg[, fisherTest.OR := fisherTest.OR]
genes$genes.tsg[, fisherTest.groupPropDiff := groupPropDiff]

## TSG VOLCANO PLOT: Log Odds Ratio
ind <- genes$genes.tsg[!is.infinite(fisherTest.OR), which =TRUE]
ind.labels <- genes$genes.tsg[-log10(fisherTest.p) > 5 & log(fisherTest.OR) > 2, which=TRUE]
genes$genes.tsg[, genesToLabel := ""]
genes$genes.tsg[ind.labels, genesToLabel := `Gene Symbol`]
genes$genes.tsg[`Gene Symbol` %in% c("AR","MUC1","IL7R","IDH2",""), genesToLabel := `Gene Symbol`]
gp <- ggplot(genes$genes.tsg[ind], aes(x=log(fisherTest.OR), y=-log10(fisherTest.p))) +
  geom_point(aes(color=fisherTest.sig)) +
  scale_color_manual(values=c(`TRUE`="blue", `FALSE`="lightgrey")) +
  xlab("Log odds ratio") + xlim(-1,5) +
  scale_y_continuous(name="-log_10(p-value)", breaks=c(0,3,6,9), labels=c("0.0","3.0","6.0","9.0")) + #ylab("-log_10(p-value)") +
  geom_text_repel(aes(label = genesToLabel), box.padding = 0.10, max.overlaps = Inf, force = 5) +
  theme_bw() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        strip.background = element_blank())
ggplot2::ggsave("TD_Volcano_logOR_tsg.pdf", gp, width=5, height=4)



##################################################
## plot tdCounts by oncogene ideograms from karyoplotR ##

## breaking of gene fully by TD event - TDP SAMPLES
maxCount <- max(genes$genes.tsg$tdCounts)
mcols(genes$genes.tsg.gr)$y <- genes$genes.tsg$tdCounts / maxCount * 1.10 # offset value for gene labels
ind <- genes$genes.tsg[fisherTest.sig == TRUE & !is.infinite(fisherTest.OR), which = TRUE]
#ind.labels <- genes$genes.tsg[tdCounts/NON_TDP_SAMPLE_SIZE > 0.01, which=TRUE]
yScaleFactor <- 0.5
plotFile <- "TD_hotspot_tsg.pdf"
pdf(plotFile, width=14, height=2.5)
kp <- plotKaryotype(plot.type = 4, main = "")
kpDataBackground(kp, col="white")
kpBars(kp, genes$genes.tsg.gr, y1 = genes$genes.tsg[, tdCounts/TDP_SAMPLE_SIZE/yScaleFactor], border="lightgrey")
kpBars(kp, genes$genes.tsg.gr[ind], y1 = genes$genes.tsg[ind, tdCounts/TDP_SAMPLE_SIZE/yScaleFactor], border="blue")
#ind <- genes.gr$y > 15
#kpBars(kp, genes.gr, y1 = genes.gr$y[ind]/maxCount, border="black")
kpText(kp, data = genes$genes.tsg.gr[ind], y = genes$genes.tsg[ind, tdCounts/TDP_SAMPLE_SIZE/yScaleFactor],
       labels = genes$genes.tsg[ind, `Gene Symbol`], col="black", srt=45)
kpAxis(kp, data.panel = 1, ymin=0, ymax=yScaleFactor)
dev.off()



## breaking of gene fully by TD event - NON-TDP SAMPLES
maxCount <- max(genes.nonTDP$genes.tsg$tdCounts)
mcols(genes.nonTDP$genes.tsg.gr)$y <- genes.nonTDP$genes.tsg$tdCounts / maxCount * 1.10 # offset value for gene labels
ind <- genes$genes.tsg[fisherTest.sig == TRUE & !is.infinite(fisherTest.OR), which = TRUE]
#ind.labels <- genes.nonTDP$genes.tsg[tdCounts/NON_TDP_SAMPLE_SIZE > 0.02, which=TRUE]
yScaleFactor <- 0.5
plotFile <- "TD_hotspot_tsg_noTDP.pdf"
pdf(plotFile, width=14, height=2.5)
kp <- plotKaryotype(plot.type = 4, main = "")
kpDataBackground(kp, col="white")
kpBars(kp, genes.nonTDP$genes.tsg.gr, y1 = genes.nonTDP$genes.tsg[, tdCounts/NON_TDP_SAMPLE_SIZE/yScaleFactor], border="lightgrey")
kpBars(kp, genes.nonTDP$genes.tsg.gr[ind], y1 = genes.nonTDP$genes.tsg[ind, tdCounts/NON_TDP_SAMPLE_SIZE/yScaleFactor], border="blue")
#ind <- genes.gr$y > 15
#kpBars(kp, genes.gr, y1 = genes.gr$y[ind]/maxCount, border="black")
kpText(kp, data = genes.nonTDP$genes.tsg.gr[ind], y = genes.nonTDP$genes.tsg[ind, tdCounts/NON_TDP_SAMPLE_SIZE/yScaleFactor],
       labels = genes.nonTDP$genes.tsg[ind, `Gene Symbol`], col="black", srt=45)
kpAxis(kp, data.panel = 1, ymin=0, ymax=yScaleFactor)
dev.off()


## print out gene-level file - oncogenes
genes$genes.oncogene[, tdCounts_nonTDP := genes.nonTDP$genes.oncogene$tdCounts]
outFile <- "TD_hotspots_oncogenes.txt"
fwrite(genes$genes.oncogene, outFile, row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
#outFile <- "TD_hotspots_oncogenes_nonTDP.txt"
#fwrite(genes.nonTDP$genes.oncogene, outFile, row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")

outFile <- "TD_hotspots_tsg.txt"
fwrite(genes$genes.tsg, outFile, row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
outFile <- "TD_hotspots_tsg_nonTDP.txt"
fwrite(genes.nonTDP$genes.tsg, outFile, row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")

