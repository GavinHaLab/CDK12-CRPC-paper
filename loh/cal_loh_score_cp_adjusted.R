#####
#####Percent of genome affected by LOH
#####
# With 3 types of LOH definition
# 	1) Delete greater than 90%/75% (90% stands for arm level loss)
# 	2) Delete greater than 90%/75% & only keep greater than 10MB
#   3) Delete greater than 90%/75% & only keep greater than 15MB
#
# How to run? 
# Rscript ./script/cal_loh_score_cp_adjusted.R CANCER_TYPE DIR_TO_TITAN
# Rscript ./script/cal_loh_score_cp_adjusted.R BRCA /fh/fast/ha_g/TCGA_analysis/TitanCNA

options(stringsAsFactors = FALSE, scipen=999)

library(data.table)
library(GenomicRanges)
library(dplyr)
library(crossval)
#library(GenomeInfoDb)
suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))

cal_LOH_score_without_segment_threshold <- function(pct, LOH_index, segment){
	#delete segments when LOH segment size is greater than thr*arm length
	del.id <- unlist(lapply(LOH_index, function(x){
		
		loh_seg_len <- segment[x, "Length"]
		arm_len <- segment[x,"arm_len"]
		#arm_len <- cytoband_arm_combined[ chr == segment[x,"Chromosome"] & arm == segment[x,"arm"],"len"]

		if(loh_seg_len >= pct*arm_len ) {
			x
		}
		#print(paste0(loh_seg_len, ",",pct*arm_len))
	}))
	#Record  deletion information (chrosomome arm info with segment legnth info)
	deleted_chr_info <- paste(segment[del.id, paste0(Chromosome, arm, "_", Length)], collapse=',')
	write.table(cbind(sample, deleted_chr_info), sprintf("./%s/LOH_%.2f_chr_del_info_%s.txt", type, pct, type), sep="\t", append=T, quote=F, col.names=F, row.names=F)

	if( !is.null(del.id) ){
		filt.seg <- seg[-del.id,]
	} else {
		filt.seg <- seg
	}

	#calculate LOH score after we delete events above
	total_length_sum <- sum(filt.seg[ ,"Length"])
	loh_length_sum <- sum(as.matrix(filt.seg[Corrected_MinorCN==0,"Length"])[,1])

	loh_score <- loh_length_sum / total_length_sum
	write.table(cbind(sample, loh_score), sprintf("./%s/LOH_pct_in_genome_%.2f_%s.txt",type, pct, type), sep="\t", append=T, quote=F, col.names=F, row.names=F)
	
	#number of LOH events
	num_events <- nrow(filt.seg[Corrected_MinorCN==0,])
	write.table(cbind(sample, num_events), sprintf("./%s/LOH_%.2f_num_events_%s.txt",type, pct, type), sep="\t", append=T, quote=F, col.names=F, row.names=F)

	loh_events_per_chr_arm <- filt.seg[Corrected_MinorCN==0, list(Freq =.N), by=list(Chromosome, arm)]
	loh_events_per_chr_arm_out <- paste(loh_events_per_chr_arm[, paste0(Chromosome, arm, "_", Freq)],collapse=',')
	write.table(cbind(sample, loh_events_per_chr_arm_out), sprintf("./%s/LOH_%.2f_events_per_arm_%s.txt",type, pct, type), sep="\t", append=T, quote=F, col.names=F, row.names=F)

	#Record arm segment fraction for LOH
	#filt.seg[Corrected_MinorCN==0, arm_seg_frac := round(Length/arm_len, digit=3)]

	loh_seg_sum_per_arm <- filt.seg[Corrected_MinorCN==0, sum(Length), by=list(Chromosome, arm)]
	arm_len_per_arm <- filt.seg[Corrected_MinorCN==0, unique(arm_len), by=list(Chromosome, arm)]
	
	loh_seg_sum_per_arm[, arm_seg_frac := round(loh_seg_sum_per_arm$V1/arm_len_per_arm$V1, digit=3)]
	
	loh_seg_fraction_per_arm <-	paste(loh_seg_sum_per_arm[, paste0(Chromosome, arm, "_", arm_seg_frac)], collapse=',')
	write.table(cbind(sample, loh_seg_fraction_per_arm), sprintf("./%s/LOH_%.2f_seg_fraction_per_arm_%s.txt",type, pct, type), sep="\t", append=T, quote=F, col.names=F, row.names=F)

	#For each cancer type, what is the distribution of the average LOH segment size (length of LOH segment)
	loh_mean_seg_size <- mean(as.matrix(filt.seg[Corrected_MinorCN == 0, "Length"])[,1])
	loh_median_seg_size <- median(as.matrix(filt.seg[Corrected_MinorCN == 0, "Length"])[,1])
	write.table(cbind(sample, loh_mean_seg_size), sprintf("./%s/LOH_%.2f_mean_seg_size_%s.txt",type, pct, type), sep="\t", append=T, quote=F,col.names=F, row.names=F)
	write.table(cbind(sample, loh_median_seg_size), sprintf("./%s/LOH_%.2f_median_seg_size_%s.txt",type, pct, type), sep="\t", append=T, quote=F,col.names=F, row.names=F)		

}

cal_LOH_score_with_segment_threshold <- function(pct, LOH_index, segment, seg_len_thr){
	#delete segments when LOH segment size is greater than thr*arm length (to remove chromosomal arm events) but keep segment size greater than seg_len_thr (10MB,15MB)
	del.id <- unlist(lapply(LOH_index, function(x){

		loh_seg_len <- segment[x,"Length"]
		arm_len <- segment[x,"arm_len"]
#		arm_len <- cytoband_arm_combined[chr==segment[x,"Chromosome"] & arm==segment[x,"arm"],"len"]
		
		#HRD-LOH score was defined as the number of LOH regions longer than 15 Mb but shorter than the whole chromosome (Abkevich et al, 2012)
		#So delete segmenets when it's greater than certain percentage of arm length + also delete when it's shorter thatn 10/15MB threshold.
		if(loh_seg_len >= pct*arm_len | loh_seg_len <= seg_len_thr) {
			x
		}
	}))
	#print( paste(segment[del.id, paste0(Chromosome, arm, Length)],collapse=','))
	
	deleted_chr_info <- paste(segment[del.id, paste0(Chromosome, arm, "_", Length)], collapse=',')
	write.table(cbind(sample, deleted_chr_info), sprintf("./%s/LOH_%.2f_with_%d_chr_del_info_%s.txt",type, pct, seg_len_thr, type), sep="\t", append=T, quote=F, col.names=F, row.names=F)

	if(!is.null(del.id)){
		filt.seg <- seg[-del.id,]
	} else {
		filt.seg <- seg
	}

	total_length_sum <- sum(filt.seg[ ,"Length"])
	loh_length_sum <- sum(as.matrix(filt.seg[Corrected_MinorCN==0, "Length"])[,1])
	loh_score <- loh_length_sum / total_length_sum
	write.table(cbind(sample, loh_score), sprintf("./%s/LOH_pct_in_genome_%.2f_with_%d_%s.txt",type, pct, seg_len_thr, type), sep="\t", append=T, quote=F, col.names=F, row.names=F)
	
	#number of LOH events
	num_events <- nrow(filt.seg[Corrected_MinorCN==0,])
	write.table( cbind(sample, num_events), sprintf("./%s/LOH_%.2f_with_%d_num_events_%s.txt",type, pct, seg_len_thr, type), sep="\t", append=T, quote=F, col.names=F, row.names=F)

	loh_events_per_chr_arm <- filt.seg[Corrected_MinorCN==0, list(Freq =.N), by=list(Chromosome,arm)]
	loh_events_per_chr_arm_out <- paste(loh_events_per_chr_arm[ ,paste0(Chromosome, arm, "_", Freq)],collapse=',')
	write.table( cbind(sample, loh_events_per_chr_arm_out), sprintf("./%s/LOH_%.2f_with_%d_events_per_arm_%s.txt",type, pct, seg_len_thr, type), sep="\t", append=T, quote=F, col.names=F, row.names=F)


	#Record arm segment fraction for all LOH events
	# filt.seg[Corrected_MinorCN==0, arm_seg_frac := round(Length/arm_len, digit=3)]
	# loh_seg_fraction_per_arm <-	paste(filt.seg[Corrected_MinorCN==0, paste0(Chromosome, arm, "_", arm_seg_frac)], collapse=',')
	# write.table(cbind(sample, loh_seg_fraction_per_arm), sprintf("LOH_%.2f_with_%d_seg_fraction_per_arm_%s.txt",pct, seg_len_thr, type), sep="\t", append=T, quote=F, col.names=F, row.names=F)

	loh_seg_sum_per_arm <- filt.seg[Corrected_MinorCN==0, sum(Length), by=list(Chromosome, arm)]
	arm_len_per_arm <- filt.seg[Corrected_MinorCN==0, unique(arm_len), by=list(Chromosome, arm)]
	
	loh_seg_sum_per_arm[, arm_seg_frac := round(loh_seg_sum_per_arm$V1/arm_len_per_arm$V1, digit=3)]
	
	loh_seg_fraction_per_arm <-	paste(loh_seg_sum_per_arm[, paste0(Chromosome, arm, "_", arm_seg_frac)], collapse=',')
	write.table(cbind(sample, loh_seg_fraction_per_arm), sprintf("./%s/LOH_%.2f_with_%d_seg_fraction_per_arm_%s.txt",type, pct, seg_len_thr, type), sep="\t", append=T, quote=F, col.names=F, row.names=F)


	#For each cancer type, what is the distribution of the average LOH segment size (length of LOH segment)
	loh_mean_seg_size <- mean(as.matrix(filt.seg[ Corrected_MinorCN==0, "Length" ])[,1])
	loh_median_seg_size <- median(as.matrix(filt.seg[ Corrected_MinorCN==0, "Length" ])[,1])
	write.table( cbind(sample, loh_mean_seg_size), sprintf("./%s/LOH_%.2f_with_%d_mean_seg_size_%s.txt",type, pct, seg_len_thr, type), sep="\t", append=T, quote=F,col.names=F, row.names=F)
	write.table( cbind(sample, loh_median_seg_size), sprintf("./%s/LOH_%.2f_with_%d_median_seg_size_%s.txt",type, pct, seg_len_thr, type), sep="\t", append=T, quote=F,col.names=F, row.names=F)		

}

args <- commandArgs(trailingOnly = TRUE)

type <- args[1] #Truseq_ICE
dir_to_titan <- args[2] #/fh/fast/ha_g/projects/SU2C/wes/Analysis/TitanCNA_Truseq_ICE/

print(type)

#output will be stored at this dir
dir.create(type, showWarnings = FALSE)

genome_hg38 <- BSgenome.Hsapiens.UCSC.hg38
genome_hg38_info <- seqinfo(genome_hg38)[seqlevels(genome_hg38)[1:22]]

file_cytoband <- "http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/cytoBand.txt.gz"

cytoband <- fread(file_cytoband,col.names=c("chr", "start", "end", "band", "stain"))

centromere <- subset(cytoband, stain=="acen" & chr!="chrY" & chr!="chrX", select=-stain)
centromere.gr <- GRanges(ranges=IRanges(start=centromere$start, 
							end=centromere$end),
		               	   	strand=rep('*', nrow(centromere)),
		                	seqnames=centromere$chr)

centromere.gr <- sortSeqlevels(centromere.gr)
seqinfo(centromere.gr) <- genome_hg38_info

cytoband_without_centromere <- subset(cytoband, stain!="acen", select=-stain)
cytoband_arm_minStart <- cytoband_without_centromere[ , .(start = min(start)), by = .(chr, arm = substring(band, 1, 1)) ]
cytoband_arm_maxEnd <- cytoband_without_centromere[ , .(end = max(end)), by = .(chr, arm = substring(band, 1, 1)) ]

cytoband_arm_combined <- cbind(cytoband_arm_minStart, cytoband_arm_maxEnd$end)
colnames(cytoband_arm_combined)[4] <- "end" 
cytoband_arm_combined[ ,arm_len := end-start]

cytoband_arm.gr <- GRanges(seqnames = cytoband_arm_combined$chr, 
								ranges = IRanges(start = cytoband_arm_combined$start, 
								end = cytoband_arm_combined$end), 
								arm = cytoband_arm_combined$arm,
								arm_len = cytoband_arm_combined$arm_len)
# cytoband_arm.gr = sortSeqlevels(cytoband_arm.gr)
# cytoband_arm.gr = keepSeqlevels(cytoband_arm.gr, seqlevels(genome_hg38_info),pruning.mode="coarse")

#############
#load seg file 
#############
##PATH TO DIR (TitanCNA run)
Opt_Sol_Dir <- sprintf("%s/titan/hmm/optimalClusterSolution", dir_to_titan)

segFiles <- list.files(Opt_Sol_Dir, pattern="titan.ichor.seg.txt", full.name = TRUE)
# cat(segFiles, "\n")
names(segFiles) <- unlist(lapply(strsplit(basename(segFiles),"_"),'[[',1))
# cat(names(segFiles), "\n")

l <- lapply(segFiles, fread)
segFiles.dt <- rbindlist(l) 

#only consider chr1-22
titan_seg <- segFiles.dt[Chromosome != "chrX",]
samples <- unique(titan_seg$Sample)

#cellular fraction threshold
cp_threshold <- 0.8


for (i in 1:length(samples)){
#for (i in 1:10){	
	print(i)

	sample <- samples[i]
	#print(sample)
	seg <- titan_seg[Sample == sample,]
	seg[ ,Length := End-Start]
	
	# stick to clonal copy number events only
	# Remove segments where cellular prevalence less than 0.8 (CP < 0.8)
	# However keep cf=NA (NEUT) and also Keep when cn >=4 
	subclonal_id <- which(seg$Cellular_Prevalence < cp_threshold & seg$Corrected_Copy_Number <4)
	#print(length(subclonal_id))

	if( length(subclonal_id) != 0 ){
		seg <- seg[-subclonal_id,]
	}
	

	seg.gr <- makeGRangesFromDataFrame(seg, keep.extra.columns=TRUE)

	#remove centromere -- centromere has already removed from titan
	#print(table(seg.gr %over% centromere.gr))
	seg.gr <- seg.gr[!(seg.gr %over% centromere.gr)]

	##annotate p and q arm info from cytoband
	#!!!!!!!!!seg where spans p and q arm will be ignored for now!!!!!!!!!
	query_id <- queryHits(findOverlaps(seg.gr, cytoband_arm.gr, type="within"))
	subject_id <- subjectHits(findOverlaps(seg.gr, cytoband_arm.gr, type="within"))

	#annotate seg with cytoband p,q arm and arm length information
	seg[query_id,"arm"] <- cytoband_arm_combined[subject_id,"arm"]
	seg[query_id,"arm_len"] <- cytoband_arm_combined[subject_id,"arm_len"]
	
	#restrict to seg to have unique p,q arm (don't allow segment have both p/q arm)
	#seg = seg[!is.na(arm),]

	# Exclude LOH segments that are greater than threshold of the chr arm
	#Corrected_MinorCN = 0 is LOH
	LOH_idx <- which(seg$Corrected_MinorCN == 0)
	#for(thr in c(0.25, 0.5, 0.75, 0.90)){

	cal_LOH_score_without_segment_threshold(0.75, LOH_idx, seg)
	cal_LOH_score_with_segment_threshold(0.75, LOH_idx, seg, 10000000) #10MB
	cal_LOH_score_with_segment_threshold(0.75, LOH_idx, seg, 15000000) #15MB
	cal_LOH_score_without_segment_threshold(0.9, LOH_idx, seg)
	cal_LOH_score_with_segment_threshold(0.9, LOH_idx, seg, 10000000) #10MB
	cal_LOH_score_with_segment_threshold(0.9, LOH_idx, seg, 15000000) #15MB
}



	