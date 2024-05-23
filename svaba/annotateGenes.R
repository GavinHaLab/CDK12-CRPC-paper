# 20231123 Armand Bankhead
# annotate breakpoints with gene information

options(stringsAsFactors=F)

library(optparse)
library(GenomicRanges)
library(dplyr)

## TPERSSE CHANGES:
## changed input and output to be dirs, not files, to allow for multiple files to be processed at once


# ** deal with arguments **
option_list <- list(
	make_option(c("--inputSVdir"), type = "character", help = "output dir containing output files from combineSVABAandTitan.R script"),
	make_option(c("--outputSVdir"), type = "character", help = "name of updated dir where updated files with gene annotations will be written")
)
parseobj <- OptionParser(option_list=option_list, usage = "usage: Rscript %prog [options]")
opt <- parse_args(parseobj)

##

# input dir is full of dirs, with each dirname being the sample name
# looks like this: {inputSVdir}/sampleName/sampleName.svabaTitan.sv.txt
# output dir will just be the files, no nested dirs
# looks like this: {outputSVdir}/sampleName.svabaTitan.genes.sv.txt

# get list of sample names
sampleNames = list.files(opt$inputSVdir, full.names = FALSE, recursive = FALSE) # change this to list.files if you want to include files in the top level dir
# cat("sampleNames: ", sampleNames, "\n")
# cat("inputSVdir: ", opt$inputSVdir, "\n")

# for each sample, read in the file, annotate the genes, and write out the new file
for (sampleName in sampleNames) {
    inFile1 = paste0(opt$inputSVdir, sampleName)
    # outFile1 = paste0(opt$outputSVdir, "/", sampleName, ".genes.sv.txt")
    # need to change so that it creates output name based on input name, but split on extension, and with genes added before .txt
    outFile1 = paste0(opt$outputSVdir, "/", sub("\\.txt$", ".genes.txt", sampleName))
    message(sampleName)

    # inFile1 <- opt$inputSVFile
    # #inFile1 = "intermediate/04/PCFHX_2001B_P29_CL.svabaTitan.sv.txt"
    # outFile1 <- opt$outputSVFile
    #outFile1 = "intermediate/04/PCFHX_2001B_P29_CL.svabaTitan.genes.sv.txt"
    gene_annotation = "/fh/fast/ha_g/WCDT_CRPC_WGS/Analysis/matrices_CNA/gene/GRCh38.p12.ensembl.gene.annotations.sorted.include.AREnhancer.txt"
    # **

    data1 = read.delim(inFile1)

    # ** read gene annotation **
    genomeStyle = "UCSC"
    chrs = c(1:22,"X")
    seqlevelsStyle(chrs) = genomeStyle
    message('Reading gene annotation')
    genes = read.delim(gene_annotation)
    genes$Chromosome=paste0('chr',genes$Chr)
    genes$Start = genes$cdsStart
    genes$End = genes$cdsEnd
    genes$length= genes$End - genes$Start

    # keep longest transript
    genes = genes[order(genes$length,decreasing=T),]
    genes = genes[!duplicated(genes$Ensg_ID) & genes$Chromosome %in% chrs,]
    genes = select(genes,ensembl_gene = Ensg_ID,Chromosome,Start,End,gene=HGNC_symbol)

    # ** sourced from /fh/fast/ha_g/user/gha/software/git/TitanCNA/R/utils.R **
    getOverlap <- function(svs, genes) {
    gene_list = rep(NA,nrow(svs))	   
    message('Finding overlaps')
    x <- as(svs, "GRanges")
    y <- as(genes, "GRanges")
    hits <- findOverlaps(query = x, subject = y, type = 'any')

    # polish a data frame
    hits = as.data.frame(hits)
    hits$sv_idx = hits$queryHits

    # get genes and remove the blanks
    hits$gene = genes$gene[hits$subjectHits]
    hits = hits[hits$gene != '',]

    # get rid of dups by ordering
    hits = hits[order(hits$gene),]
    hits = hits[!duplicated(hits$queryHits),]
    message('Found ', nrow(hits), ' overlaps')
    gene_list[hits$queryHits] = hits$gene
    return(gene_list)
    }
    # **

    # create overlap queries for first and second partner
    ## do some quick cleaning of chromosome_1 and chromosome_2 to make sure there are no NAs
    ## also, just check start_1 and start_2 to be safe
    data1 = data1[!is.na(data1$chromosome_1) & !is.na(data1$chromosome_2) & !is.na(data1$start_1) & !is.na(data1$start_2),]

    svs1 = select(data1,id=SV.id,Chromosome = chromosome_1,Start = start_1)
    svs1$End=svs1$Start+1
    svs2 = select(data1,id=SV.id,Chromosome = chromosome_2,Start = start_2)
    svs2$End=svs2$Start+1

    # get overlapping genes
    data1$gene1 = getOverlap(svs1,genes)
    data1$gene2 = getOverlap(svs2,genes)

    print(paste0('writing ', outFile1, ' ...'))
    write.table(data1,outFile1,quote=F,row.names=F,sep="\t")
}
