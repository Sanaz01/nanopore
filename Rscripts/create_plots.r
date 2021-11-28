#!/usr/bin/env Rscript

##Find average methylation frequency in genomic regions and write them in a table
##Downstream processing after Guppy basecalling and modified base identification

##Created by Sanaz Agarwal
##14 October 2021
##Email: sanaz.iitd@gmail.com
##For use in Rstudio environment

#install Genomic Ranges, IRanges
lib = list("IRanges", "GenomicRanges", "S4Vectors", "dplyr")
for (i in lib) { if (!require(i)) install.packages(i); library(i) }

#Input files: bedfile  - BED file from methylation calling 
#             reference.gtf - GTF file of reference genome

##If operating from terminal
#args <- commandArgs(trailingOnly = TRUE)
#gtf = read.delim("args[1]", header=F)
#bed = read.delim("args[2]", header=F) 

##set a working directory, containing reference.gtf and bedfile 
#setwd('working/directory/')

bd = read.delim("bedfile", header=FALSE)
gtf = read.delim("reference.gtf", header=FALSE)

##remove CpG sites with zero methylation frequency
bed = bd[which(bd$V11!="0"),]

##Extract genomic regions from GTF file
gene=gtf[which(gtf$V3=="gene"),]
intron=gtf[which(gtf$V3=="intron"),]
exon=gtf[which(gtf$V3=="exon"),]


name = list("bed", "exon", "intron", "gene")
for (i in name){

find_gr <- function(arg){

#make GRanges
gr_bed = GRanges(seqnames = bed$V1, range = IRanges(start=bed$V2, end=bed$V3))
gr_exon = GRanges(seqnames = exon$V1, range = IRanges(start=exon$V4, end=exon$V5))
gr_intron = GRanges(seqnames = intron$V1, range = IRanges(start=intron$V4, end=intron$V5))
gr_gene = GRanges(seqnames = gene$V1, range = IRanges(start=gene$V4, end=gene$V5))
	
#finding methlation frequency per Cpg site in each genomic region and find average over that genomic region.	
ov_meth2exon = findOverlaps(query = gr_bed, subject = gr_exon)
ov_m2e = as.data.frame(ov_meth2exon)
split_ovm2e = split(ov_m2e, ov_m2e$subjectHits)
meth_exonwise = sapply(split_ovm2e, function(x) mean(bed$V11[x$queryHits]))

ov_meth2intron = findOverlaps(query = gr_bed, subject = gr_intron)
ov_m2i = as.data.frame(ov_meth2intron)
split_ovm2i = split(ov_m2i, ov_m2i$subjectHits)
meth_intronwise = sapply(split_ovm2i, function(x) mean(bed$V11[x$queryHits]))

ov_meth2gene = findOverlaps(query = gr_bed, subject = gr_gene)
ov_m2g = as.data.frame(ov_meth2gene)
split_ovm2g = split(ov_m2g, ov_m2g$subjectHits)
meth_genewise = sapply(split_ovm2g, function(x) mean(bed$V11[x$queryHits]))

###Making boxplots
#regions = list("meth_exonwise", meth_intronwise")
#boxplot(regions, names=c("exon", "intron"), ylab="Average methylation of CpG sites", outline=F)
##for manipulation use -  pars=list(par(mar=c(8,5,3,2))), cex.lab=0.75

###Making violin plots
#df_exon = data.frame(value=meth_exonwise, variable="exon")
#df_intron = data.frame(value=meth_intronwise, variable="intron")
#ggplot(rbind(df_exon, df_intron), aes(x=variable, y=value)) + geom_violin(scale="width", adjust = 1, width=0.5) + geom_boxplot(width=0.03)

###Making histograms
#hist(v2, xlab = "Methylation frequencies", col = "green", border = "black", ylim=c(0,50000))

##making cumulative histogram
#h = hist(v2, xlab = "Methylation frequencies", col = "green", border = "black", ylim=c(0,50000))
#h$counts = cumsum(h$counts)
#plot(h, xlab ="Methylation frequencies", ylab ="No of CpG sites", main = "Cumulative histogram of non-zero methylation frequencies")

	

###For finding intergenic regions
##gr_gene
##gr_len  - name, start and end of reference genome fasta sequences length)
#len = read.delim("len_reference_fasta", header=F)
#gr_len = GRanges(seqnames = len$V1, range = IRanges(start=len$V2, end=len$V3))
#gr_intergenic = setdiff(gr_len, gr_gene)

###for finding promoter regions
##gr_gene
#pos = gene[gene$V7=="+",]
#neg = gene[gene$V7=="-",]
#pos$V5 = pos$V4
#pos$V4 = pos$V5-1000
#neg$V4 = neg$V5
#neg$V5 = neg$V4+1000
#promoters = rbind(pos, neg)
#gr_promoters = GRanges(seqnames = promoters$V1, range = IRanges(start=promoters$V4, end=promoters$V5))
