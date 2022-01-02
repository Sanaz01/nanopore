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
gene=gtf[which(gtf$V3=="gene"),]                        #extract genes from GTF file
intron=gtf[which(gtf$V3=="intron"),]                    #extract introns from GTF file
exon=gtf[which(gtf$V3=="exon"),]                        #extract exons from GTF file


#####name = list("bed", "exon", "intron", "gene")
#####for (i in name){
######find_gr <- function(arg){

#make GRanges
gr_bed = GRanges(seqnames = bed$V1, range = IRanges(start=bed$V2, end=bed$V3))
gr_exon = GRanges(seqnames = exon$V1, range = IRanges(start=exon$V4, end=exon$V5))
gr_intron = GRanges(seqnames = intron$V1, range = IRanges(start=intron$V4, end=intron$V5))
gr_gene = GRanges(seqnames = gene$V1, range = IRanges(start=gene$V4, end=gene$V5))
	
#finding methlation frequency per Cpg site in each genomic region and find average over that genomic region.	
ov_meth2exon = findOverlaps(query = gr_bed, subject = gr_exon)                                    #finding overlaps between GRanges
ov_m2e = as.data.frame(ov_meth2exon)                                                              #converting to data frame
split_ovm2e = split(ov_m2e, ov_m2e$subjectHits)                                                   #split meth freq values per genomic region 
meth_exonwise = sapply(split_ovm2e, function(x) mean(bed$V11[x$queryHits]))                       #find average of meth freq values associated with each genomic region

ov_meth2intron = findOverlaps(query = gr_bed, subject = gr_intron)
ov_m2i = as.data.frame(ov_meth2intron)
split_ovm2i = split(ov_m2i, ov_m2i$subjectHits)
meth_intronwise = sapply(split_ovm2i, function(x) mean(bed$V11[x$queryHits]))

ov_meth2gene = findOverlaps(query = gr_bed, subject = gr_gene)
ov_m2g = as.data.frame(ov_meth2gene)
split_ovm2g = split(ov_m2g, ov_m2g$subjectHits)
meth_genewise = sapply(split_ovm2g, function(x) mean(bed$V11[x$queryHits]))



#Merging Avg. meth. freq. per genomic region with the regions extracted from GTF file
e = data.frame(exon)
e$seq = seq_len(length(e$V1))
me = data.frame(meth_exonwise)
me$num = rownames(me)
me$num <- as.numeric(me$num)
se <- S4Vectors::merge(x=e,y=me,by.x="seq",by.y="num", all.x=TRUE, sort=F) %>% arrange(factor(seq, levels = e$seq))
df_exon  = subset(se, select=-c(seq))
write.table(df_exon, "Meth freq in exons", sep="\t", row.names = FALSE, col.names = FALSE, quote=FALSE)

i = data.frame(intron)
i$seq = seq_len(length(i$V1))
mi = data.frame(meth_intronwise)
mi$num = rownames(mi)
mi$num <- as.numeric(mi$num)
si <- S4Vectors::merge(x=i,y=mi,by.x="seq",by.y="num", all.x=TRUE, sort=F) %>% arrange(factor(seq, levels = i$seq))
df_intron  = subset(si, select=-c(seq))
write.table(df_intron, "Meth freq in introns", sep="\t", row.names = FALSE, col.names = FALSE, quote=FALSE)

g = data.frame(gene)
g$seq = seq_len(length(g$V1))
mg = data.frame(meth_genewise)
mg$num = rownames(mg)
mg$num <- as.numeric(mg$num)
sg <- S4Vectors::merge(x=g,y=mg,by.x="seq",by.y="num", all.x=TRUE, sort=F) %>% arrange(factor(seq, levels = g$seq))
df_gene  = subset(sg, select=-c(seq))
write.table(df_gene, "Meth freq in genes", sep="\t", row.names = FALSE, col.names = FALSE, quote=FALSE)

	

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
