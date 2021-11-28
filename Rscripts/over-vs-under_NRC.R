#!/usr/bin/env Rscript

##Find average methylation frequency in genomic regions and write them in a table
##Downstream processing after Guppy basecalling and modified base identification

##Created by Sanaz Agarwal
##14 October 2021
##Email: sanaz.iitd@gmail.com
##For use in Rstudio environment

#Input files: bedfile  - BED file from methylation calling 
#             nrc_dep - scaffold name, start and end location of under represented regions in NRC
#             nrc_enr - scaffold name, start and end location of over represented regions in NRC

##If operating from terminal
#args <- commandArgs(trailingOnly = TRUE)
#bed = read.delim("args[1]", header=F)
#dep = read.delim("args[2]", header=F)

##set a working directory, containing reference.gtf and bedfile 
#setwd('working/directory/')

#install necessary libraries
lib = list("IRanges", "GenomicRanges", "ggplot2", "ggstatsplot")
for (i in lib) { if (!require(i)) install.packages(i); library(i) }

bd = read.delim("bedfile", header=FALSE)            #load BED file in Rstudio
bed = bd[which(bd$V11!="0"),]                       #remove CpG sites with zero methylation frequency
dep = read.delim("nrc_dep", header=F)               #NRC under represented regions
enr = read.delim("nrc_enr", header=F)               #NRC over represented regions

##making GRanges
gr_bed = GRanges(seqnames = bed$V1, range = IRanges(start=bed$V2, end=bed$V3))
gr_dep = GRanges(seqnames = dep$V1, range = IRanges(start=dep$V2, end=dep$V3))
gr_enr = GRanges(seqnames = enr$V1, range = IRanges(start=enr$V2, end=enr$V3))

##finding methlation frequency per Cpg site in each genomic region and find average over that genomic region.
ov_meth2dep = findOverlaps(query = gr_bed, subject = gr_dep)
ov_dep = as.data.frame(ov_meth2dep)
split_dep = split(ov_dep, ov_dep$subjectHits)
meth_dep = sapply(split_dep, function(x) mean(bed$V11[x$queryHits]))

ov_meth2enr = findOverlaps(query = gr_bed, subject = gr_enr)
ov_enr = as.data.frame(ov_meth2enr)
split_enr = split(ov_enr, ov_enr$subjectHits)
meth_enr = sapply(split_enr, function(x) mean(bed$V11[x$queryHits]))

##preparation of data for plotting
df_enr= data.frame(value=meth_enr, variable="Over rep.")
df_dep= data.frame(value=meth_dep, variable="Under rep.")
df = rbind(df_enr,df_dep)
colnames(df) = c("Frequency", "NRC_regions")

##Boxplot of methylation in over and under rep. NRC regions
ggplot(df) + aes(x=NRC_regions, y=Frequency) + geom_boxplot(fill = "#FFDB6D", color = "#C4961A") + theme_minimal()

##Test for normal distribution
hist(subset(df, NRC_regions=="Over rep.")$Frequency, main="Average methylaion in Over Rep. NRC", xlab="Methylation frequency", ylab="No. of regions")
hist(subset(df, NRC_regions=="Under rep.")$Frequency, main="Average methylaion in Under Rep. NRC", xlab="Methylation frequency", ylab="No. of regions")

shapiro.test(subset(df, NRC_regions == "Over rep.")$Frequency)
shapiro.test(subset(df, NRC_regions == "Under rep.")$Frequency)

##Wilcox t-test
wilcox.test(df$Frequency ~ df$NRC_regions)

##Plot comparision of over and under NRC regions
ggbetweenstats( # independent samples
  data = df,
  x = NRC_regions,
  y = Frequency,
  plot.type = "box", # for boxplot
  type = "nonparametric", # for wilcoxon
  centrality.plotting = FALSE # remove median
)

##Plot distribution histograms of over and under rep. NRCs
gghistostats(data=df_enr, x=value, xlab="Methylation frequency", title="Average methylaion in Over Rep. NRC")
gghistostats(data=df_dep, x=value, xlab="Methylation frequency", title="Average methylaion in Under Rep. NRC")
