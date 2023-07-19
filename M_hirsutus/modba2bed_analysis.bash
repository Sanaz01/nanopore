#!/bin/bash

###Methylation feature extraction using modbam2bed program (Feature of EPI2ME, by Oxford Nanopore Technologies)
##NOTE: Works only if --ref_align and --bam_out options are enabled while Guppy basecalling

##Created by Sanaz Agarwal
##15 September 2021
##Email: sanaz.iitd@gmail.com

##Software and version used
#modbam2bed - 0.3.1                    "https://github.com/epi2me-labs/modbam2bed.git"
#samtools - 1.13

basecall_dir=$1     #directory where output files of guppy basecalling is present
ref_gen=$2          #reference genome to be used for alignment against fastq reads
bedfile=$3          #outfile in BED format containing identified CpG sites with their methylation status (frequency)

#echo $PWD
#echo "Making a new directory mod_base/"
#mkdir mod_base && cd mod_base

##sort - already sorted bam files
## combine bam files
samtools merge merged.bam ${basecall_dir}/pass/*.bam      #merging BAM alignemnt files of only passed fastq reads.

#If samtools merge does not work, use 'samtools cat *_0.bam -o concat.bam'

#modbam2bed -e -m 5mC --cpg -t 4 -r chr20 ${ref_gen} $data_folder/*.bam > $bedmethyl
modbam2bed -e -m 5mC --cpg -t 12 ${ref_gen} merged.bam > ${bedfile}
#cd ../


