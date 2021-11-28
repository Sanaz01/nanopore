#!/bin/bash

###Extracting 5mC DNA methylation frequencies at CpG sites using Nanopolish software 

##Created by Sanaz Agarwal
##9 September 2021
##Email: sanaz.iitd@gmail.com

##Software and version used
#minimap2 - 2.22-r1105-dirty
#samtools - 1.9
#nanopolish - 0.13.3

echo $PWD
echo "Making a new directory nanopolish/"
mkdir nanopolish && cd nanopolish

wd=$PWD
basecall_dir=$1      #directory containing basecalling results
fast5_dir=$2         #directory containing fast5 files generated from nanopore sequencing
ref_gen=$3           #reference genome to be used for alignment against fastq reads

#Input files: Fast5 files: ./data/fast5_files/
#             Reference genome: ./data/reference_genome.fasta
#             Guppy basecalled fastq files: 

#Output directory - /nfs_master/guest2/sanaz/analysis/basecalling_output
#Copy all fastq files (after decompressing) from pass/ and fail/ into /nfs_master/guest2/sanaz/analysis/fastq_basecalled


cd ${basecall_dir}/pass/                       #change directory to where passed fastq files are
cat *.fastq > combined_basecalled.fastq		   #combine all fastq files into a single file					
mv combined_basecalled.fastq ${wd}             #move combined file to original directory

###########Data preprocessing
nanopolish index -d ${fast5_dir} ${wd}/combined_basecalled.fastq

##Output files - combined_basecalled.fastq.index
#                combined_basecalled.fastq.index.gzi
#                combined_basecalled.fastq.index.readdb
#                combined_basecalled.fastq.index.fai 

###########Aligning reads to the reference genome
minimap2 -a -x map-ont ${ref_gen} ${wd}/combined_basecalled.fastq > minimap2_output.sam   #Output in SAM format, index file created - reference_genome.fasta.fai
samtools view -S -b minimap2_output.sam > minimap2_output.bam                                      #SAM to BAM conversion
samtools sort minimap2_output.bam -o minimap2_output.sorted.bam                                    #Sort BAM file 
samtools index minimap2_output.sorted.bam

###########Calling methylation
nanopolish call-methylation -t 12 -r ${wd}/combined_basecalled.fastq -b minimap2_output.sorted.bam -g ${ref_gen} --progress -q cpg > methylation_calls.tsv

#./calculate_methylation_frequency.py methylation_calls.tsv > methylation_frequency.tsv              #helper script calculate methylation frequencies; present in nanopolish/scripts/
#./calculate_methylation_frequency.py -s methylation_calls.tsv > methylation_frequency_split.tsv

##remove extrafiles generated
#rm combined_basecalled* minimap2*