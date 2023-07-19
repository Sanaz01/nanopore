#!/bin/bash

###Wrapper function for modified basecalling, getting methylation frequencies, plotting methylation status

##Created by Sanaz Agarwal
##25 November 2021
##Email: sanaz.iitd@gmail.com

config_file=$1        ##onfiguration file used. See list of configuration files available as 'ls /ont-guppy/data/*.cfg'
fast5_dir=$2          ##directory containing fast5 files generated from Nanopore Sequencing, fast5 files will be taken recursively within this directory
save_dir=$3           ##output directory for saving fastq files   
ref_fasta=$4          ##reference genome FASTA file used for alignment
bedfile=$5            ##name of output BED file to be generated
ref_gtf=$6            ##reference genome GTF file used for analysis and plotting
nrc_enr=$7            ##tab delim file with scaffold name, start and end as 3 columns for NRC over represented regions
nrc_dep=$8            ##tab delim file with scaffold name, start and end as 3 columns for NRC over represented regions

##basecalling Guppy
sh ./guppy_basecalling.sh ${config_file} ${fast5_dir} ${save_dir}

##reading BAM files to get methylation frequencies
sh ./modbam2bed_analysis.sh ${save_dir} ${ref_genome} ${bedfile}

##Plot the differential methylation status in over and under represented NRC regions
Rscipt ../Rscipts/over-vs-under_NRC.R ${bedfile} ${ref_gtf} ${nrc_dep} ${nrc_enr}

#For further analysis
#python3 ../helper_scripts/fasta_length.py > len_reference_fasta
#Rscipt ../Rscipts/intergenic-distribution_NRC.R ${bedfile} ${ref_gtf} ${nrc_dep} ${nrc_enr} len_reference_fasta
