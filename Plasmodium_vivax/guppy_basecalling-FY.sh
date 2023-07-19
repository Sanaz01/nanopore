#!/bin/bash

###Performing modified basecalling using Guppy basecalling software developed by Oxford Nanopore Technologies

##Created by Sanaz Agarwal
##20 February 2022
##Email: sanaz.iitd@gmail.com

##Version: Guppy - 6.0.0 (GPU based)


config_file=$1        ##onfiguration file used. See list of configuration files available as 'ls /ont-guppy/data/*.cfg' (We choose "dna_r9.4.1_450bps_modbases_5mc_hac.cfg")
fast5_dir=$2          ##directory containing fast5 files generated from Nanopore Sequencing, fast5 files will be taken recursively within this directory
save_dir=$3           ##output directory for saving fastq files   
ref_gen=$4            ##reference genome FASTA file for alignment  
bedfile=$5              ##Name of final BED file containing methylation frequency of each CG site identified in barcoding 1


##GPU enabled basecalling command. Using V100 Tesla GPU CUDA.
guppy_basecaller -c $config_file -i $fast5_dir -s $save_dir --align_ref $ref_gen --device cuda:0,1:60% --bam_out --trim_adapters > nohup

modbam2bed -e -m 5mC --cpg -t 8 $ref_gen $save_dir/pass/*.bam > $bedfile

