#!/bin/bash

###Performing modified basecalling using Guppy basecalling software developed by Oxford Nanopore Technologies

##Created by Sanaz Agarwal
##24 February 2022
##Email: sanaz.iitd@gmail.com

##Version: Guppy - 6.0.0 (GPU based)

config_file=$1        ##onfiguration file used. See list of configuration files available as 'ls /ont-guppy/data/*.cfg' (We choose "dna_r9.4.1_450bps_modbases_5mc_hac.cfg")
fast5_dir=$2          ##directory containing fast5 files generated from Nanopore Sequencing, fast5 files will be taken recursively within this directory
save_dir=$3           ##output directory for saving fastq files   
ref_gen=$4            ##reference genome FASTA file for alignment  
bar_kit=$5            ##Name of barcoding kit used for sequencing/multiplexing. We used "EXP-NBD104"
bed1=$6               ##Name of final BED file containing methylation frequency of each CG site identified in barcoding 1
bed2=$7               ##Name of final BED file containing methylation frequency of each CG site identified in barcoding 2


#Optional parameters: --compress_fastq: .gz compress the fastq files generated 
#					  --fast5_out     : allow generation of fast5 files with additional data stored
#					  --align_ref     : provide refernce genome for internal aligment with generated fastq reads
#					  --bam_out       : output bam alignment files in case aligment option is enabled


                                           
##GPU enabled basecalling command. Using V100 Tesla GPU CUDA.
guppy_basecaller -c $config_file -i $fast5_dir -s $save_dir --align_ref $ref_gen --barcode_kits $bar_kit --device cuda:0,1:80% --bam_out --trim_adapters > nohup

For making BED files
modbam2bed -e -m 5mC --cpg -t 8 $ref_gen $save_dir/pass/barcode01/*.bam > $bed1
modbam2bed -e -m 5mC --cpg -t 8 $ref_gen $save_dir/pass/barcode02/*.bam > $bed2



#barcoding
#guppy_barcoder --input_path /nfs_master/guest2/sanaz/bits/EBd/basecalled_data2/pass/fastq  --save_path /nfs_master/guest2/sanaz/bits/EBd/basecalled_data2/barcoding --config configuration.cfg --barcode_kits EXP-NBD104

