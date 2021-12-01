#!/bin/bash

###Performing modified basecalling using Guppy basecalling software developed by Oxford Nanopore Technologies

##Created by Sanaz Agarwal
##25 August 2021
##Email: sanaz.iitd@gmail.com

##Version: Guppy - 5.0.14 (GPU based)

config_file=$1          ##onfiguration file used. See list of configuration files available as 'ls /ont-guppy/data/*.cfg'
fast5_dir=$2          ##directory containing fast5 files generated from Nanopore Sequencing, fast5 files will be taken recursively within this directory
save_dir=$3           ##output directory for saving fastq files   

##Inside the selected config file:
# Basecalling.
#model_file                          = template_r9.4.1_450bps_modbases_5mc_hac.jsn
#chunk_size                          = 2000
#gpu_runners_per_device              = 4
#chunks_per_runner                   = 512
#chunks_per_caller                   = 10000
#overlap                             = 50
#qscore_offset                       = -3.202
#qscore_scale                        = 1.572
#builtin_scripts                     = 1

##GPU enabled basecalling command. Using V100 Tesla GPU at CUDA 0.
guppy_basecaller -c ${config_file} -i ${fast5_dir} -s ${save_dir} --recursive\
				 --num_callers 8 --gpu_runners_per_device 6 --chunks_per_runner 512 --chunk_size 3000 --chunks_per_caller 10000 --device cuda:0:50%      #enhance computation power                                                        

#Optional parameters: --compress_fastq: .gz compress the fastq files generated 
#					  --fast5_out     : allow generation of fast5 files with additional data stored
#					  --align_ref     : provide refernce genome for internal aligment with generated fastq reads
#					  --bam_out       : output bam alignment files in case aligment option is enabled

