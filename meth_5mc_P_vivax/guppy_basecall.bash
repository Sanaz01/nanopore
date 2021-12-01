#!/bin/bash

##Created by Sanaz Agarwal
##25 August 2021
##Email: sanaz.iitd@gmail.com


###Performing modified basecalling using Guppy basecalling software developed by Oxford Nanopore Technologies
###Used to detect 5mC modifed CpG sites in Plasmodium vivax sal 1

##Version: Guppy - 5.0.14 (GPU based)

config_file=$1        ##onfiguration file used. See list of configuration files available as 'ls /ont-guppy/data/*.cfg'
fast5_dir=$2          ##directory containing fast5 files generated from Nanopore Sequencing, fast5 files will be taken recursively within this directory
save_dir=$3           ##output directory for saving fastq files   
ref_gen=$4            ##reference genome FASTA file for alignment  

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
guppy_basecaller -c ${config_file} -i ${fast5_dir} -s ${save_dir} --recursive ----align_ref ${ref_gen} --bam_out\
				 --num_callers 8 --gpu_runners_per_device 6 --chunks_per_runner 512 --chunk_size 3000 --chunks_per_caller 10000 --device cuda:0:50%      #enhance computation power                                                        

#Optional parameters: --compress_fastq: .gz compress the fastq files generated 
#					  --fast5_out     : allow generation of fast5 files with additional data stored
#					  --align_ref     : provide refernce genome for internal aligment with generated fastq reads
#					  --bam_out       : output bam alignment files in case aligment option is enabled


#  from FYd1

#data/  from FYd1 (using Rerio config file)

#data2/   from EBd_ELd2
#nohup guppy_basecaller -c dna_r9.4.1_450bps_modbases_5mc_hac.cfg -i /nfs_master/guest2/sanaz/bits/data2  -s /nfs_master/guest2/sanaz/bits/basecalled_data2 --align_ref reference_genome.fasta --recursive --num_callers 10 --gpu_runners_per_device 4 --chunks_per_runner 512 --chunk_size 2000 --chunks_per_caller 10000 --device cuda:0:95% --bam_out > nohup_basecalling_data2 &

#barcoding
#guppy_barcoder --input_path /nfs_master/guest2/sanaz/bits/EBd/basecalled_data2/pass/fastq  --save_path /nfs_master/guest2/sanaz/bits/EBd/basecalled_data2/barcoding --config configuration.cfg --barcode_kits EXP-NBD104


#data2/   from EBd_ELd2 (Rerio config file)
#nohup guppy_basecaller -c res_dna_r941_min_modbases_5mC_5hmC_v001.cfg -i /nfs_master/guest2/sanaz/bits/EBd/data2 -s /nfs_master/guest2/sanaz/bits/EBd/basecalled_rerio --align_ref reference_genome.fasta --recursive --num_callers 10 --gpu_runners_per_device 4 --chunks_per_runner 512 --chunk_size 2000 --chunks_per_caller 10000 --device cuda:0:95% --bam_out > nohup_basecalling_rerio &

##sort - already sorted bam files
## combine bam files
#samtools merge merged.bam *.bam     ( may not be working)
#[W::bam_merge_core2] No @HD tag found.
#samtools cat *_0.bam -o concat.bam

#modbam2bed -e -m 5mC --cpg -t 4 -r chr20 $reference_genome $data_folder/*.bam > $bedmethyl
#modbam2bed -e -m 5mC --cpg -t 12 reference_genome.fasta merged.bam > bedmethyl.cpg

