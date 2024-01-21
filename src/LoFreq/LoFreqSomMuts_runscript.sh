#!/bin/bash

################ This script calls somatic mutations using LoFreq package in tumor-normal pair mode on whole-genome sequencing (WGS) samples
################ LoFreq v2.1.5 is in usage, 
################ NOTE: you need to install LoFreq package in Linux environment. Also all scripts were written in R, so have R and installed optparse package.

# srun -N 1 --cpus-per-task 16 --mem-per-cpu=16gb -t 5:00:00 -p interactive --pty bash

################ R-version 4.0.4
module load R/4.0.4

################ python 3.6.3
module load python/3.6.3

################ samtools (v1.17)
module load samtools/1.17

################ htslib
module load htslib/1.9

################ bcftools (v1.6)
module load bcftools/1.6

export script_dir=/scratch.global/balay011/WGS_skChordomas_Bai/SomaticSNV_calling/scripts
export DBSNP_FILE=/scratch.global/balay011/WGS_skChordomas_Bai/reference/All_20180418.vcf.gz
export SAMTOOLS_PATH=/common/software/install/manual/samtools/1.17_gcc-7.2.0/bin/samtools
export REF_FILE=/scratch.global/balay011/WGS_skChordomas_Bai/reference/Homo_sapiens_assembly38.fasta
export ANNO_FILE=/scratch.global/balay011/WGS_skChordomas_Bai/scripts/WGS_chordomaBai_anno.txt
export BED_FILE=/scratch.global/balay011/WGS_skChordomas_Bai/reference/Homo_sapiens_assembly38_chrs.bed
export LOFREQ_PATH=/home/aventeic/balay011/lofreq/bin/lofreq
export WORK_DIR=/scratch.global/balay011/WGS_skChordomas_Bai/SomaticSNV_calling/work/lofreq
export CALLER_FILE=/scratch.global/balay011/WGS_skChordomas_Bai/SomaticSNV_calling/scripts/LoFreqSomaticCaller.r
export MIN_COV=7
export NTHREADS=16

Rscript $script_dir/runLoFreqSomaticMutationCalling.r \
--work-dir=$WORK_DIR \
--caller-file=$CALLER_FILE \
--anno-file=$ANNO_FILE \
--ref-file=$REF_FILE \
--bed-file=$BED_FILE \
--cpus=$NTHREADS \
--dbsnp=$DBSNP_FILE \
--lofreq-path=$LOFREQ_PATH \
--min-cov=$MIN_COV \
--samtools-path=$SAMTOOLS_PATH

################ End of LoFreq run

################ Extract somatic mutations from LoFreq run
export ANNO_FILE=/home/aventeic/balay011/scripts/WGS_chordomaBai_sommutsfilter.txt
export OUTPUT_DIR=/home/aventeic/balay011/WGS_skChordomas_Bai/SomaticSNV_calling/lofreq

Rscript $script_dir/Call_sommutsfiltering.r \
--script-dir  $script_dir \
--anno-file $ANNO_FILE \
--mode lofreq \
--output-dir $OUTPUT_DIR \
--output-vcf 0
  
