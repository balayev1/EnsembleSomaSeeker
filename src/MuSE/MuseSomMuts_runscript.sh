#!/bin/bash

################ This script calls somatic mutations using Muse package in tumor-normal pair mode on whole-genome sequencing (WGS) samples
################ Muse v2.0 is in usage, 
################ NOTE: you need to install Muse package in Linux environment. Also all scripts were written in R, so have R and installed optparse package.

################ R-version 4.0.4
module load R/4.0.4

################ htslib 
module load htslib 

################ bcftools (v1.6)
module load bcftools/1.6 

export script_dir=/scratch.global/balay011/WGS_skChordomas_Bai/SomaticSNV_calling/scripts
export REF_FILE=/scratch.global/balay011/WGS_skChordomas_Bai/reference/Homo_sapiens_assembly38.fasta
export ANNO_FILE=/scratch.global/balay011/WGS_skChordomas_Bai/scripts/WGS_chordomaBai_anno.txt
export DBSNP_FILE=/scratch.global/balay011/WGS_skChordomas_Bai/reference/All_20180418.vcf.gz
export BED_FILE=/scratch.global/balay011/WGS_skChordomas_Bai/reference/Homo_sapiens_assembly38_chrs.bed
export NTHREADS=50
export MUSE_PATH=/home/aventeic/balay011/MuSE/MuSE
export WORK_DIR=/scratch.global/balay011/WGS_skChordomas_Bai/SomaticSNV_calling/work/muse
export CALLER_FILE=/scratch.global/balay011/WGS_skChordomas_Bai/SomaticSNV_calling/scripts/MuseSomaticCaller.r

Rscript $script_dir/runMuseSomaticMutationCalling.r \
--work-dir $WORK_DIR \
--caller-file $CALLER_FILE \
--anno-file $ANNO_FILE \
--ref-file $REF_FILE \
--cpus $NTHREADS \
--dbsnp $DBSNP_FILE \
--muse-path $MUSE_PATH

################ End of Muse run

################ Extract somatic mutations from Muse run
export ANNO_FILE=/home/aventeic/balay011/scripts/WGS_chordomaBai_sommutsfilter.txt
export OUTPUT_DIR=/home/aventeic/balay011/WGS_skChordomas_Bai/SomaticSNV_calling/muse
export PICARD_PATH=/home/aventeic/balay011/.conda/envs/gatk_env/bin/picard

Rscript $script_dir/Call_sommutsfiltering.r \
--script-dir  $script_dir \
--anno-file $ANNO_FILE \
--mode muse \
--output-dir $OUTPUT_DIR \
--picard-path $PICARD_PATH \
--output-vcf 0
