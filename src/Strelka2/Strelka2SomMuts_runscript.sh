#!/bin/bash

################ This script calls somatic mutations using Strelka package in tumor-normal pair mode on whole-genome sequencing (WGS) samples
################ Strelka v2.9.10 is in usage, 
################ NOTE: you need to install Strelka package in Linux environment. Also all scripts were written in R, so have R and installed optparse package.

# srun -N 1 --cpus-per-task 4 --mem-per-cpu=4gb -t 8:00:00 -p interactive --x11 --pty bash

################ I'm using conda environment to use Strelka
module load conda
source activate strelka_env

################ R-version 4.0.4
module load R/4.0.4

################ bcftools v1.6
module load bcftools/1.6

export script_dir=/scratch.global/balay011/WGS_skChordomas_Bai/SomaticSNV_calling/scripts
export WORK_DIR=/scratch.global/balay011/WGS_skChordomas_Bai/SomaticSNV_calling/work/strelka2
export CALLER_FILE=/scratch.global/balay011/WGS_skChordomas_Bai/SomaticSNV_calling/scripts/Strelka2SomaticCaller.r
export ANNO_FILE=/scratch.global/balay011/WGS_skChordomas_Bai/scripts/WGS_chordomaBai_anno.txt
export DBSNP_FILE=/scratch.global/balay011/WGS_skChordomas_Bai/reference/All_20180418.vcf.gz
export REF_FILE=/scratch.global/balay011/WGS_skChordomas_Bai/reference/Homo_sapiens_assembly38.fasta
export BED_FILE=/scratch.global/balay011/WGS_skChordomas_Bai/reference/Homo_sapiens_assembly38_chrs.bed.gz
export STRELKA_PATH=/home/aventeic/balay011/.conda/envs/strelka_env/bin
export SAMTOOLS_PATH=/common/software/install/manual/samtools/1.17_gcc-7.2.0/bin/samtools
export NTHREADS=16

Rscript $script_dir/runStrelka2SomaticMutationCalling.r \
--work-dir $WORK_DIR \
--caller-file $CALLER_FILE \
--anno-file $ANNO_FILE \
--ref-file $REF_FILE \
--bed-file $BED_FILE \
--cpus $NTHREADS \
--strelka-path $STRELKA_PATH/configureStrelkaSomaticWorkflow.py

################ End of Strelka run

################ Extract somatic mutations from Strelka run
export ANNO_FILE=/home/aventeic/balay011/scripts/WGS_chordomaBai_sommutsfilter.txt
export OUTPUT_DIR=/home/aventeic/balay011/WGS_skChordomas_Bai/SomaticSNV_calling/strelka2

Rscript $script_dir/Call_sommutsfiltering.r \
--script-dir  $script_dir \
--anno-file $ANNO_FILE \
--mode strelka2 \
--output-dir $OUTPUT_DIR \
--output-vcf 0




