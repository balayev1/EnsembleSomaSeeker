#!/bin/bash

################ This script merges somatic mutation from different callers using SomaticSeq package in tumor-normal pair mode on whole-genome sequencing (WGS) samples
################ SomaticSeq v3.7.4 is in usage, 
################ NOTE: you need to install SomaticSeq package in Linux environment. Also all scripts were written in R, so have R and installed optparse package.

# srun -N 1 --cpus-per-task 4 --mem-per-cpu=16gb -t 5:00:00 -p interactive --pty bash

################ I'm using conda environment to use VarScan
module load conda
source activate somaticseq_env

################ R-version 4.0.4
module load R/4.0.4

################ Set required variables
export script_dir=/scratch.global/balay011/WGS_skChordomas_Bai/SomaticSNV_calling/scripts
export ANNO_FILE=/home/aventeic/balay011/scripts/WGS_chordomaBai_somaticseq.txt
export DBSNP_FILE=/home/aventeic/balay011/references/All_20180418.vcf.gz
export COSMIC_FILE=/home/aventeic/balay011/references/Cosmic_allvariants.vcf.gz
export REF_FILE=/scratch.global/balay011/WGS_skChordomas_Bai/reference/Homo_sapiens_assembly38.fasta
export BED_FILE=/scratch.global/balay011/WGS_skChordomas_Bai/reference/Homo_sapiens_assembly38_chrs.bed
export SOMATICSEQ_PATH=/home/aventeic/balay011/somaticseq/somaticseq/somaticseq_parallel.py
export WORK_DIR=/scratch.global/balay011/WGS_skChordomas_Bai/SomaticSNV_calling/work/somaticseq
export NTHREADS=16
export CALLER_FILE=/scratch.global/balay011/WGS_skChordomas_Bai/SomaticSNV_calling/scripts/SomaticSeqMerger.r

Rscript $script_dir/runSomaticSeqMerging.r \
--work-dir $WORK_DIR \
--caller-file $CALLER_FILE \
--anno-file $ANNO_FILE \
--ref-file $REF_FILE \
--inclusion-region $BED_FILE \
--dbsnp-vcf $DBSNP_FILE \
--cosmic-vcf $COSMIC_FILE \
--threads $NTHREADS \
--somaticseq-path $SOMATICSEQ_PATH

################ End of SomaticSeq run