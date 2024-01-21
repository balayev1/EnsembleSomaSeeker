#!/bin/bash

################ This script calls somatic mutations using VarScan package in tumor-normal pair mode on whole-genome sequencing (WGS) samples
################ VarScan v2.4.6 is in usage, 
################ NOTE: you need to install VarScan package in Linux environment. Also all scripts were written in R, so have R and installed optparse package.

# srun -N 1 --cpus-per-task 4 --mem-per-cpu=16gb -t 5:00:00 -p interactive --pty bash

################ I'm using conda environment to use VarScan
module load conda
source activate varscan_env

################ R-version 4.0.4
module load R/4.0.4

################ htslib 
module load htslib 

################ samtools (v1.17)
module load samtools/1.17

################ bcftools (v1.6)
module load bcftools/1.6 

################ Step 1: to pileup somatic mutations for each sample for each chromosome using samtools mpileup toolset
export script_dir=/scratch.global/balay011/WGS_skChordomas_Bai/SomaticSNV_calling/scripts
export SAMTOOLS_PATH=/common/software/install/manual/samtools/1.17_gcc-7.2.0/bin/samtools
export ANNO_FILE=/scratch.global/balay011/WGS_skChordomas_Bai/scripts/WGS_chordomaBai_anno.txt
export REF_FILE=/scratch.global/balay011/WGS_skChordomas_Bai/reference/Homo_sapiens_assembly38.fasta
export BED_FILE=/scratch.global/balay011/WGS_skChordomas_Bai/reference/Homo_sapiens_assembly38_chrs.bed
export VARSCAN_PATH=/home/aventeic/balay011/.conda/envs/varscan_env/bin/varscan
export WORK_DIR=/scratch.global/balay011/WGS_skChordomas_Bai/SomaticSNV_calling/work/varscan
export NTHREADS=8
export CALLER_FILE=/scratch.global/balay011/WGS_skChordomas_Bai/SomaticSNV_calling/scripts/VarScanSomaticCaller.r

Rscript $script_dir/runVarScanSomaticMutationCalling.r \
--work-dir $WORK_DIR \
--caller-file $CALLER_FILE \
--anno-file $ANNO_FILE \
--ref-file $REF_FILE \
--bed-file $BED_FILE \
--cpus $NTHREADS \
--varscan-path $VARSCAN_PATH \
--samtools-path $SAMTOOLS_PATH

################ End of VarScan run

################ Extract somatic mutations from VarScan run
export ANNO_FILE=/home/aventeic/balay011/scripts/WGS_chordomaBai_sommutsfilter.txt
export OUTPUT_DIR=/home/aventeic/balay011/WGS_skChordomas_Bai/SomaticSNV_calling/varscan

Rscript $script_dir/Call_sommutsfiltering.r \
--script-dir  $script_dir \
--anno-file $ANNO_FILE \
--mode varscan \
--output-dir $OUTPUT_DIR \
--output-vcf 0