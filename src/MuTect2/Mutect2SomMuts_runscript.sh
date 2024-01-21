#!/bin/bash

################ This script calls and assigns filters to somatic mutations using Mutect2 and other functions from GATK package in tumor-normal pair mode on whole-genome sequencing (WGS) samples
################ Mutect2 (GATK v4.4.0.0) is in usage, 
################ NOTE: you need to install VarScan package in Linux environment. Also all scripts were written in R, so have R and installed optparse package.

# srun -N 1 --cpus-per-task 4 --mem-per-cpu=4gb -t 8:00:00 -p interactive --pty bash

################ I'm using conda environment to use GATK functions
module load conda
source activate gatk_env

################ R-version 4.0.4
module load R/4.0.4

################ bcftools v1.6
module load bcftools/1.6

################ set required variables
export script_dir=/scratch.global/balay011/WGS_skChordomas_Bai/SomaticSNV_calling/scripts
export ANNO_FILE=/scratch.global/balay011/WGS_skChordomas_Bai/scripts/WGS_chordomaBai_anno.txt
export REF_FILE=/scratch.global/balay011/WGS_skChordomas_Bai/reference/Homo_sapiens_assembly38.fasta
export BED_FILE=/scratch.global/balay011/WGS_skChordomas_Bai/reference/Homo_sapiens_assembly38_chrs.bed
export GATK_PATH=/home/aventeic/balay011/.conda/envs/gatk_env/bin/gatk
export WORK_DIR=/scratch.global/balay011/WGS_skChordomas_Bai/SomaticSNV_calling/work/mutect2
export TMPDIR=/scratch.global/balay011/tmp
export CALLER_FILE=/scratch.global/balay011/WGS_skChordomas_Bai/SomaticSNV_calling/scripts/Mutect2SomaticCaller.r
export OMP_NUM_THREADS=8

################ Step 1: call somatic mutations using Mutect2
Rscript $script_dir/runMutect2SomaticMutationCalling.r \
--work-dir $WORK_DIR \
--caller-file $CALLER_FILE \
--anno-file $ANNO_FILE \
--ref-file $REF_FILE \
--bed-file $BED_FILE \
--cpus $OMP_NUM_THREADS \
--mitochondria TRUE \
--tmp-dir $TMPDIR \
--gatk-path $GATK_PATH 

################ End of Mutect2 somatic mutation calling run

################ set required variables
export OMP_NUM_THREADS=1
export CALLER_FILE=/scratch.global/balay011/WGS_skChordomas_Bai/SomaticSNV_calling/scripts/Mutect2FilterSomatic.r
export BIALLELIC_FILE=/home/aventeic/balay011/references/common_all_biallelic.vcf.gz

################ Step 2: assign filters to somatic mutations
Rscript $script_dir/runMutect2Filtering.r \
--work-dir $WORK_DIR \
--caller-file $CALLER_FILE \
--anno-file $ANNO_FILE \
--ref-file $REF_FILE \
--biallelic-file $BIALLELIC_FILE \
--cpus $OMP_NUM_THREADS \
--mitochondria TRUE \
--gatk-path $GATK_PATH 

################ End of Mutect2 somatic mutation filtering run

################ Extract somatic mutations from Mutect2 run
export ANNO_FILE=/home/aventeic/balay011/scripts/WGS_chordomaBai_sommutsfilter.txt
export OUTPUT_DIR=/home/aventeic/balay011/WGS_skChordomas_Bai/SomaticSNV_calling/mutect2
export PICARD_PATH=/home/aventeic/balay011/.conda/envs/gatk_env/bin/picard

Rscript $script_dir/Call_sommutsfiltering.r \
--script-dir  $script_dir \
--anno-file $ANNO_FILE \
--mode mutect2 \
--output-dir $OUTPUT_DIR \
--picard-path $PICARD_PATH \
--output-vcf 0
