#!/usr/bin/env bash

############ Script implements NGSCheckMate software package to identify next generation sequencing (NGS) data files from the same individual.

## create conda environment
# mamba create --copy --prefix /projects/standard/aventeic/balay011/.conda/envs/ngscheckmate_env --file ngscheckmate_env.yml

## activate conda environment
source activate /projects/standard/aventeic/balay011/.conda/envs/ngscheckmate_env

## path to directory to store temporary files
TMPDIR="/scratch.global/balay011/TMPDIR"

## path to human reference genome
REF="/projects/standard/aventeic/balay011/references/reference_genome/GRCh38.primary_assembly.genome.fa"

## path to SNP file
SNP="/projects/standard/aventeic/balay011/.conda/envs/ngscheckmate_env/NGSCheckMate/SNP/SNP_GRCh38_hg38_wChr.bed"

## path to input BAM files
BAMDIR="/scratch.global/balay011/WGS_skChordomas_Bai/SomaticSNV_calling/work/bwamem"

## path to working directory
WORKDIR="/scratch.global/balay011/NGSCheckMate/CHORDBAI"
mkdir -p "$WORKDIR"/{mpileup_out,ngscheckmate_out,scripts}

############ Maximize speed to process multiple large BAM files by extracting mutations frim BAM to VCF using samtools mpileup & bcftools call
BAM_LIST_FILE="$WORKDIR/scripts/bam_files.txt"
find "$BAMDIR" -type f -name "*.bam" ! -name "*.bai" | sort > "$BAM_LIST_FILE"
NUM_BAMS=$(wc -l < "$BAM_LIST_FILE")
if [ "$NUM_BAMS" -eq 0 ]; then
    echo "ERROR: No BAM files found in $BAMDIR. Exiting."
    exit 1
fi

echo "Found $NUM_BAMS BAM files."

# Write the SLURM Array Mpileup Script to File
MPILEUP_SCRIPT="$WORKDIR/scripts/mpileup_array.sh"
cat << EOF > "$MPILEUP_SCRIPT"
#!/bin/bash

#SBATCH --job-name=NCM_VCF_Extract
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --output=$WORKDIR/scripts/NCM_VCF.%A_%a.out
#SBATCH --error=$WORKDIR/scripts/NCM_VCF.%A_%a.err
#SBATCH --array=1-$NUM_BAMS
#SBATCH --partition=asvnode1,msilarge,msismall,msibigmem
#SBATCH --account=aventeic
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=balay011@umn.edu

# Define Variables
BAM_LIST_PATH="$BAM_LIST_FILE"
REF_FASTA="$REF"
SNP_BED="$SNP"
OUT_BASE_DIR="$WORKDIR/mpileup_out"
TMP_DIR="$TMPDIR"

# Select BAM file
BAM_FILE=\$(sed -n "\${SLURM_ARRAY_TASK_ID}p" "\$BAM_LIST_PATH")
SAMPLE_NAME=\$(basename "\$BAM_FILE" .bam)
OUT_VCF="\$OUT_BASE_DIR/\$SAMPLE_NAME.vcf"

echo "Processing Sample: \$SAMPLE_NAME (Task ID: \$SLURM_ARRAY_TASK_ID of $NUM_BAMS)"

# Sort BAM by coordinates
echo "Sorting BAM file to temporary location: \$TMP_DIR"
TEMP_SORT_BAM="\$TMP_DIR/\${SAMPLE_NAME}.sorted.bam"
samtools sort -o "\$TEMP_SORT_BAM" "\$BAM_FILE"
if [ $? -ne 0 ]; then
    echo "Error: samtools sort failed for \$BAM_FILE"
    exit 1
fi
echo "Sorting complete."

# Index BAM
if [ ! -f "\$TEMP_SORT_BAM.bai" ]; then
    echo "Indexing BAM file..."
    samtools index "\$TEMP_SORT_BAM"
    if [ \$? -ne 0 ]; then
        echo "Error: samtools index failed."
        exit 1
    fi
fi
echo "Indexing complete."

# Run samtools mpileup piped to bcftools call
echo "Starting VCF generation..."
samtools mpileup \
    -f "\$REF_FASTA" \
    -l "\$SNP_BED" \
    "\$TEMP_SORT_BAM" | \
    bcftools call -c - > "\$OUT_VCF"

if [ \$? -eq 0 ]; then
    echo "SUCCESS: VCF generation complete for \$SAMPLE_NAME."

    echo "Cleaning up temporary files..."
    rm -f "\$TEMP_SORT_BAM"
    rm -f "\$TEMP_SORT_BAM.bai"
    echo "Cleanup complete."
else
    echo "ERROR: VCF generation failed for \$SAMPLE_NAME."
    exit 1
fi
EOF

echo "Submitting mpileup array job to SLURM..."
STEP1_JOB_ID=$(sbatch "$MPILEUP_SCRIPT" | awk '{print $NF}')

echo "SLURM job submitted. Job details will be in $WORKDIR/scripts/"
echo "The list of VCF files will be generated in $WORKDIR/mpileup_out/"

############ NGSCheckMate
# Write the SLURM Array NGSCheckMate Script to File
NGSCHECKMATE_SCRIPT="$WORKDIR/scripts/ngscheckmate_array.sh"
cat << EOF > "$NGSCHECKMATE_SCRIPT"
#!/bin/bash

#SBATCH --job-name=NGSCheckMate_Run
#SBATCH --time=06:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --output=$WORKDIR/scripts/%x.%j.out
#SBATCH --error=$WORKDIR/scripts/%x.%j.err
#SBATCH --partition=asvnode1,msismall,msibigmem,msilarge
#SBATCH --account=aventeic
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=balay011@umn.edu

# Define Variables
SNP_BED="$SNP"

VCF_INPUT_DIR="$WORKDIR/mpileup_out"
OUTPUT_DIR="$WORKDIR/ngscheckmate_out"

# Full path to ncm.py 
NCM_SCRIPT_PATH="/projects/standard/aventeic/balay011/.conda/envs/ngscheckmate_env/NGSCheckMate/ncm.py"

echo "Executing ngscheckmate (ncm.py)..."
echo "VCF Input Directory: \$VCF_INPUT_DIR"
echo "Output Directory: \$OUTPUT_DIR"
echo "NCM Script Path: \$NCM_SCRIPT_PATH"

python "\$NCM_SCRIPT_PATH" \
    --VCF \
    --dir \$VCF_INPUT_DIR \
    --bedfile "\$SNP" \
    --outdir "\$OUTPUT_DIR" \
    --outfilename "ngscheckmate_out"

# Check the exit status of ncm.py
if [ \$? -eq 0 ]; then
    echo "NGSCheckMate run successful!"
    echo "Results saved to: \$OUTPUT_DIR/ngscheckmate_out.txt"
else
    echo "NGSCheckMate run failed."
    exit 1
fi

echo "NGSCheckMate execution complete."
EOF

echo "Submitting NGSCheckMate array job to SLURM - executed once mpileup run is finished"
sbatch --depend=afterok:$STEP1_JOB_ID "$NGSCHECKMATE_SCRIPT"

echo "SLURM job submitted. Job details will be in $WORKDIR/scripts/"
echo "NGSCheckMate output file will be in $WORKDIR/ngscheckmate_out/"