#!/bin/bash

################# This script takes tab-/comma-delimited file with sample information, and executes EnsembleSomaSeeker.sh as sbatch job for each sample
################# Tab-/comma-delimited file should have the following columns: Path/to/tumor.bam, Path/to/normal.bam, Tumor_ID, Normal_ID, Subject ID (optional)

## REQUIRED INPUTS (Change these paths as needed)
export MANIFEST="/path/to/manifest.tsv"
export REF_GENOME="/projects/standard/aventeic/balay011/references/reference_genome/GRCh38.primary_assembly.genome.fa"
export BIALLELIC_VCF="/projects/standard/aventeic/balay011/references/common_all_biallelic.vcf.gz"
export OUT_BASE="/scratch.global/balay011/EnsembleSomaSeeker_output"

# Check required variables
if [[ -z "$MANIFEST" || -z "$REF_GENOME" || -z "$BIALLELIC_VCF" || -z "$OUT_BASE" ]]; then
    echo "Error: Missing required environment variables."
    exit 1
fi

## OPTIONAL FILE ARGUMENTS (Change these paths as needed or leave empty "")
export INTERVALS_BED="/projects/standard/aventeic/balay011/references/reference_genome/Homo_sapiens_assembly38_chrs.bed"
export GERMLINE_RESOURCE="/users/1/balay011/references/af-only-gnomad.hg38.vcf.gz"
export PON_FILE="/users/1/balay011/references/pon_hiseqx.vcf.gz"
export DBSNP_FILE="/projects/standard/aventeic/balay011/references/dbSNP_files/All_20180418.vcf.gz"
export COSMIC_FILE="/projects/standard/aventeic/balay011/references/Cosmic_files/Cosmic_allvariants.vcf.gz"
export TMP_DIR="/scratch.global/balay011/TMPDIR"

## COMPUTE RESOURCES (Change as needed)
export SNAKEMAKE_CORES=16
export MUTECT2_CORES=8
export MUSE_CORES=16
export STRELKA2_CORES=4
export VARSCAN2_CORES=1
export LOFREQ_CORES=4
export SOMATICSEQ_CORES=8

## READ MANIFEST
DELIM=$'\t'
[[ "$MANIFEST" == *.csv ]] && DELIM=","

# Process manifest (Order: T_BAM, N_BAM, T_ID, N_ID, S_ID)
while IFS="$DELIM" read -r TBAM NBAM TID NID SID || [ -n "$TBAM" ]; do
    
    TBAM=$(echo "$TBAM" | xargs); NBAM=$(echo "$NBAM" | xargs)
    TID=$(echo "$TID" | xargs); NID=$(echo "$NID" | xargs); SID=$(echo "$SID" | xargs)

    [[ -z "$SID" ]] && SID=$(basename "${TBAM%.*}")

    JOB_DIR="$OUT_BASE/$SID"
    mkdir -p "$JOB_DIR/scripts"
    
    # Construct optional file flags
    FILE_ARGS=""
    [[ -n "$INTERVALS_BED" ]]    && FILE_ARGS="$FILE_ARGS -L $INTERVALS_BED"
    [[ -n "$GERMLINE_RESOURCE" ]] && FILE_ARGS="$FILE_ARGS -g $GERMLINE_RESOURCE"
    [[ -n "$PON_FILE" ]]         && FILE_ARGS="$FILE_ARGS -p $PON_FILE"
    [[ -n "$DBSNP_FILE" ]]       && FILE_ARGS="$FILE_ARGS -s $DBSNP_FILE"
    [[ -n "$COSMIC_FILE" ]]      && FILE_ARGS="$FILE_ARGS -v $COSMIC_FILE"
    [[ -n "$TMP_DIR" ]]          && FILE_ARGS="$FILE_ARGS -d $TMP_DIR"

    SBATCH_FILE="$JOB_DIR/scripts/submit_${SID}.sbatch"

    # Generate SBATCH script
    cat << EOF > "$SBATCH_FILE"
#!/bin/bash
#SBATCH --job-name=ESS_$SID
#SBATCH --output=$JOB_DIR/scripts/ESS_%j.out
#SBATCH --error=$JOB_DIR/scripts/ESS_%j.err
#SBATCH --partition=asvnode1,msibigmem,msilarge,msismall
#SBATCH --cpus-per-task=$SNAKEMAKE_CORES
#SBATCH --mem-per-cpu=16G
#SBATCH --time=96:00:00
#SBATCH --account=aventeic
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=balay011@umn.edu

# Run pipeline script
# Note: Ensure EnsembleSomaSeeker.sh is in your PATH or current directory
bash EnsembleSomaSeeker.sh \\
    -t "$TBAM" -n "$NBAM" -T "$TID" -N "$NID" -S "$SID" \\
    -R "$REF_GENOME" -b "$BIALLELIC_VCF" -o "$JOB_DIR" \\
    -j "$SNAKEMAKE_CORES" \\
    -C "$MUTECT2_CORES" -M "$MUSE_CORES" -K "$STRELKA2_CORES" \\
    -V "$VARSCAN2_CORES" -Q "$LOFREQ_CORES" -F "$SOMATICSEQ_CORES" \\
    $FILE_ARGS
EOF

    echo "Submitting job for: $SID"
    sbatch "$SBATCH_FILE"

done < "$MANIFEST"