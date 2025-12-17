#!/bin/bash

usage() {
    echo "Usage: $(basename "$0") -t <tumor_bam> -n <normal_bam> -R <ref_genome> -o <dir> [options]"
    echo ""
    echo "Required Arguments:"
    echo "  -t <file>  : Path to the Tumor BAM file."
    echo "  -n <file>  : Path to the Normal BAM file."
    echo "  -R <file>  : Path to the Reference Genome FASTA file."
    echo "  -o <dir>   : Output directory for VCF file."
    echo "  -b <file>  : Path to biallelic VCF file."
    echo ""
    echo "Optional Arguments:"
    echo "  -L <file>  : BED file of coordinates over which to operate."
    echo "  -g <file>  : Germline resource VCF file."
    echo "  -p <file>  : Panel of Normals (PoN) VCF file."
    echo "  -d <dir>   : Temporary directory to use."
    echo "  -S <id>    : Subject/Patient ID used in vcf file name [default: basename tumor_bam]."
    echo "  -C <int>   : Cores for Mutect2 [default: 1]."
    echo "  -M <int>   : Cores for MuSE [default: 1]."
    echo "  -S2 <int>  : Cores for Strelka2 [default: 1]."
    echo "  -V <int>   : Cores for VarScan2 [default: 1]."
    echo "  -Lq <int>  : Cores for Lofreq [default: 1]."
    echo "  -h         : Show help text."
    exit 1
}

## Reset state variables for robust parsing
unset OPTARG
unset OPTIND

## Set locale to operate on ASCII values for characters (disable UTF8)
export LC_ALL=C

## Define variables
TUMOR_BAM=""
NORMAL_BAM=""
REF_GENOME=""
OUTPUT_DIR="" 
BIALLELIC_VCF=""
INTERVALS_BED=""
GERMLINE_RESOURCE=""
PON_FILE=""
TMP_DIR=""
SUBJECT_ID=""
MUTECT2_CORES=1
MUSE_CORES=1
STRELKA2_CORES=1
VARSCAN2_CORES=1
LOFREQ_CORES=1

## getopts processes the script's arguments based on the short options
while getopts ":t:n:R:o:b:L:g:p:C:M:S2:V:Lq:md:S:h" option; do
    case "$option" in
        t) TUMOR_BAM=$OPTARG ;;
        n) NORMAL_BAM=$OPTARG ;;
        R) REF_GENOME=$OPTARG ;;
        o) OUTPUT_DIR=$OPTARG ;;
        b) BIALLELIC_VCF=$OPTARG ;;
        L) INTERVALS_BED=$OPTARG ;;
        g) GERMLINE_RESOURCE=$OPTARG ;;
        p) PON_FILE=$OPTARG ;;
        d) TMP_DIR=$OPTARG ;;
        S) SUBJECT_ID=$OPTARG ;;
        C) MUTECT2_CORES=$OPTARG ;;
        M) MUSE_CORES=$OPTARG ;;
        S2) STRELKA2_CORES=$OPTARG ;;
        V) VARSCAN2_CORES=$OPTARG ;;
        Lq) LOFREQ_CORES=$OPTARG ;;
        h) usage ;;
        \?) printf "Error: Illegal option: -%s\n" "$OPTARG" >&2
            usage
            ;;
    esac
done

## Shifts all arguments processed by getopts to first argument
shift $((OPTIND - 1))

set -ueo pipefail

## Validate required arguments
if [ -z "$TUMOR_BAM" ]
then
    echo "No tumor BAM file specified. To run EnsembleSomaSeeker, tumor BAM file has to be provided."
    exit 1
elif [ ! -f "$TUMOR_BAM" ]
then
    echo "Tumor BAM file $TUMOR_BAM does not exist."
    exit 1
fi

if [ -z "$NORMAL_BAM" ]
then
    echo "No normal BAM file specified. To run EnsembleSomaSeeker, normal BAM file has to be provided."
    exit 1
elif [ ! -f "$NORMAL_BAM" ]
then
    echo "Normal BAM file $NORMAL_BAM does not exist."
    exit 1
fi

if [ -z "$REF_GENOME" ]
then
    echo "No reference genome FASTA file specified. To run EnsembleSomaSeeker, reference genome FASTA file has to be provided."
    exit 1
elif [ ! -f "$REF_GENOME" ]
then
    echo "Reference genome FASTA $REF_GENOME does not exist."
    exit 1
fi

if [ -z "$OUTPUT_DIR" ]
then
    echo "No output directory for VCF file specified. To run EnsembleSomaSeeker, output directory has to be provided."
    exit 1
elif [ ! -d "$OUTPUT_DIR" ]
then
    echo "Output directory $OUTPUT_DIR does not exist."
    exit 1
fi

## Define subject ID if not specified
if [ -z "$SUBJECT_ID" ]; then
    SUBJECT_ID=$(basename "${TUMOR_BAM%.*}")
fi

## Define tumor ID and normal ID
TUMOR_ID=$(basename "${TUMOR_BAM%.*}")
NORMAL_ID=$(basename "${NORMAL_BAM%.*}")

echo "Generating temporary configuration file: $CONFIG_FILE"
cat > "$CONFIG_FILE" << EOF
# Configuration generated for subject: $SUBJECT_ID

reference_genome: "$REF_GENOME"
outputdir: "$OUTPUT_DIR"
tmpdir: "$TMP_DIR"
subject_id: "$SUBJECT_ID"

rule_cores:
    mutect2: $MUTECT2_CORES
    muse: $MUSE_CORES
    strelka2: $STRELKA2_CORES
    varscan2: $VARSCAN2_CORES
    lofreq: $LOFREQ_CORES


sample_data:
    tumor_bam_path: "$TUMOR_BAM"
    normal_bam_path: "$NORMAL_BAM"
    tumor_sample_id: "$TUMOR_ID"
    normal_sample_id: "$NORMAL_ID"

mutect2_params:
    intervals_bed: "$INTERVALS_BED"
    germline_resource: "$GERMLINE_RESOURCE"
    panel_of_normals: "$PON_FILE"

filter_mutect2_params:
    biallelic_resource: "$BIALLELIC_VCF"
EOF

SCRIPT=$(readlink -f "$0")
export ENSEMBLESOMASEEKER=$(dirname "$SCRIPT")

# Final VCF for the single subject ID
TARGET_VCF="$OUTPUT_DIR/mutect2_filtered/${SUBJECT_ID}_filtered.vcf.gz"

echo "Running snakemake pipeline:"
command="snakemake $TARGET_VCF \
    --use-conda --conda-prefix $ENSEMBLESOMASEEKER/envs/conda \
    --cores 1 --configfile $CONFIG \
    --snakefile $ENSEMBLESOMASEEKER/Snakefile -p"

echo -e $command

eval $command


# --- 3. Cleanup ---

# Clean up the single generated config file
rm "$CONFIG_FILE"

echo "Pipeline finished for $SUBJECT_ID. Output VCF: $TARGET_VCF"
exit 0