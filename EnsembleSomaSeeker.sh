#!/bin/bash

usage() {
    echo "Usage: $(basename "$0") -t <tumor_bam> -n <normal_bam> -R <ref_genome> -o <dir> [options]"
    echo ""
    echo "Required Arguments:"
    echo "  -t <file>  : Path to the Tumor BAM file."
    echo "  -n <file>  : Path to the Normal BAM file."
    echo "  -T <id>    : Tumor Sample ID."
    echo "  -N <id>    : Normal Sample ID."
    echo "  -R <file>  : Path to the Reference Genome FASTA file."
    echo "  -o <dir>   : Output directory for VCF file."
    echo "  -b <file>  : Path to biallelic VCF file."
    echo ""
    echo "Optional Arguments:"
    echo "  -L <file>  : Path to BED file of coordinates over which to operate."
    echo "  -g <file>  : Path to Germline resource VCF file (used only for Mutect2)."
    echo "  -p <file>  : Path to Panel of Normals (PoN) VCF file (used only for Mutect2)."
    echo "  -P <file>  : Path to high-confidence PoN VCF file for variant exclusion (e.g. 1000G PoN see GATK resource bundle)."
    echo "  -G <file>  : Path to Germline SNV VCF from Normal BAM file for variant exclusion."
    echo "  -H <file>  : Path to Germline INDEL VCF from Normal BAM file for variant exclusion."
    echo "  -s <file>  : Path to dbSNP VCF file."
    echo "  -v <file>  : Path to COSMIC VCF file (used only for SomaticSeq)."
    echo "  -d <dir>   : Temporary directory to use."
    echo "  -S <id>    : Subject/Patient ID used in vcf file name [default: basename tumor_bam]."
    echo "  -j <int>   : Cores for snakemake run [default: 1]."
    echo "  -C <int>   : Cores for Mutect2 [default: 1]."
    echo "  -M <int>   : Cores for MuSE [default: 1]."
    echo "  -K <int>   : Cores for Strelka2 [default: 1]."
    echo "  -V <int>   : Cores for VarScan2 [default: 1]."
    echo "  -Q <int>   : Cores for Lofreq [default: 1]."
    echo "  -F <int>   : Cores for SomaticSeq [default: 1]."
    echo "  -A <int>   : Min normal depth for variant exclusion [default: 10]"
    echo "  -B <int>   : Min tumor depth for variant exclusion [default: 15]"
    echo "  -E <float> : Min tumor allele frequency for variant exclusion [default: 0.05]"
    echo "  -I <float> : Max normal allele frequency for variant exclusion [default: 0.02]"
    echo "  -J <int>   : Min tumor ALT concordant reads for variant exclusion [default: 5]"
    echo "  -X <int>   : Max normal ALT concordant reads for variant exclusion [default: 2]"
    echo "  -Y <int>   : Min mean mapping quality for variant exclusion [default: 23]"
    echo "  -Z <int>   : Min mean base quality for variant exclusion [default: 23]"
    echo "  -W <int>   : Min tumor ALT mapping quality for variant exclusion [default: 30]"
    echo "  -U <int>.  : Min number of supporting callers for variant exclusion [default: 2]"
    echo "  -h         : Show help text."
    exit 1
}

## Reset state variables for robust parsing
unset OPTARG
unset OPTIND

## Set locale to operate on ASCII values for characters (disable UTF8)
export LC_ALL=C

## Initialize variables with default values
TUMOR_BAM=""
NORMAL_BAM=""
TUMOR_ID=""
NORMAL_ID=""
REF_GENOME=""
OUTPUT_DIR="" 
BIALLELIC_VCF=""
INTERVALS_BED=""
GERMLINE_RESOURCE=""
PON_FILE=""
HC_VCF=""
GERMLINE_SNV_VCF=""
GERMLINE_INDEL_VCF=""
DBSNP_FILE=""
COSMIC_FILE=""
TMP_DIR=""
SUBJECT_ID=""
SNAKEMAKE_CORES=1
MUTECT2_CORES=1
MUSE_CORES=1
STRELKA2_CORES=1
VARSCAN2_CORES=1
LOFREQ_CORES=1
SOMATICSEQ_CORES=1
MIN_N_DP=10
MIN_T_DP=15
MIN_T_AF=0.05
MAX_N_AF=0.02
MIN_T_CONC=5
MAX_N_CONC=2
MIN_MQ=23
MIN_BQ=23
MIN_T_ALT_MQ=30
MIN_CALLERS=2

## Updated getopts string to match the CASE logic below
# Removed 'c', 'm', 'k', 'q', 'f' and replaced with 'C', 'M', 'K', 'Q', 'F'
while getopts ":t:n:T:N:R:o:b:L:g:p:P:G:H:s:v:d:S:j:C:M:K:V:Q:F:A:B:E:I:J:X:Y:Z:W:U:h" option; do
    case "$option" in
        t) TUMOR_BAM=$OPTARG ;;
        n) NORMAL_BAM=$OPTARG ;;
        T) TUMOR_ID=$OPTARG ;;
        N) NORMAL_ID=$OPTARG ;;
        R) REF_GENOME=$OPTARG ;;
        o) OUTPUT_DIR=$OPTARG ;;
        b) BIALLELIC_VCF=$OPTARG ;;
        L) INTERVALS_BED=$OPTARG ;;
        g) GERMLINE_RESOURCE=$OPTARG ;;
        p) PON_FILE=$OPTARG ;;
        P) HC_VCF=$OPTARG ;;
        G) GERMLINE_SNV_VCF=$OPTARG ;;
        H) GERMLINE_INDEL_VCF=$OPTARG ;;
        s) DBSNP_FILE=$OPTARG ;;
        v) COSMIC_FILE=$OPTARG ;;
        d) TMP_DIR=$OPTARG ;;
        S) SUBJECT_ID=$OPTARG ;;
        j) SNAKEMAKE_CORES=$OPTARG ;;
        C) MUTECT2_CORES=$OPTARG ;;
        M) MUSE_CORES=$OPTARG ;;
        K) STRELKA2_CORES=$OPTARG ;;
        V) VARSCAN2_CORES=$OPTARG ;;
        Q) LOFREQ_CORES=$OPTARG ;;
        F) SOMATICSEQ_CORES=$OPTARG ;;
        A) MIN_N_DP=$OPTARG ;;
        B) MIN_T_DP=$OPTARG ;;
        E) MIN_T_AF=$OPTARG ;;
        I) MAX_N_AF=$OPTARG ;;
        J) MIN_T_CONC=$OPTARG ;;
        X) MAX_N_CONC=$OPTARG ;;
        Y) MIN_MQ=$OPTARG ;;
        Z) MIN_BQ=$OPTARG ;;
        W) MIN_T_ALT_MQ=$OPTARG ;;
        U) MIN_CALLERS=$OPTARG ;;
        h) usage ;;
        \?) printf "Error: Illegal option: -%s\n" "$OPTARG" >&2
            usage
            ;;
    esac
done

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

if [ -z "$TUMOR_ID" ]
then
    echo "No tumor ID specified. To run EnsembleSomaSeeker, tumor ID has to be provided."
    exit 1
fi

if [ -z "$NORMAL_ID" ]
then
    echo "No normal ID specified. To run EnsembleSomaSeeker, normal ID has to be provided."
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

if [ -z "$BIALLELIC_VCF" ]
then
    echo "No biallelic VCF file $BIALLELIC_VCF for Mutect2 specified."
    exit 1
elif [ ! -f "$BIALLELIC_VCF" ]
then
    echo "Biallelic VCF file $BIALLELIC_VCF for Mutect2 does not exist."
    exit 1
fi

## Validate optional file arguments if not empty
if [ ! -z "$INTERVALS_BED" ] && [ ! -f "$INTERVALS_BED" ]
then
    echo "Intervals BED file $INTERVALS_BED specified but does not exist."
    exit 1
fi

if [ ! -z "$DBSNP_FILE" ] && [ ! -f "$DBSNP_FILE" ]
then
    echo "dbSNP VCF file $DBSNP_FILE specified but does not exist."
    exit 1
fi

if [ ! -z "$COSMIC_FILE" ] && [ ! -f "$COSMIC_FILE" ]
then
    echo "COSMIC VCF file $COSMIC_FILE specified but does not exist."
    exit 1
fi

if [ ! -z "$GERMLINE_RESOURCE" ] && [ ! -f "$GERMLINE_RESOURCE" ]
then
    echo "Germline resource VCF file $GERMLINE_RESOURCE for Mutect2 specified but does not exist."
    exit 1
fi

if [ ! -z "$PON_FILE" ] && [ ! -f "$PON_FILE" ]
then
    echo "PoN VCF file $PON_FILE for Mutect2 specified but does not exist."
    exit 1
fi

## Remove any trailing slashes from output directory
OUTPUT_DIR=$(echo "$OUTPUT_DIR" | sed 's:/*$::')

## Define subject ID if not specified
if [ -z "$SUBJECT_ID" ]; then
    SUBJECT_ID=$(basename "${TUMOR_BAM%.*}")
fi

## Define current path
SCRIPT=$(readlink -f "$0")
export ENSEMBLESOMASEEKER=$(dirname "$SCRIPT")

## Define path to config file
CONFIG_FILE="$ENSEMBLESOMASEEKER/config_${SUBJECT_ID}.yaml"

echo "Generating temporary configuration file: $CONFIG_FILE"
cat > "$CONFIG_FILE" << EOF
# Configuration generated for subject: $SUBJECT_ID

reference_genome: "$REF_GENOME"
outputdir: "$OUTPUT_DIR"
tmpdir: "$TMP_DIR"
subject_id: "$SUBJECT_ID"
intervals_bed: "$INTERVALS_BED"
dbsnp_resource: "$DBSNP_FILE"
cosmic_resource: "$COSMIC_FILE"
hc_pon: "$HC_VCF"
germline_snv_resource: "$GERMLINE_SNV_VCF"
germline_indel_resource: "$GERMLINE_INDEL_VCF"

rule_cores:
    mutect2: $MUTECT2_CORES
    muse: $MUSE_CORES
    strelka2: $STRELKA2_CORES
    varscan2: $VARSCAN2_CORES
    lofreq: $LOFREQ_CORES
    somaticseq: $SOMATICSEQ_CORES

sample_data:
    tumor_bam_path: "$TUMOR_BAM"
    normal_bam_path: "$NORMAL_BAM"
    tumor_sample_id: "$TUMOR_ID"
    normal_sample_id: "$NORMAL_ID"


mutect2_params:
    germline_resource: "$GERMLINE_RESOURCE"
    panel_of_normals: "$PON_FILE"


filter_mutect2_params:
    biallelic_resource: "$BIALLELIC_VCF"

variant_exclusion_params:
    n_dp: $MIN_N_DP
    t_dp: $MIN_T_DP
    t_af: $MIN_T_AF
    n_af: $MAX_N_AF
    t_conc: $MIN_T_CONC
    n_conc: $MAX_N_CONC
    mq: $MIN_MQ
    bq: $MIN_BQ
    t_alt_mq: $MIN_T_ALT_MQ
    min_callers: $MIN_CALLERS
EOF

# Final VCFs for subject ID
TARGET_SNV_VCF="${OUTPUT_DIR}/final_vcf/${SUBJECT_ID}_somatic_SNV.vcf.gz"
TARGET_INDEL_VCF="${OUTPUT_DIR}/final_vcf/${SUBJECT_ID}_somatic_INDEL.vcf.gz"

echo "Running snakemake pipeline:"
command="snakemake $TARGET_SNV_VCF $TARGET_INDEL_VCF \
    --use-conda --conda-prefix $ENSEMBLESOMASEEKER/envs/conda \
    --cores $SNAKEMAKE_CORES --configfile $CONFIG_FILE \
    --snakefile $ENSEMBLESOMASEEKER/Snakefile -p"

echo -e $command

eval $command

# Clean up the config file
rm "$CONFIG_FILE"

echo "EnsembleSomaSeeker finished for $SUBJECT_ID."
echo " Somatic SNVs: $TARGET_SNV_VCF"
echo " Somatic INDELs: $TARGET_INDEL_VCF"

exit 0