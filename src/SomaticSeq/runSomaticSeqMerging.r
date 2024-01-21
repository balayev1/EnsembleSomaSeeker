############################################### runSomaticSeqMerging.r 
################ This R script executes SomaticSeqMerger.r file to merge somatic mutation from different callers using SomaticSeq package in tumor-normal pair mode on whole-genome sequencing (WGS) samples
################ Inputs:
################ @param work-dir: directory where all output of SomaticSeq is going to be saved
################ @param caller-file: Full path to SomaticSeqMerger.r file
################ @param anno-file: full path to annotation file in TXT format which contains full path to tumor BAM file, normal BAM file and VCF files if available (otherwise set 0)
################ @param ref-file: Full path to reference genome GRCh38 FASTA file
################ @param inclusion-region: Full path to bed file with included genomic coordinates which contains chromosome name, start, end (optional)
################ @param exclusion-region: Full path to bed file with excluded genomic coordinates which contains chromosome name, start, end (optional)
################ @param pass-threshold: SCORE for PASS (default: 0.5)
################ @param lowqual-threshold: SCORE for PASS (default: 0.1)
################ @param dbsnp-vcf: Full path to dbSNP VCF file
################ @param cosmic-vcf: Full path to COSMIC VCF file
################ @param threads: Number of threads
################ @param somaticseq-path: Full path to SomaticSeq executable file
################ Outputs: Outputs: see VarScan output files (see https://varscan.sourceforge.net/somatic-calling.html) 

### install required packages
suppressPackageStartupMessages(require(optparse))

### set the arguments
option_list = list(
    make_option(c("-w", "--work-dir"), type="character", 
              help="Full path to working directory", metavar="character"),
    make_option(c("-c", "--caller-file"), type="character", 
              help="Full path to SomaticSeqMerger.r file", metavar="character"),
    make_option(c("-a", "--anno-file"), type="character", 
              help="Full path to annotation file of BAM and VCF files of samples", metavar="character"),
    make_option(c("-r", "--ref-file"), type="character", 
              help="Full path to reference genome GRCh38 FASTA file", metavar="character"),
    make_option(c("-E", "--inclusion-region"), type="character", default=NA, 
              help="Full path to bed file with included genomic coordinates", metavar="character"),
    make_option(c("-e", "--exclusion-region"), type="character", default=NA, 
              help="Full path to bed file with excluded genomic coordinates", metavar="character"),
    make_option(c("-p", "--pass-threshold"), type="numeric", default=0.5, 
              help="SCORE for PASS (default: 0.5)", metavar="numeric"),  
    make_option(c("-q", "--lowqual-threshold"), type="numeric", default=0.1, 
              help="SCORE for LowQual (default: 0.1)", metavar="numeric"),    
    make_option(c("-d", "--dbsnp-vcf"), type="character", default=NA, 
              help="Full path to dbSNP VCF file", metavar="character"),
    make_option(c("-C", "--cosmic-vcf"), type="character", default=NA, 
              help="Full path to COSMIC VCF file", metavar="character"),                
    make_option(c("-t", "--threads"), type="numeric", default=1, 
              help="Number of threads", metavar="numeric"),
    make_option(c("-O", "--somaticseq-path"), type="character", 
              help="Full path to SomaticSeq executable file", metavar="character")); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)
print(opt)
### set working directory
WORK_DIR <- opt$'work-dir'
if (!dir.exists(WORK_DIR)){
    cat("Creating working directory\n")
    dir.create(WORK_DIR)
}

### create script directory
SCRIPT_DIR <- file.path(WORK_DIR, 'scripts')
if (!dir.exists(SCRIPT_DIR)){
    cat("Creating script directory\n")
    dir.create(SCRIPT_DIR)
}

### create sample output directory
SAMPLE_DIR <- file.path(WORK_DIR, 'samples')
if (!dir.exists(SAMPLE_DIR)){
    cat("Creating sample output directory\n")
    dir.create(SAMPLE_DIR)
}

### open the annotation file
ANNO_FILE <- read.table(opt$'anno-file', sep="\t", header=TRUE)
colnames(ANNO_FILE) <- c("subject_id", "tumor_bam", "normal_bam", "mutect2_vcf", "varscan_snv", "varscan_indel", "jsm_vcf", 
    "somaticsniper_vcf", "vardict_vcf", "muse_vcf", "lofreq_snv", "lofreq_indel", "scalpel_vcf", 
    "strelka_snv", "strelka_indel", "arbitrary_snvs", "arbitrary_indels")

for (index in 1:length(unique(ANNO_FILE$subject_id))){
    subj <- unique(ANNO_FILE$subject_id)[index]

    anno_sub <- ANNO_FILE[ANNO_FILE$subject_id == subj, ]
    for (j in 1:length(colnames(anno_sub))){
        if (anno_sub[,j] == 0){
            anno_sub[,j] = NULL
        }
    }

    ### create tumor-normal pair subdirectory in sample output directory
    OUTPUT_DIR <- file.path(SAMPLE_DIR, subj)
    if (!dir.exists(OUTPUT_DIR)){
        cat("Creating output directory for subject", subj, "\n")
        dir.create(OUTPUT_DIR)
    }

    ### write the script for SomaticSeq run
    #### specify job parameters
    cmd.out <- NULL
    cmd.out <- paste0(cmd.out, "#!/bin/bash\n\n")
    cmd.out <- paste0(cmd.out, "#SBATCH --time=10:00:00 --cpus-per-task=", opt$threads, " --mem=128G ", 
        " --partition ag2tb", " --mail-type='END,FAIL' --mail-user balay011@umn.edu -A aventeic --job-name=", file.path(SCRIPT_DIR, paste0(subj, "_somaticseq")), 
        " -o ", file.path(SCRIPT_DIR, paste0(subj, "_somaticseq.o%J")), " -e ", file.path(SCRIPT_DIR, paste0(subj, "_somaticseq.e%J\n\n")))
    

    #### write the rest of the script 
    cmd.out <- paste0(cmd.out, "Rscript ", opt$'caller-file', " --tumor-bam ", anno_sub$tumor_bam, 
            " --normal-bam ", anno_sub$normal_bam,
            " --mutect2-vcf ", anno_sub$mutect2_vcf, " --varscan-snv ", anno_sub$varscan_snv,
            " --varscan-indel ", anno_sub$varscan_indel, " --jsm-vcf ", anno_sub$jsm_vcf, " --somaticsniper-vcf ", anno_sub$somaticsniper_vcf,
            " --vardict-vcf ", anno_sub$vardict_vcf, " --muse-vcf ", anno_sub$muse_vcf, 
            " --lofreq-snv ", anno_sub$lofreq_snv, " --lofreq-indel ", anno_sub$lofreq_indel, 
            " --scalpel-vcf ", anno_sub$scalpel_vcf, " --strelka-snv ", anno_sub$strelka_snv,
            " --strelka-indel ", anno_sub$strelka_indel, " --arbitrary-snvs ", anno_sub$arbitrary_snvs,
            " --arbitrary-indels ", anno_sub$arbitrary_indels, " --output-directory ", OUTPUT_DIR, 
            " --genome-reference ", opt$'ref-file', " --somaticseq-path ", opt$'somaticseq-path', " --threads ", opt$'threads')
    
    #### add included or excluded region file if provided       
    if (!is.na(opt$'inclusion-region') & file.exists(opt$'inclusion-region')){
        cmd.out <- paste0(cmd.out, " --inclusion-region ", opt$'inclusion-region')
    } else {
        cmd.out <- paste0(cmd.out, " --inclusion-region ", NULL)
    }

    if (!is.na(opt$'exclusion-region') & file.exists(opt$'exclusion-region')){
        cmd.out <- paste0(cmd.out, " --exclusion-region ", opt$'exclusion-region')
    }
    
    # add dbSNP and COSMIC files if provided
    if (!is.na(opt$'dbsnp-vcf') & file.exists(opt$'dbsnp-vcf')){
        cmd.out <- paste0(cmd.out, " --dbsnp-vcf ", opt$'dbsnp-vcf')
    }

    if (!is.na(opt$'cosmic-vcf') & file.exists(opt$'cosmic-vcf')){
        cmd.out <- paste0(cmd.out, " --cosmic-vcf ", opt$'cosmic-vcf')
    }

    ### write the script to shell file
    cat(cmd.out, file = file.path(SCRIPT_DIR, paste0(subj, "_somaticseq.sh")))

    ### execute the job using sbatch 
    # system(paste0("sbatch ", file.path(SCRIPT_DIR, paste0(subj, "_somaticseq.sh"))))
}
