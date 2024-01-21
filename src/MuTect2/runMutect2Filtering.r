############################################### runMutect2Filtering.r 
################ This R script executes Mutect2FilterSomaticVar.r file to assign filters to somatic mutations using several functions from GATK package in tumor-normal pair mode on whole-genome sequencing (WGS) samples
################ Inputs:
################ @param work-dir: directory where all output of Mutect2 filtering is going to be saved
################ @param caller-file: full path to Mutect2FilterSomaticVar.r
################ @param anno-file: full path to annotation file in TXT format which contains Run ID, tumor status (Yes/No), subject ID and full path to BAM file
################ @param ref-file: full path to reference genome GRCh38 FASTA file
################ @param biallelic-file: full path to file with biallelic loci
################ @param cpus: number of cpus
################ @param mitochondria: whether to filter mitochondrial mutations
################ @param gatk-path: Full path to GATK executable file
################ Output: see GATK MuTect2 output files (see https://gatk.broadinstitute.org/hc/en-us/articles/360036713111-Mutect2#--intervals) 

### install required packages
require("optparse")

### set the arguments
option_list = list(
    make_option(c("-w", "--work-dir"), type="character", default=NULL, 
              help="full path to working directory", metavar="character"),
    make_option(c("-c", "--caller-file"), type="character", default=NULL, 
              help="full path to Mutect2FilterSomaticVar.r file", metavar="character"),
    make_option(c("-a", "--anno-file"), type="character", default=NULL, 
              help="full path to annotation file of samples", metavar="character"),
    make_option(c("-r", "--ref-file"), type="character", default=NULL, 
              help="full path to reference genome GRCh38 FASTA file", metavar="character"),
    make_option(c("-b", "--biallelic-file"), type="character", default=NULL, 
              help="full path to file with biallelic loci", metavar="character"),
    make_option(c("-z", "--cpus"), type="character", default=NULL, 
              help="number of cpus", metavar="character"),
    make_option(c("-m", "--mitochondria"), type="logical", default=TRUE, 
              help="whether to filter mitochondrial mutations", metavar="logical"),
    make_option(c("-g", "--gatk-path"), type="character", default=NULL, 
              help="Full path to GATK executable file", metavar="character")); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

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
ANNO_FILE <- read.table(opt$'anno-file', sep="\t")
colnames(ANNO_FILE) <- c("run_id", "is_tumor", "subject_id", "bam_path")

for (index in 1:length(unique(ANNO_FILE$subject_id))){
    subj <- unique(ANNO_FILE$subject_id)[index]

    anno_sub <- ANNO_FILE[ANNO_FILE$subject_id == subj, ]

    tumor_bam <- anno_sub$bam_path[anno_sub$is_tumor == "Yes"]

    normal_bam <- anno_sub$bam_path[anno_sub$is_tumor == "No"]

    ### create tumor-normal pair subdirectory in sample output directory
    OUTPUT_DIR <- file.path(SAMPLE_DIR, subj)
    if (!dir.exists(OUTPUT_DIR)){
        cat("Creating output directory for subject", subj, "\n")
        dir.create(OUTPUT_DIR)
    }

    ### write the script for Mutect2 filtering run
    #### specify job parameters
    cmd.out <- NULL
    cmd.out <- paste0(cmd.out, "#!/bin/bash\n\n")
    cmd.out <- paste0(cmd.out, "#SBATCH --time=2:00:00 --cpus-per-task=", opt$cpus, " --mem=24G ", 
        " --partition ag2tb", " --mail-type='END,FAIL' --mail-user balay011@umn.edu -A aventeic --job-name=", file.path(SCRIPT_DIR, paste0(subj, "_mutect2filter")), 
        " -o ", file.path(SCRIPT_DIR, paste0(subj, "_mutect2filter.o%J")), " -e ", file.path(SCRIPT_DIR, paste0(subj, "_mutect2filter.e%J\n\n")))
    
    #### write the rest of the script 
    cmd.out <- paste0(cmd.out, "Rscript ", opt$'caller-file', " --tumor-bam ", tumor_bam, " --normal-bam ", normal_bam,
        " --ref-file ", opt$'ref-file', " --biallelic-file ", opt$'biallelic-file',
        " --mitochondria ", opt$'mitochondria',
        " --output-dir ", OUTPUT_DIR,
        " --gatk-path ", opt$'gatk-path')
    
    ### write the script to shell file
    cat(cmd.out, file = file.path(SCRIPT_DIR, paste0(subj, "_mutect2filter.sh")))

    ### execute the job using sbatch 
    # system(paste0("sbatch ", file.path(SCRIPT_DIR, paste0(subj, "_mutect2filter.sh"))))
}