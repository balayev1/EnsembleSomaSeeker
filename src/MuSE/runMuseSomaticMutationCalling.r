############################################### runMuseSomaticMutationCalling.r 
################ This R script executes MuseSomaticCaller.r file to call somatic mutations using Muse package in tumor-normal pair mode on whole-genome sequencing (WGS) samples
################ Inputs:
################ @param work-dir: directory where all output of Muse is going to be saved
################ @param caller-dir: full path to MuseSomaticCaller.r
################ @param anno-dir: full path to annotation file in TXT format which contains Run ID, tumor status (Yes/No), subject ID and full path to BAM file
################ @param ref-dir: full path to reference genome GRCh38 FASTA file
################ @param cpus: number of cpus
################ @param dbsnp: full path to dbsnp file
################ @param muse-path: full path to Muse executable file
################ Outputs: Outputs: see Muse output files (see https://bioinformatics.mdanderson.org/public-software/muse/) 

### install required packages
require("optparse")

### set the arguments
option_list = list(
    make_option(c("-w", "--work-dir"), type="character", default=NULL, 
              help="full path to working directory", metavar="character"),
    make_option(c("-c", "--caller-file"), type="character", default=NULL, 
              help="full path to MuseSomaticCaller.r file", metavar="character"),
    make_option(c("-a", "--anno-file"), type="character", default=NULL, 
              help="full path to annotation file of samples", metavar="character"),
    make_option(c("-r", "--ref-file"), type="character", default=NULL, 
              help="full path to reference genome GRCh38 FASTA file", metavar="character"),
    make_option(c("-z", "--cpus"), type="character", default=NULL, 
              help="number of cpus", metavar="character"),
    make_option(c("-d", "--dbsnp"), type="character", default=NULL, 
              help="full path to dbsnp file", metavar="character"),
    make_option(c("-p", "--muse-path"), type="character", default=NULL, 
              help="full path to Muse executable file", metavar="character")); 

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
    ### specify job commands
    cmd.out <- NULL
    cmd.out <- paste0(cmd.out, "#!/bin/bash\n\n")
    cmd.out <- paste0(cmd.out, "#SBATCH --time=24:00:00 --cpus-per-task=", opt$cpus, " --mem-per-cpu=16G ", 
        " --partition ag2tb", " --mail-type='END,FAIL' --mail-user balay011@umn.edu -A aventeic --job-name=", 
        file.path(SCRIPT_DIR, paste0(subj, "_muse")), 
        " -o ", file.path(SCRIPT_DIR, paste0(subj, "_muse.o%J")), " -e ", file.path(SCRIPT_DIR, paste0(subj, "_muse.e%J\n\n")))

    ### write the rest of the script
        cmd.out <- paste0(cmd.out, "Rscript ", opt$'caller-file', " --tumor-bam ", tumor_bam, " --normal-bam ", normal_bam,
    " --ref-file ", opt$'ref-file', " --cpus ", opt$cpus, " --output-dir ", OUTPUT_DIR, " --dbsnp ", opt$dbsnp, 
        " --muse-path ", opt$'muse-path', "\n")

    ### write the script to shell file
    cat(cmd.out, file = file.path(SCRIPT_DIR, paste0(subj, "_muse.sh")))

    ### execute the job using sbatch 
    system(paste0("sbatch ", file.path(SCRIPT_DIR, paste0(subj, "_muse.sh"))))
}