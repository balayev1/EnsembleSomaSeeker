############################################### runMutect2SomaticMutationCalling.r 
################ This R script executes Mutect2SomaticCaller.r file to call somatic mutations using Mutect2 function from GATK package in tumor-normal pair mode on whole-genome sequencing (WGS) samples
################ Inputs:
################ @param work-dir: directory where all output of Mutect2 is going to be saved
################ @param caller-file: full path to Mutect2SomaticCaller.r
################ @param anno-file: full path to annotation file in TXT format which contains Run ID, tumor status (Yes/No), subject ID and full path to BAM file
################ @param ref-file: full path to reference genome GRCh38 FASTA file
################ @param bed-file: full path to bed file with genome coordinates of interest which contains chromosome name, start, end (optional)
################ @param: germ-file: full path to germline resource file
################ @param pon: full path to panel of normals file
################ @param cpus: number of cpus
################ @param mitochondria: whether to call mitochondrial mutations
################ @param tmp-dir: temporary directory
################ @param gatk-path: Full path to GATK executable file
################ Output: see GATK MuTect2 output files (see https://gatk.broadinstitute.org/hc/en-us/articles/360036713111-Mutect2#--intervals) 

### install required packages
require("optparse")

### set the arguments
option_list = list(
    make_option(c("-w", "--work-dir"), type="character", default=NULL, 
              help="full path to working directory", metavar="character"),
    make_option(c("-c", "--caller-file"), type="character", default=NULL, 
              help="full path to Mutect2SomaticCaller.r file", metavar="character"),
    make_option(c("-a", "--anno-file"), type="character", default=NULL, 
              help="full path to annotation file of samples", metavar="character"),
    make_option(c("-r", "--ref-file"), type="character", default=NULL, 
              help="full path to reference genome GRCh38 FASTA file", metavar="character"),
    make_option(c("-b", "--bed-file"), type="character", default=NULL, 
        help="full path to bed file with genomic coordinates of interest", metavar="character"),
    make_option(c("-s", "--germ-file"), type="character", default=NULL, 
              help="full path to germline resource file", metavar="character"),
    make_option(c("-d", "--pon"), type="character", default=NULL, 
              help="full path to panel of normals file", metavar="character"),
    make_option(c("-z", "--cpus"), type="integer", default=1, 
              help="number of cpus", metavar="integer"),
    make_option(c("-m", "--mitochondria"), type="logical", default=TRUE, 
              help="whether to call mitochondrial mutations", metavar="logical"),
    make_option(c("-p", "--tmp-dir"), type="character", default=NULL, 
              help="temporary directory", metavar="character"),
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

    tumor_id <- anno_sub$run_id[anno_sub$is_tumor == "Yes"]

    normal_id <- anno_sub$run_id[anno_sub$is_tumor == "No"]

    tumor_bam <- anno_sub$bam_path[anno_sub$is_tumor == "Yes"]

    normal_bam <- anno_sub$bam_path[anno_sub$is_tumor == "No"]

    ### create tumor-normal pair subdirectory in sample output directory
    OUTPUT_DIR <- file.path(SAMPLE_DIR, subj)
    if (!dir.exists(OUTPUT_DIR)){
        cat("Creating output directory for subject", subj, "\n")
        dir.create(OUTPUT_DIR)
    }

    ### check if tmpdir exists
    if (!is.null(opt$'tmp-dir')){
        if (file.exists(opt$'tmp-dir')){
            tmpdir <- file.path(opt$'tmp-dir', paste(tumor_id, normal_id, sep = "_"))
            dir.create(tmpdir)
        } else {
            tmpdir <- NULL
            cat("TMPDIR does not exist\n")
        }
    }

    ### write the script for Mutect2 run
    #### specify job parameters
    cmd.out <- NULL
    cmd.out <- paste0(cmd.out, "#!/bin/bash\n\n")
    cmd.out <- paste0(cmd.out, "#SBATCH --time=96:00:00 --cpus-per-task=", opt$cpus, " --mem=256G ", 
        " --partition ag2tb", " --mail-type='END,FAIL' --mail-user balay011@umn.edu -A aventeic --job-name=", file.path(SCRIPT_DIR, paste0(subj, "_mutect2")), 
        " -o ", file.path(SCRIPT_DIR, paste0(subj, "_mutect2.o%J")), " -e ", file.path(SCRIPT_DIR, paste0(subj, "_mutect2.e%J\n\n")))
    
    #### write the rest of the script 
    cmd.out <- paste0(cmd.out, "Rscript ", opt$'caller-file', " --tumor-bam ", tumor_bam, " --normal-bam ", normal_bam,
        " --tumor-id ", tumor_id, " --normal-id ", normal_id, " --cpus ", opt$'cpus', " --mitochondria ", opt$'mitochondria', 
        " --tmp-dir ", tmpdir, 
        " --ref-file ", opt$'ref-file', " --output-dir ", OUTPUT_DIR,
        " --gatk-path ", opt$'gatk-path')

    if (!is.null(opt$'bed-file')){
        cmd.out <- paste(cmd.out, "--bed-file", opt$'bed-file', sep = " ")
    }

    if (!is.null(opt$'germ-file')){
        cmd.out <- paste(cmd.out, "--germ-file", opt$'germ-file', sep = " ")
    }

    if (!is.null(opt$'pon')){
        cmd.out <- paste(cmd.out, "--pon", opt$'pon', sep = " ")
    }
    
    ### write the script to shell file
    cat(cmd.out, file = file.path(SCRIPT_DIR, paste0(subj, "_mutect2.sh")))

    ### execute the job using sbatch 
    # system(paste0("sbatch ", file.path(SCRIPT_DIR, paste0(subj, "_mutect2.sh"))))
}









