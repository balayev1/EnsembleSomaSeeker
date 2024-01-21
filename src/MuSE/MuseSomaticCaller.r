############################################### MuseSomaticCaller.r 
################ This R script contains a function to call somatic mutations using Muse package in tumor-normal pair mode on whole-genome sequencing (WGS) samples
################ Inputs:
################ @param tumor: full path to tumor BAM file
################ @param normal: full path to normal BAM file
################ @param ref: full path to reference genome GRCh38 FASTA file
################ @param outdir: output directory
################ @param cpus: number of cpus
################ @param dbsnp: full path to dbsnp file
################ @param muse_path: full path to MuSE executable file
################ Outputs: see Muse output files (see https://bioinformatics.mdanderson.org/public-software/muse/) 


### install required packages
require("optparse")

### set the arguments
option_list = list(
    make_option(c("-t", "--tumor-bam"), type="character", default=NULL,  
              help="full path to tumor BAM file", metavar="character"),
    make_option(c("-n", "--normal-bam"), type="character", default=NULL, 
              help="full path to normal BAM file", metavar="character"),
    make_option(c("-r", "--ref-file"), type="character", default=NULL, 
              help="full path to reference genome GRCh38 FASTA file", metavar="character"),
    make_option(c("-z", "--cpus"), type="numeric", default=1, 
              help="number of cpus", metavar="numeric"),
    make_option(c("-o", "--output-dir"), type="character", default=NULL, 
              help="output directory", metavar="character"),
    make_option(c("-d", "--dbsnp"), type="character", default=NULL, 
              help="full path to dbsnp file", metavar="character"),
    make_option(c("-p", "--muse-path"), type="character", default=NULL, 
              help="full path to Muse executable file", metavar="character")); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

### muse_somatic_caller function
muse_somatic_caller <- function(tumor, normal, ref, cpus, outdir, dbsnp, muse_path){
    if (is.null(tumor) || is.null(normal) || is.null(ref) || is.null(outdir) || is.null(muse_path)) {
        stop("Missing required arguments. Please provide all required options.")
    }

    ### Check if output directory exists
    if (!file.exists(outdir)) {
        cat("Output directory does not exist. Creating ...\n")
        dir.create(outdir)
    }

    cat("============================================================\n")
    cat(format(Sys.time(),usetz = TRUE), "\n")
    cat("Processing sample pair:", tumor, " ", normal, "\n")
    cat("Running somatic mutation calling using Muse\n")

    ### Write the script for Muse somatic mutation caller execution
    cmd_call <- paste(muse_path, "call -f", ref, "-O", file.path(outdir, basename(outdir)), "-n", cpus, tumor, normal, sep = " ")
    cmd_sump <- paste(muse_path, "sump -I", file.path(outdir, paste0(basename(outdir), ".MuSE.txt")), 
        "-O", file.path(outdir, paste0(basename(outdir),".MuSE.vcf")), "-G -n", cpus, "-D", dbsnp, sep = " ")
    
    ### Error handler script
    tryCatch({
        system(cmd_call)
    }, error = function(e) {
        stop("Muse call command failed.\n")
    })
    
    tryCatch({
        system(cmd_sump)
    }, error = function(e) {
        stop("Muse sump command failed.\n")
    }) 

    cat(format(Sys.time(),usetz = TRUE), "\n")
    cat("Done. Finished somatic mutation calling succesfully!!!")
}

muse_somatic_caller(tumor = opt$'tumor-bam',
    normal = opt$'normal-bam',
    ref = opt$'ref-file',
    cpus = opt$'cpus',
    outdir = opt$'output-dir',
    dbsnp = opt$'dbsnp',
    muse_path = opt$'muse-path')

