############################################### Strelka2SomaticCaller.r 
################ This R script contains a function to call somatic mutations using Strelka package in tumor-normal pair mode on whole-genome sequencing (WGS) samples
################ Inputs:
################ @param tumor: full path to tumor BAM file
################ @param normal: full path to normal BAM file
################ @param ref: full path to reference genome GRCh38 FASTA file
################ @param bed_file: full path to bed file with genome coordinates of interest which contains chromosome name, start, end (optional)
################ @param cpus: number of cpus
################ @param outdir: output directory
################ @param strelka2_path: full path to Strelka2 somatic mutation caller file configureStrelkaSomaticWorkflow.py
################ Outputs: see Strelka2 configureStrelkaSomaticWorkflow.py output files (see https://github.com/Illumina/strelka/tree/v2.9.x/docs/userGuide) 

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
    make_option(c("-b", "--bed-file"), type="character", default=NA, 
              help="full path to bed file with genomic coordinates of interest", metavar="character"),
    make_option(c("-c", "--cpus"), type="character", default=1, 
              help="number of cpus", metavar="character"),
    make_option(c("-o", "--output-dir"), type="character", default=NULL, 
              help="output directory", metavar="character"),
    make_option(c("-p", "--strelka-path"), type="character", default=NULL, 
              help="full path to Strelka2 configureStrelkaSomaticWorkflow.py file", metavar="character")); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)


### strelka2_somatic_caller function
strelka2_somatic_caller <- function(tumor_bam, normal_bam, ref_file, bed_file, cpus, outdir, strelka2_path){
    if (is.null(tumor_bam) || is.null(normal_bam) || is.null(ref_file) || is.null(outdir) || is.null(strelka2_path)) {
        stop("Missing required arguments. Please provide all required options.")
    }

    ### Check if output directory exists
    if (!file.exists(outdir)) {
        stop("Output directory does not exist\n")
    }

    cat("============================================================\n")
    cat(format(Sys.time(),usetz = TRUE), "\n")
    start_time <- Sys.time()
    cat("Processing sample pair:", tumor_bam, " ", normal_bam, "\n")
    cat("Running somatic mutation calling using Strelka2\n")

    ### Set strelka2 options
    strelka2.cmdI.options <- paste("--referenceFasta", ref_file, "--runDir", outdir, sep=" ")
    strelka2.cmdII.options <- paste("--mode local --jobs", cpus, "--memGb 64", sep=" ")

    ### add bed file if included
    if (!is.na(bed_file)) {
        paste(strelka2.cmdI.options, "--callRegions", bed_file, sep = " ")
    }

    ### Set strelka2 commands
    strelka2_cmdI <- paste("python", strelka2_path, "--normalBam", normal_bam, "--tumorBam", tumor_bam, strelka2.cmdI.options, sep=" ")
    strelka2_cmdII <- paste("python", file.path(outdir, "runWorkflow.py"), strelka2.cmdII.options, sep=" ")

    ### Execute the strelka2 commands
    tryCatch({
        system(strelka2_cmdI)
        system(strelka2_cmdII)
    }, error = function(e) {
        stop("Error: strelka2 run has failed.\n")
    })

    cat("Finished strelka2 run succesfully\n")
    end_time <- Sys.time()

    ### Print the elapsed time for strelka2 run
    cat("Time elapsed: ", difftime(end_time, start_time, units = "hours"), "hours\n")
}

strelka2_somatic_caller(tumor_bam=opt$'tumor-bam', 
    normal_bam=opt$'normal-bam', 
    ref_file=opt$'ref-file', 
    bed_file=opt$'bed-file',
    cpus = opt$'cpus', 
    outdir=opt$'output-dir', 
    strelka2_path=opt$'strelka-path')
