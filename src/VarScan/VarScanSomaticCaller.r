############################################### VarScanSomaticCaller.r 
################ This R script contains a function to call somatic mutations using VarScan package in tumor-normal pair mode on whole-genome sequencing (WGS) samples
################ Inputs:
################ @param tumor-bam: Full path to tumor BAM file
################ @param normal-bam: Full path to normal BAM file
################ @param tumor-id: Tumor sample ID
################ @param normal-id: Normal sample ID
################ @param: ref-file: Full path to reference genome FASTA file
################ @param: bed-file: Full path to bed file with genomic coordinates of interest
################ @param cpus: number of cpus
################ @param outdir: Output directory
################ @param samtools-path: Full path to samtools executable file
################ @param varscan_path: Full path to VarScan executable file
################ Outputs: see VarScan output files (see https://varscan.sourceforge.net/somatic-calling.html) 

### install required packages
require("optparse")

### set the arguments
option_list = list(
    make_option(c("-T", "--tumor-bam"), type="character", default=NULL, 
              help="Full path to tumor BAM file", metavar="character"),
    make_option(c("-N", "--normal-bam"), type="character", default=NULL, 
              help="Full path to normal BAM file", metavar="character"),
    make_option(c("-t", "--tumor-id"), type="character", default=NULL, 
              help="Tumor sample ID", metavar="character"),
    make_option(c("-n", "--normal-id"), type="character", default=NULL, 
              help="Normal sample ID", metavar="character"),
    make_option(c("-r", "--ref-file"), type="character", default=NULL, 
              help="Full path to reference genome FASTA file", metavar="character"),
    make_option(c("-e", "--bed-file"), type="character", default=NA, 
        help="Full path to bed file with genomic coordinates of interest", metavar="character"),
    make_option(c("-z", "--cpus"), type="character", default=1, 
              help="number of cpus", metavar="character"),
    make_option(c("-o", "--output-dir"), type="character", default=NULL, 
              help="Output directory", metavar="character"),
    make_option(c("-s", "--samtools-path"), type="character", default=NULL, 
              help="Full path to samtools executable file", metavar="character"),
    make_option(c("-v", "--varscan-path"), type="character", default=NULL, 
              help="Full path to VarScan executable file", metavar="character")); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

### varscan_somatic_caller function
varscan_somatic_caller <- function(tumor_bam, normal_bam, tumor_id, normal_id, ref_file, bed_file, cpus, outdir, samtools_path, varscan_path){
    if (is.null(tumor_bam) || is.null(normal_bam) || is.null(tumor_id) || is.null(normal_id) || is.null(ref_file) || is.null(outdir) || is.null(samtools_path) || is.null(varscan_path)) {
        stop("Missing required arguments. Please provide all required options.")
    }

    ### Check if output directory exists
    if (!file.exists(outdir)) {
        cat("Output directory does not exist. Creating ...\n")
        dir.create(outdir)
    }

    cat("============================================================\n")
    cat(format(Sys.time(),usetz = TRUE), "\n")
    start_time <- Sys.time()
    cat("Processing sample pair:", tumor_id, " ", normal_id, "\n")
    cat("Beginning to run mpileup from samtools\n")

    ### Set mpileup options
    mpileup.output.file <- file.path(outdir, paste0(tumor_id, "_", normal_id, ".mpileup.gz"))
    mpileup.options <- paste("--fasta-ref", ref_file, "--min-MQ 1", "--no-BAQ", sep = " ")

    ### add bed file if included
    if (!is.na(bed_file)) {
        paste(mpileup.options, "--positions", bed_file, sep = " ")
    }

    ### Set mpileup command
    mpileup_cmd <- paste(samtools_path, "mpileup", mpileup.options, normal_bam, tumor_bam, "| bgzip --threads",
        cpus, ">", mpileup.output.file, 
        sep = " ")

    ### Execute the mpileup command
    tryCatch({
        system(mpileup_cmd)
    }, error = function(e) {
        stop("Error: mpileup run has failed.\n")
    })

    cat("Finished mpileup run succesfully\n")
    end_time <- Sys.time()

    ### Print the elapsed time for mpileup run
    cat("Time elapsed: ", difftime(end_time, start_time, units = "hours"), "hours\n")

    start_time <- Sys.time()
    ### Set Varscan options
    varscan.options <- paste("--mpileup 1 --min-coverage 1 --min-var-freq 0.001 --output-vcf 1", sep = " ")

    ### Set Varscan command
    varscan_cmd <- paste(varscan_path, "somatic", paste0("<(zcat ", mpileup.output.file, ")"), file.path(outdir, basename(outdir)),
        varscan.options, sep = " ")

    ### Compress output vcf files from Varscan and index them
    output.snv.vcf <- file.path(outdir, paste(basename(outdir), ".snp.vcf"))
    output.indel.vcf <- file.path(outdir, paste(basename(outdir), ".indel.vcf"))
    
    tryCatch({
        system2("/bin/bash", c("-c", shQuote(varscan_cmd)))
    }, error = function(e) {
        stop("Error: Varscan somatic mutation calling run has failed.\n")
    })

    cat(format(Sys.time(),usetz = TRUE), "\n")
    cat("Done. Finished somatic mutation calling succesfully!!!\n")
    end_time <- Sys.time()

    ### Print the elapsed time for Varscan run
    cat("Time elapsed: ", difftime(end_time, start_time, units = "hours"), "hours\n")
}

varscan_somatic_caller(tumor_bam = opt$'tumor-bam',
    normal_bam = opt$'normal-bam',
    tumor_id = opt$'tumor-id',
    normal_id = opt$'normal-id',
    ref_file = opt$'ref-file',
    bed_file = opt$'bed-file',
    cpus = opt$'cpus',
    outdir = opt$'output-dir',
    samtools_path = opt$'samtools-path',
    varscan_path = opt$'varscan-path')

