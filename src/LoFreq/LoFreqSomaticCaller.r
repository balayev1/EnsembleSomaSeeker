############################################### LoFreqSomaticCaller.r 
################ This R script contains a function to call somatic mutations using somatic function in LoFreq package in tumor-normal pair mode 
################ on whole-genome sequencing (WGS) and whole-exome sequencing (WES) samples
################ Inputs:
################ @param tumor: full path to tumor BAM file
################ @param normal: full path to normal BAM file
################ @param ref: full path to reference genome GRCh38 FASTA file
################ @param bed_file: full path to bed file with genome coordinates of interest which contains chromosome name, start, end (optional)
################ @param outdir: output directory
################ @param cpus: number of cpus
################ @param dbsnp: full path to dbsnp file
################ @param lofreq_path: full path to LoFreq executable file
################ @param min-cov: minimum coverage for somatic calls
################ @param samtools-path: full path to Samtools executable file
################ Outputs: see 'lofreq somatic' output files (see https://csb5.github.io/lofreq/commands/) 


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
    make_option(c("-o", "--output-dir"), type="character", default=NULL, 
              help="output directory", metavar="character"),
    make_option(c("-c", "--cpus"), type="numeric", default=1, 
              help="number of cpus", metavar="numeric"),
    make_option(c("-d", "--dbsnp"), type="character", default=NA, 
              help="full path to dbsnp file", metavar="character"),
    make_option(c("-p", "--lofreq-path"), type="character", default=NULL, 
              help="full path to LoFreq executable file", metavar="character"),
    make_option(c("-m", "--min-cov"), type="numeric", default=7, 
              help="minimum coverage for somatic calls", metavar="numeric"),
    make_option(c("-s", "--samtools-path"), type="character", default=NULL, 
              help="full path to samtools executable file", metavar="character")); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

### lofreq_somatic_caller function
lofreq_somatic_caller <- function(tumor, normal, ref, bed_file, outdir, dbsnp, cpus, lofreq_path, min_cov, samtools_path){
    if (is.null(tumor) || is.null(normal) || is.null(ref) || is.null(outdir) || is.null(lofreq_path) || is.null(samtools_path)) {
        stop("Missing required arguments. Please provide all required options.")
    }

    cat("============================================================\n")
    cat("Processing sample pair:", tumor, " ", normal, "\n")
    cat(format(Sys.time(),usetz = TRUE), "\n")
    cat("Running indel realignment using LoFreq\n")

    ### Due to extra indel realingment step, new tumor and normal BAM files will be generated, setting paths
    new_normal_path <- file.path(outdir, paste0(strsplit(basename(normal), split=".bam")[[1]][1], "_BD_BI.bam"))
    new_tumor_path <- file.path(outdir, paste0(strsplit(basename(tumor), split=".bam")[[1]][1], "_BD_BI.bam"))

    ### Write the script for LoFreq somatic mutation caller execution
    tryCatch({
        system2(lofreq_path, args = c("indelqual", "--dindel", "--ref", ref, "--out", new_normal_path, normal))
        system2(lofreq_path, args = c("indelqual", "--dindel", "--ref", ref, "--out", new_tumor_path, tumor))
        system2(samtools_path, args = c("index", "-o", paste0(strsplit(new_normal_path, split=".bam")[[1]][1], ".bai"), new_normal_path))
        system2(samtools_path, args = c("index", "-o", paste0(strsplit(new_tumor_path, split=".bam")[[1]][1], ".bai"), new_tumor_path))
    }, error = function(e) {
        stop("Error executing commands for indel realignment and indexing.")
    })

    cat("Done. Finished running indel realignment using LoFreq\n")
    cat("Running somatic mutation calling using LoFreq\n")

    outdir <- paste0(outdir, "/")

    ### Add lofreq options
    lofreq.options <- paste("--outprefix", outdir, "--ref", ref, "--threads", cpus, 
        "--call-indels", "--min-cov", min_cov, sep = " ")

    ### Add bed file if specified
    if (!is.na(bed_file)){
        lofreq.options <- paste(lofreq.options, "--bed", bed_file, sep = " ")
    }

    ### Add dbsnp file with germline variants
    if (!is.na(dbsnp)){
        lofreq.options <- paste(lofreq.options, "--dbsnp", dbsnp, sep = " ")
    }

    ### Specify lofreq command
    lofreq_cmd <- paste(lofreq_path, "somatic", "--normal", new_normal_path, "--tumor", new_tumor_path, lofreq.options)

    ### Execute lofreq command
    tryCatch({
        system(lofreq_cmd)

        snv.output.vcf <- file.path(outdir, "somatic_final.snvs.vcf.gz")
        indel.output.vcf <- file.path(outdir, "somatic_final.indels.vcf.gz")
        
        if (file.exists(snv.output.vcf) && file.exists(indel.output.vcf) && file.size(snv.output.vcf) != 0L && file.size(indel.output.vcf) != 0L) {
            cat(format(Sys.time(), usetz = TRUE), "\n")
            cat("Done. Finished somatic mutation calling successfully!!!")
        } else {
            stop("Somatic mutation calling failed or produced empty output files.")
        }
    }, error = function(e) {
        stop("Error executing somatic mutation calling commands.")
    })
}

lofreq_somatic_caller(tumor=opt$'tumor-bam', 
    normal=opt$'normal-bam', 
    ref=opt$'ref-file', 
    bed_file=opt$'bed-file', 
    outdir=opt$'output-dir', 
    cpus=opt$'cpus', 
    dbsnp=opt$'dbsnp',
    lofreq_path=opt$'lofreq-path',
    min_cov=opt$'min-cov',
    samtools_path=opt$'samtools-path')
