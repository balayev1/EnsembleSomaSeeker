############################################### Mutect2FilterSomatic.r 
################ This R script contains a function to assign filters to somatic mutations using several functions from GATK package in tumor-normal pair mode on whole-genome sequencing (WGS) samples
################ Inputs:
################ @param tumor-bam: Full path to tumor BAM file
################ @param normal-bam: Full path to normal BAM file
################ @param: ref-file: Full path to reference genome FASTA file
################ @param biallelic-file: Full path to file with biallelic loci
################ @param mitochondria: whether to filter mitochondrial mutations
################ @param outdir: Output directory
################ @param gatk-path: Full path to GATK executable file
################ Output: see GATK MuTect2 output files (see https://gatk.broadinstitute.org/hc/en-us/articles/360036713111-Mutect2#--intervals) 

### install required packages
require("optparse")

### set the arguments
option_list = list(
    make_option(c("-T", "--tumor-bam"), type="character", default=NULL, 
              help="Full path to tumor BAM file", metavar="character"),
    make_option(c("-N", "--normal-bam"), type="character", default=NULL, 
              help="Full path to normal BAM file", metavar="character"),
    make_option(c("-r", "--ref-file"), type="character", default=NULL, 
              help="Full path to reference genome FASTA file", metavar="character"),
    make_option(c("-b", "--biallelic-file"), type="character", default=NULL, 
              help="Full path to file with biallelic loci", metavar="character"),
    make_option(c("-m", "--mitochondria"), type="logical", default=TRUE, 
              help="whether to filter mitochondrial mutations", metavar="logical"),
    make_option(c("-o", "--output-dir"), type="character", default=NULL, 
              help="Output directory", metavar="character"),
    make_option(c("-g", "--gatk-path"), type="character", default=NULL, 
              help="Full path to GATK executable file", metavar="character")); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

### filter_mutect2somatic function
filter_mutect2somatic <- function(tumor_bam, normal_bam, ref_file, biallelic_file, outdir, gatk_path, mitochondria){
    if (is.null(tumor_bam) || is.null(normal_bam) || is.null(ref_file) || is.null(biallelic_file) || is.null(outdir) || is.null(gatk_path)) {
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
    cat("Beginning to run somatic mutation filtering pipeline in GATK platform\n")

    ### Step 1: Summarize read counts that support reference, alternate and other alleles for given sites 
    ### Set GetPileupSummaries command for both tumor and matched normal sample
    pileupsum_tumor_cmd <- paste(gatk_path, "GetPileupSummaries", "-I", tumor_bam,
        "-L", biallelic_file, "-V", biallelic_file, 
        "-O", file.path(outdir, paste0(basename(outdir), "_tumor.pileups.table")), sep = " ")

    pileupsum_normal_cmd <- paste(gatk_path, "GetPileupSummaries", "-I", normal_bam,
        "-L", biallelic_file, "-V", biallelic_file, 
        "-O", file.path(outdir, paste0(basename(outdir), "_normal.pileups.table")), sep = " ")

    ### Execute the GetPileupSummaries command
    tryCatch({
        system(pileupsum_tumor_cmd)
        system(pileupsum_normal_cmd)
    }, error = function(e) {
        stop("Error: GetPileupSummaries run has failed.\n")
    })

    cat("Finished GetPileupSummaries run succesfully\n")
    end_time <- Sys.time()

    ### Print the elapsed time for GetPileupSummaries run
    cat("Time elapsed for GetPileupSummaries run: ", difftime(end_time, start_time, units = "hours"), "hours\n")

    ### Step 2: Calculate the fraction of reads coming from cross-sample contamination
    start_time <- Sys.time()
    ### Set CalculateContamination command
    calccont_cmd <- paste(gatk_path, "CalculateContamination", "-I", file.path(outdir, paste0(basename(outdir), "_tumor.pileups.table")),
        "--matched-normal", file.path(outdir, paste0(basename(outdir), "_normal.pileups.table")),
        "-O", file.path(outdir, paste0(basename(outdir), "_contamination.table")),
        "--tumor-segmentation", file.path(outdir, paste0(basename(outdir), "_segments.table")), sep = " ")

    ### Execute the CalculateContamination command
    tryCatch({
        system(calccont_cmd)
    }, error = function(e) {
        stop("Error: CalculateContamination run has failed.\n")
    })

    cat("Finished CalculateContamination run succesfully\n")
    end_time <- Sys.time()

    ### Print the elapsed time for CalculateContamination run
    cat("Time elapsed for CalculateContamination run: ", difftime(end_time, start_time, units = "hours"), "hours\n")

    ### Step 3: Learn probability of read orientation artifacts 
    start_time <- Sys.time()
    ### Set LearnReadOrientationModel command
    learn_readori_cmd <- paste(gatk_path, "LearnReadOrientationModel", 
        "-I", file.path(outdir, paste0(basename(outdir), ".f1r2.tar.gz")), 
        "-O", file.path(outdir, paste0(basename(outdir), ".artifact_prior_tables.tar.gz")), sep = " ")

    ### Execute the LearnReadOrientationModel command
    tryCatch({
        system(learn_readori_cmd)
    }, error = function(e) {
        stop("Error: LearnReadOrientationModel run has failed.\n")
    })

    cat("Finished LearnReadOrientationModel run succesfully\n")
    end_time <- Sys.time()

    ### Print the elapsed time for FilterMutectCalls run
    cat("Time elapsed for LearnReadOrientationModel: ", difftime(end_time, start_time, units = "hours"), "hours\n")


    ### Step 4: Filter Mutect2 calls
    start_time <- Sys.time()
    ### Set FilterMutectCalls command
    filtersom_cmd <- paste(gatk_path, "FilterMutectCalls", "-R", ref_file, 
        "-V", file.path(outdir, paste0(basename(outdir), "_unfiltered.vcf.gz")), 
        "--orientation-bias-artifact-priors", file.path(outdir, paste0(basename(outdir), ".artifact_prior_tables.tar.gz")),
        "--contamination-table", file.path(outdir, paste0(basename(outdir), "_contamination.table")),
        "--tumor-segmentation", file.path(outdir, paste0(basename(outdir), "_segments.table")),
        "-O", file.path(outdir, paste0(basename(outdir), "_filtered.vcf.gz")), sep = " ")

    ### Execute the FilterMutectCalls command
    tryCatch({
        system(filtersom_cmd)
    }, error = function(e) {
        stop("Error: FilterMutectCalls run has failed.\n")
    })

    cat("Finished FilterMutectCalls run succesfully\n")
    end_time <- Sys.time()

    ### Print the elapsed time for FilterMutectCalls run
    cat("Time elapsed for FilterMutectCalls: ", difftime(end_time, start_time, units = "hours"), "hours\n")

    if (mitochondria == TRUE){
        cat("Mitochondrial run was requested ...\n")
        start_time <- Sys.time()
        ### Set LearnReadOrientationModel command
        learn_readori_chrM_cmd <- paste(gatk_path, "LearnReadOrientationModel", 
            "-I", file.path(outdir, paste0(basename(outdir), "_chrM.f1r2.tar.gz")), 
            "-O", file.path(outdir, paste0(basename(outdir), "_chrM.artifact_prior_tables.tar.gz")), sep = " ")
        ### Set FilterMutectCalls command
        filtersom_chrM_cmd <- paste(gatk_path, "FilterMutectCalls", "-R", ref_file, 
            "-V", file.path(outdir, paste0(basename(outdir), "_chrM_unfiltered.vcf.gz")), 
            "--mitochondria-mode true", 
            "--orientation-bias-artifact-priors", file.path(outdir, paste0(basename(outdir), "_chrM.artifact_prior_tables.tar.gz")),
            "--contamination-table", file.path(outdir, paste0(basename(outdir), "_contamination.table")),
            "--tumor-segmentation", file.path(outdir, paste0(basename(outdir), "_segments.table")),
            "-O", file.path(outdir, paste0(basename(outdir), "_chrM_filtered.vcf.gz")), sep = " ")
        tryCatch({
        system(learn_readori_chrM_cmd)
        system(filtersom_chrM_cmd)
        }, error = function(e) {
            stop("Error: Filtering run of mitochondrial mutations has failed.\n")
        })

        cat("Finished filtering run of mitochondrial mutations succesfully\n")
        end_time <- Sys.time()
        ### Print the elapsed time for FilterMutectCalls mitochondria run
        cat("Time elapsed for filtering run of mitochondrial mutations: ", difftime(end_time, start_time, units = "hours"), "hours\n")

        ### Gather vcf files from specified chromosomes and chrM together
        start_time <- Sys.time()
        gather_vcf_cmd <- paste(gatk_path, "GatherVcfs", "-I", file.path(outdir, paste0(basename(outdir), "_filtered.vcf.gz")),
            "-I", file.path(outdir, paste0(basename(outdir), "_chrM_filtered.vcf.gz")),
            "-O", file.path(outdir, paste0(basename(outdir), "_allfiltered.vcf.gz")))
        tryCatch({
        system(gather_vcf_cmd)
        }, error = function(e) {
            stop("Error: Gathering vcf and stats files has failed.\n")
        })

        cat("Finished gathering vcf and stats files succesfully\n")
        end_time <- Sys.time()
        ### Print the elapsed time for gathering vcf and stats files
        cat("Time elapsed for gathering vcf and stats files: ", difftime(end_time, start_time, units = "hours"), "hours\n")

        system(paste("mv", file.path(outdir, paste0(basename(outdir), "_allfiltered.vcf.gz")), 
            file.path(outdir, paste0(basename(outdir), "_filtered.vcf.gz")), sep = " "))
        
        ### Remove all intermediate files
        system(paste("rm", file.path(outdir, paste0(basename(outdir), "_filtered.vcf.gz.filteringStats.tsv")), 
            file.path(outdir, paste0(basename(outdir), "_chrM_filtered.vcf.gz*")), sep = " "))
    }
}

filter_mutect2somatic(tumor_bam = opt$'tumor-bam',
    normal_bam = opt$'normal-bam',
    ref_file = opt$'ref-file',
    biallelic_file = opt$'biallelic-file',
    mitochondria = opt$'mitochondria',
    outdir = opt$'output-dir',
    gatk_path = opt$'gatk-path')