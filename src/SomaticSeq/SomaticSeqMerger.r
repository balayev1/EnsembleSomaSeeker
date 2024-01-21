############################################### SomaticSeqMerger.r 
################ This R script contains a function to merge somatic mutation from different callers using SomaticSeq package in tumor-normal pair mode on whole-genome sequencing (WGS) samples
################ Inputs:
################ @param tumor-bam: Full path to tumor BAM file
################ @param normal-bam: Full path to normal BAM file
################ @param mutect2-vcf: Full path to Mutect2 VCF file
################ @param varscan-snv: Full path to Varscan snv VCF file
################ @param varscan-indel: Full path to Varscan indel VCF file
################ @param jsm-vcf: Full path to JSM VCF file
################ @param somaticsniper-vcf: Full path to Somaticsniper VCF file
################ @param vardict-vcf: Full path to Vardict VCF file
################ @param muse-vcf: Full path to MuSE VCF file
################ @param lofreq-snv: Full path to Lofreq snv VCF file
################ @param lofreq-indel: Full path to Lofreq indel VCF file
################ @param scalpel-vcf: Full path to Scalpel VCF file
################ @param strelka-snv: Full path to Strelka2 snv VCF file
################ @param strelka-indel: Full path to Strelka2 indel VCF file
################ @param arbitrary-snvs: Full path to arbitrary snv VCF file(s)
################ @param arbitrary-indels: Full path to arbitrary indel VCF file(s)
################ @param output-directory: Output directory
################ @param genome-reference: Full path to reference genome FASTA file
################ @param inclusion-region: Full path to bed file with included genomic coordinates
################ @param exclusion-region: Full path to bed file with excluded genomic coordinates
################ @param pass-threshold: SCORE for PASS (default: 0.5)
################ @param lowqual-threshold: SCORE for PASS (default: 0.1)
################ @param dbsnp-vcf: Full path to dbSNP VCF file
################ @param cosmic-vcf: Full path to COSMIC VCF file
################ @param somaticseq-path: Full path to SomaticSeq executable file
################ @param threads: Number of threads
################ Outputs: see SomaticSeq output files (see https://github.com/bioinform/somaticseq) 


### install required packages
suppressPackageStartupMessages(require(optparse))

### set the arguments
option_list = list(
    make_option(c("-T", "--tumor-bam"), type="character", 
              help="Full path to tumor BAM file", metavar="character"),
    make_option(c("-N", "--normal-bam"), type="character", 
              help="Full path to normal BAM file", metavar="character"),
    make_option(c("-m", "--mutect2-vcf"), type="character", default=NA, 
              help="Full path to Mutect2 VCF file", metavar="character"),
    make_option(c("-V", "--varscan-snv"), type="character", default=NA, 
              help="Full path to Varscan snv VCF file", metavar="character"),
    make_option(c("-v", "--varscan-indel"), type="character", default=NA, 
              help="Full path to Varscan indel VCF file", metavar="character"),
    make_option(c("-j", "--jsm-vcf"), type="character", default=NA, 
        help="Full path to JSM VCF file", metavar="character"),
    make_option(c("-s", "--somaticsniper-vcf"), type="character", default=NA, 
              help="Full path to Somaticsniper VCF file", metavar="character"),
    make_option(c("-c", "--vardict-vcf"), type="character", default=NA, 
              help="Full path to Vardict VCF file", metavar="character"),
    make_option(c("-U", "--muse-vcf"), type="character", default=NA, 
              help="Full path to MuSE VCF file", metavar="character"),
    make_option(c("-L", "--lofreq-snv"), type="character", default=NA, 
              help="Full path to Lofreq snv VCF file", metavar="character"),
    make_option(c("-l", "--lofreq-indel"), type="character", default=NA, 
              help="Full path to Lofreq indel VCF file", metavar="character"),
    make_option(c("-S", "--scalpel-vcf"), type="character", default=NA, 
              help="Full path to Scalpel VCF file", metavar="character"),
    make_option(c("-R", "--strelka-snv"), type="character", default=NA, 
              help="Full path to Strelka2 snv VCF file", metavar="character"),
    make_option(c("-r", "--strelka-indel"), type="character", default=NA, 
              help="Full path to Strelka2 indel VCF file", metavar="character"),
    make_option(c("-A", "--arbitrary-snvs"), type="character", default=NA, 
              help="Full path to arbitrary snv VCF file(s)", metavar="character"),
    make_option(c("-a", "--arbitrary-indels"), type="character", default=NA, 
              help="Full path to arbitrary indel VCF file(s)", metavar="character"),
    make_option(c("-o", "--output-directory"), type="character", 
              help="Output directory", metavar="character"),
    make_option(c("-g", "--genome-reference"), type="character", 
              help="Full path to reference genome FASTA file", metavar="character"),
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
    make_option(c("-O", "--somaticseq-path"), type="character", 
              help="Full path to SomaticSeq executable file", metavar="character"),
    make_option(c("-t", "--threads"), type="numeric", default=1, 
              help="Number of threads", metavar="numeric"));

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)
print(opt)

### varscan_somatic_caller function
somaticseq_merger <- function(tumor_bam, normal_bam, output_directory, genome_reference, somaticseq_path, inclusion_region, 
    exclusion_region, threads=1, pass_threshold=0.5, lowqual_threshold=0.1, dbsnp, cosmic, 
    mutect2_vcf, varscan_snv, varscan_indel, jsm_vcf, 
    somaticsniper_vcf, vardict_vcf, muse_vcf, lofreq_snv, lofreq_indel, scalpel_vcf, 
    strelka_snv, strelka_indel, arbitrary_snvs, arbitrary_indels){

    if (missing(tumor_bam) || missing(normal_bam) || missing(output_directory) || missing(genome_reference) || missing(somaticseq_path)) {
        stop("Missing required arguments. Please provide all required options.")
    }

    ### Check if output directory exists
    if (!file.exists(output_directory)) {
        cat("Output directory does not exist... Creating directory\n")
        dir.create(output_directory)
    }

    cat("============================================================\n")
    cat(format(Sys.time(),usetz = TRUE), "\n")
    start_time <- Sys.time()
    cat("Processing sample pair:", tumor_bam, " ", normal_bam, "\n")
    cat("Beginning to run SomaticSeq\n")

    # Set SomaticSeq main options
    somaticseq.main.options <- NULL
    somaticseq.main.options <- paste(somaticseq.main.options, "--output-directory", output_directory, 
        "--genome-reference", genome_reference, "--threads", threads, "--pass-threshold", pass_threshold,
        "--lowqual-threshold", lowqual_threshold, sep = " ")

    # add included or excluded regions if provided
    if (!is.na(inclusion_region) & file.exists(inclusion_region)){
        somaticseq.main.options <- paste(somaticseq.main.options, "--inclusion-region", inclusion_region, sep = " ")
    }

    if (!is.na(exclusion_region) & file.exists(exclusion_region)){
        somaticseq.main.options <- paste(somaticseq.main.options, "--exclusion-region", exclusion_region, sep = " ")
    }

    # add dbSNP and COSMIC files if provided
    if (!is.na(dbsnp) & file.exists(dbsnp)){
        somaticseq.main.options <- paste(somaticseq.main.options, "--dbsnp-vcf", dbsnp, sep = " ")
    }

    if (!is.na(cosmic) & file.exists(cosmic)){
        somaticseq.main.options <- paste(somaticseq.main.options, "--cosmic-vcf", cosmic, sep = " ")
    }

    somaticseq.main.options <- paste(somaticseq.main.options, "paired", "--tumor-bam-file", tumor_bam, "--normal-bam-file", normal_bam, sep = " ")

    # Set SomaticSeq secondary options
    somaticseq.options <- NULL

    # add VCF files from mutation caller if provided
    if (!is.na(mutect2_vcf) & file.exists(mutect2_vcf)){
        somaticseq.options <- paste(somaticseq.options, "--mutect2-vcf", mutect2_vcf, sep = " ")
    }

    if (!is.na(varscan_snv) & file.exists(varscan_snv)){
        somaticseq.options <- paste(somaticseq.options, "--varscan-snv", varscan_snv, sep = " ")
    }

    if (!is.na(varscan_indel) & file.exists(varscan_indel)){
        somaticseq.options <- paste(somaticseq.options, "--varscan-indel", varscan_indel, sep = " ")
    }

    if (!is.na(jsm_vcf) & file.exists(jsm_vcf)){
        somaticseq.options <- paste(somaticseq.options, "--jsm-vcf", jsm_vcf, sep = " ")
    }

    if (!is.na(somaticsniper_vcf) & file.exists(somaticsniper_vcf)){
        somaticseq.options <- paste(somaticseq.options, "--somaticsniper-vcf", somaticsniper_vcf, sep = " ")
    }

    if (!is.na(vardict_vcf) & file.exists(vardict_vcf)){
        somaticseq.options <- paste(somaticseq.options, "--vardict-vcf", vardict_vcf, sep = " ")
    }

    if (!is.na(muse_vcf) & file.exists(muse_vcf)){
        somaticseq.options <- paste(somaticseq.options, "--muse-vcf", muse_vcf, sep = " ")
    }

    if (!is.na(lofreq_snv) & file.exists(lofreq_snv)){
        somaticseq.options <- paste(somaticseq.options, "--lofreq-snv", lofreq_snv, sep = " ")
    }

    if (!is.na(lofreq_indel) & file.exists(lofreq_indel)){
        somaticseq.options <- paste(somaticseq.options, "--lofreq-indel", lofreq_indel, sep = " ")
    }

    if (!is.na(scalpel_vcf) & file.exists(scalpel_vcf)){
        somaticseq.options <- paste(somaticseq.options, "--scalpel-vcf", scalpel_vcf, sep = " ")
    }

    if (!is.na(strelka_snv) & file.exists(strelka_snv)){
        somaticseq.options <- paste(somaticseq.options, "--strelka-snv", strelka_snv, sep = " ")
    }
    
    if (!is.na(strelka_indel) & file.exists(strelka_indel)){
        somaticseq.options <- paste(somaticseq.options, "--strelka-indel", strelka_indel, sep = " ")
    }

    if (!is.na(arbitrary_snvs)){
        arb_snv_files <- unlist(strsplit(arbitrary_snvs, ","))
        arb.options <- NULL
        num <- 0
        for (i in arb_snv_files){
            if (file.exists(i)){
                arb.options <- paste(arb.options, i, sep = " ")
                num <- num + 1
            }
        }
        if (num > 0){
            somaticseq.options <- paste(somaticseq.options, "--arbitrary-snvs", arb.options, sep = " ")
        }
    }

    if (!is.na(arbitrary_indels)){
        arb_indel_files <- unlist(strsplit(arbitrary_indels, ","))
        arb.options <- NULL
        num <- 0
        for (j in arb_indel_files){
            if (file.exists(j)){
                arb.options <- paste(arb.options, j, sep = " ")
                num <- num + 1
            }
        }
        if (num > 0){
            somaticseq.options <- paste(somaticseq.options, "--arbitrary-indels", arb.options, sep = " ")
        }
    }

    ### Set SomaticSeq command
    somaticseq_cmd <- paste(somaticseq_path, somaticseq.main.options, somaticseq.options, sep = " ")
    print(somaticseq_cmd)
    ### Execute the SomaticSeq command
    tryCatch({
        system(somaticseq_cmd)
    }, error = function(e) {
        stop("Error: SomaticSeq run has failed.\n")
    })

    cat("Finished SomaticSeq run succesfully\n")
    end_time <- Sys.time()

    ### Print the elapsed time for SomaticSeq run
    cat("Time elapsed: ", difftime(end_time, start_time, units = "hours"), "hours\n")
}

### Run somaticseq_merger function
somaticseq_merger(tumor_bam = opt$'tumor-bam',
  normal_bam = opt$'normal-bam',
  output_directory = opt$'output-directory',
  genome_reference = opt$'genome-reference',
  somaticseq_path = opt$'somaticseq-path',
  inclusion_region = opt$'inclusion-region',
  exclusion_region = opt$'exclusion-region',
  pass_threshold = opt$'pass-threshold',
  lowqual_threshold = opt$'lowqual-threshold',
  dbsnp = opt$'dbsnp-vcf',
  cosmic = opt$'cosmic-vcf',
  threads = opt$'threads',
  mutect2_vcf = opt$'mutect2-vcf',
  varscan_snv = opt$'varscan-snv',
  varscan_indel = opt$'varscan-indel',
  jsm_vcf = opt$'jsm-vcf',
  somaticsniper_vcf = opt$'somaticsniper-vcf',
  vardict_vcf = opt$'vardict-vcf',
  muse_vcf = opt$'muse-vcf',
  lofreq_snv = opt$'lofreq-snv',
  lofreq_indel = opt$'lofreq-indel',
  scalpel_vcf = opt$'scalpel-vcf',
  strelka_snv = opt$'strelka-snv',
  strelka_indel = opt$'strelka-indel',
  arbitrary_snvs = opt$'arbitrary-snvs',
  arbitrary_indels = opt$'arbitrary-indels')


