############################################### Mutect2SomaticCaller.r 
################ This R script contains a function to call somatic mutations using Mutect2 function from GATK package in tumor-normal pair mode on whole-genome sequencing (WGS) samples
################ Inputs:
################ @param tumor-bam: Full path to tumor BAM file
################ @param normal-bam: Full path to normal BAM file
################ @param tumor-id: Tumor sample ID
################ @param normal-id: Normal sample ID
################ @param: ref-file: Full path to reference genome FASTA file
################ @param: bed-file: Full path to bed file with genomic coordinates of interest
################ @param: germ-file: Full path to germline resource file
################ @param pon: Full path to panel of normals file
################ @param cpus: number of cpus
################ @param mitochondria: whether to call mitochondrial mutations
################ @param tmp-dir: temporary directory
################ @param outdir: Output directory
################ @param gatk-path: Full path to GATK executable file
################ Output: see GATK MuTect2 output files (see https://gatk.broadinstitute.org/hc/en-us/articles/360036713111-Mutect2#--intervals) 
################ ! NOTE: chromosomes should have 'chr' prefix in reference_genome.fa and in bed file if provided

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
    make_option(c("-e", "--bed-file"), type="character", default=NULL, 
        help="Full path to bed file with genomic coordinates of interest", metavar="character"),
    make_option(c("-s", "--germ-file"), type="character", default=NULL, 
              help="Full path to germline resource file", metavar="character"),
    make_option(c("-d", "--pon"), type="character", default=NULL, 
              help="Full path to panel of normals file", metavar="character"),
    make_option(c("-z", "--cpus"), type="integer", default=1, 
              help="number of cpus", metavar="integer"),
    make_option(c("-m", "--mitochondria"), type="logical", default=TRUE, 
              help="whether to call mitochondrial mutations", metavar="logical"),
    make_option(c("-p", "--tmp-dir"), type="character", default=NULL, 
              help="temporary directory", metavar="character"),
    make_option(c("-o", "--output-dir"), type="character", default=NULL, 
              help="Output directory", metavar="character"),
    make_option(c("-g", "--gatk-path"), type="character", default=NULL, 
              help="Full path to GATK executable file", metavar="character")); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

### mutect2_somatic_caller function
mutect2_somatic_caller <- function(tumor_bam, normal_bam, tumor_id, normal_id, ref_file, bed_file, germ_file, pon, cpus, outdir, gatk_path, mitochondria, tmpdir){
    if (is.null(tumor_bam) || is.null(normal_bam) || is.null(tumor_id) || is.null(normal_id) || is.null(ref_file) || is.null(outdir) || is.null(gatk_path)) {
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
    cat("Beginning to run somatic mutation calling using Mutect2\n")
    
    ### Set Mutect2 options
    mutect2.options <- paste("--native-pair-hmm-threads", cpus, "--af-of-alleles-not-in-resource 0.0000025",
        "--reference", ref_file, sep = " ")

    ### add germline file if included
    if (!is.null(germ_file)) {
        mutect2.options <- paste(mutect2.options, "--germline-resource", germ_file, sep = " ")
    }

    ### add panel of normals file if included
    if (!is.null(pon)) {
        mutect2.options <- paste(mutect2.options, "--panel-of-normals", pon, sep = " ")
    }

    ### add tmpdir if included
    if (!is.null(tmpdir)){
        mutect2.options <- paste(mutect2.options, "--tmp-dir", tmpdir, sep = " ")
    }

    ### Set Mutect2 command to process chromosome M
    if (mitochondria == TRUE){
        mutect2.chrM.options <- mutect2.options ### copy mutect2 options if mitochondria == TRUE
        mutect2_chrm_cmd <- paste(gatk_path, "Mutect2", "-I", tumor_bam, "-tumor", tumor_id, "-I", normal_bam, "-normal", normal_id,
            mutect2.chrM.options, "-O", file.path(outdir, paste0(basename(outdir), "_chrM_unfiltered.vcf.gz")), "--intervals chrM --mitochondria-mode true",
            "--f1r2-tar-gz", file.path(outdir, paste0(basename(outdir), "_chrM.f1r2.tar.gz")), sep = " ")
    }

    ### add bed file if included
    if (!is.null(bed_file)) {
        temp_bed_file <- paste(tumor_id, normal_id, "chrs1_22XY.bed", sep = "_")
        system(paste("grep -v -w 'chrM'", bed_file, ">", temp_bed_file, sep = " "))
        mutect2.options <- paste(mutect2.options, "--intervals", temp_bed_file, "--interval-padding 0", sep = " ")
    } else {
        mutect2.options <- paste(mutect2.options, "--exclude-intervals chrM", sep = " ")
    }

    ### Set Mutect2 command to process all given chromosomes excluding chrM
    mutect2_cmd <- paste(gatk_path, "Mutect2 -I", tumor_bam, "-tumor", tumor_id, "-I", normal_bam, "-normal", normal_id,
        mutect2.options, "-O", file.path(outdir, paste0(basename(outdir), "_unfiltered.vcf.gz")),  
        "--f1r2-tar-gz", file.path(outdir, paste0(basename(outdir), ".f1r2.tar.gz")), sep = " ")
        
    ### Execute the commands
    tryCatch({
        system(mutect2_cmd)
        if (mitochondria == TRUE){
            system(mutect2_chrm_cmd)
        }
    }, error = function(e) {
        stop("Error: Mutect2 run has failed.\n")
    })

    cat("Finished somatic mutation calling run using Mutect2 succesfully\n")
    end_time <- Sys.time()

    ### Print the elapsed time for Mutect2 run
    cat("Time elapsed: ", difftime(end_time, start_time, units = "hours"), "hours\n")

    ### Remove intermediate files
    system(paste("rm ", temp_bed_file))
}

mutect2_somatic_caller(tumor_bam = opt$'tumor-bam',
    normal_bam = opt$'normal-bam',
    tumor_id = opt$'tumor-id',
    normal_id = opt$'normal-id',
    ref_file = opt$'ref-file',
    bed_file = opt$'bed-file',
    germ_file = opt$'germ-file',
    pon = opt$'pon',
    cpus = opt$'cpus',
    mitochondria = opt$'mitochondria',
    tmpdir = opt$'tmp-dir',
    outdir = opt$'output-dir',
    gatk_path = opt$'gatk-path')

