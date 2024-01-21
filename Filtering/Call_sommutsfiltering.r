############################################### Call_sommutsfiltering.r
######################## This script executes files filtering somatic variants from 5 different callers: Strelka2, Mutect2,
######################## Varscan, MuSE, Lofreq
################ Inputs:
################ @param anno-file: full path to variant annotation file
################ @param mode: software used to call somatic variants (accepts one argument)
################ @param output-dir: output directory
################ @param total-rd-threshold-snv: total read depth threshold for SNVs
################ @param total-rd-threshold-indel: total read depth threshold for INDELs
################ @param tumor-rd-threshold-snv: tumor read depth threshold for SNVs
################ @param tumor-rd-threshold-indel: tumor read depth threshold for INDELs
################ @param normal-rd-threshold-snv: normal read depth threshold for SNVs
################ @param normal-rd-threshold-indel: normal read depth threshold for INDELs
################ @param mq-threshold-snv: mapping quality threshold for SNVs
################ @param mq-threshold-indel: mapping quality threshold for INDELs
################ @param tumor-afalt-snv: allele frequency threshold of alternate allele in tumor for SNVs
################ @param tumor-afalt-indel: allele frequency threshold of alternate allele in tumor for INDELs
################ @param normal-afalt-snv: allele frequency threshold of alternate allele in normal for SNVs
################ @param normal-afalt-indel: allele frequency threshold of alternate allele in normal for INDELs
################ @param tumor-rdalt-snv: read depth threshold of alternate allele in tumor for SNVs
################ @param tumor-rdalt-indel: read depth threshold of alternate allele in tumor for INDELs
################ @param normal-rdalt-snv: read depth threshold of alternate allele in normal for SNVs
################ @param normal-rdalt-indel: read depth threshold of alternate allele in normal for INDELs
################ @param output-vcf: return vcf file with variants passing filters 

### install required packages
require("optparse")

### set the arguments
option_list = list(
    make_option(c("-s", "--script-dir"), type="character", default=NULL, 
              help="full path to filtering script file", metavar="character"),
    make_option(c("-z", "--anno-file"), type="character", default=NULL, 
              help="full path to variant annotation file", metavar="character"),
    make_option(c("-f", "--mode"), type="character", default=NULL, 
              help="software used to call somatic variants (Options: strelka2, mutect2, varscan, muse, lofreq)", metavar="character"),
    make_option(c("-o", "--output-dir"), type="character", default=NULL, 
              help="output directory", metavar="character"),
    make_option(c("-g", "--total-rd-threshold-snv"), type="numeric", default=40,
              help="total read depth threshold for SNVs", metavar="numeric"),
    make_option(c("-G", "--total-rd-threshold-indel"), type="numeric", default=40,
              help="total read depth threshold for INDELs", metavar="numeric"),
    make_option(c("-a", "--tumor-rd-threshold-snv"), type="numeric", default=25,
              help="tumor read depth threshold for SNVs", metavar="numeric"),
    make_option(c("-A", "--tumor-rd-threshold-indel"), type="numeric", default=25, 
              help="tumor read depth threshold for INDELs", metavar="numeric"),
    make_option(c("-n", "--normal-rd-threshold-snv"), type="numeric", default=15,
              help="normal read depth threshold for SNVs", metavar="numeric"),
    make_option(c("-N", "--normal-rd-threshold-indel"), type="numeric", default=15,
              help="normal read depth threshold for INDELs", metavar="numeric"),
    make_option(c("-m", "--mq-threshold-snv"), type="numeric", default=30,
              help="mapping quality threshold for SNVs", metavar="numeric"),
    make_option(c("-M", "--mq-threshold-indel"), type="numeric", default=50,
              help="mapping quality threshold for INDELs", metavar="numeric"),
    make_option(c("-t", "--tumor-afalt-snv"), type="numeric", default=0.05,
              help="allele frequency threshold of alternate allele in tumor for SNVs", metavar="numeric"),
    make_option(c("-T", "--tumor-afalt-indel"), type="numeric", default=0.05,
              help="allele frequency threshold of alternate allele in tumor for INDELs", metavar="numeric"),
    make_option(c("-l", "--normal-afalt-snv"), type="numeric", default=0.02,
              help="allele frequency threshold of alternate allele in normal for SNVs", metavar="numeric"),
    make_option(c("-L", "--normal-afalt-indel"), type="numeric", default=0.02,
            help="allele frequency threshold of alternate allele in normal for INDELs", metavar="numeric"),
    make_option(c("-p", "--tumor-rdalt-snv"), type="numeric", default=5,
              help="read depth threshold of alternate allele in tumor for SNVs", metavar="numeric"),
    make_option(c("-P", "--tumor-rdalt-indel"), type="numeric", default=5,
            help="read depth threshold of alternate allele in tumor for INDELs", metavar="numeric"),
    make_option(c("-r", "--normal-rdalt-snv"), type="numeric", default=2,
              help="read depth threshold of alternate allele in normal for SNVs", metavar="numeric"),
    make_option(c("-R", "--normal-rdalt-indel"), type="numeric", default=2,
              help="read depth threshold of alternate allele in normal for INDELs", metavar="numeric"),
    make_option(c("-c", "--picard-path"), type="numeric", default=NULL, action="store",
              help="full path to picard executable file", metavar="numeric"),
    make_option(c("-d", "--output-vcf"), type="numeric", default=0, 
              help="return vcf file with variants passing filters", metavar="numeric")); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

### specify job commands
cmd.out <- NULL
cmd.out <- paste0(cmd.out, "#!/bin/bash\n\n")
cmd.out <- paste0(cmd.out, "#SBATCH --time=24:00:00 --mem-per-cpu=16G ", 
    " --partition ag2tb", " --mail-type='END,FAIL' --mail-user balay011@umn.edu -A aventeic --job-name=", file.path(opt$'output-dir', "sommutsfiltering"), 
    " -o ", file.path(opt$'output-dir', "sommutsfiltering.o%J"), " -e ", file.path(opt$'output-dir', "sommutsfiltering.e%J\n\n"))

if (length(opt$'mode') > 1){
    stop("Two arguments were provided but one required\n")
}

if (opt$'mode' == "strelka2"){
    strelka2_cmd <- paste("Rscript", file.path(opt$'script-dir', "FilterStrelka2sommuts.r"), "--anno-file", opt$'anno-file',
        "--output-dir", opt$'output-dir', 
        "--total-rd-threshold-snv", opt$'total-rd-threshold-snv', "--total-rd-threshold-indel", opt$'total-rd-threshold-indel',
        "--tumor-rd-threshold-snv", opt$'tumor-rd-threshold-snv', "--tumor-rd-threshold-indel", opt$'tumor-rd-threshold-indel',
        "--normal-rd-threshold-snv", opt$'normal-rd-threshold-snv', "--normal-rd-threshold-indel", opt$'tumor-rd-threshold-indel',
        "--mq-threshold-snv", opt$'mq-threshold-snv', "--mq-threshold-indel", opt$'mq-threshold-indel', 
        "--tumor-afalt-snv", opt$'tumor-afalt-snv', "--tumor-afalt-indel", opt$'tumor-afalt-indel',
        "--normal-afalt-snv", opt$'normal-afalt-snv', "--normal-afalt-indel", opt$'normal-afalt-indel',
        "--tumor-rdalt-snv", opt$'tumor-rdalt-snv', "--tumor-rdalt-indel", opt$'tumor-rdalt-indel',
        "--normal-rdalt-snv", opt$'normal-rdalt-snv', "--normal-rdalt-indel", opt$'normal-rdalt-indel',
        "--output-vcf", opt$'output-vcf', sep = " ")
    cmd.out <- paste(cmd.out, strelka2_cmd)

    ### write the script to shell file
    cat(cmd.out, file = file.path(opt$'output-dir', "sommutsfiltering_strelka2.sh"))

    ### execute the job using sbatch 
    # system(paste0("sbatch ", file.path(opt$'output-dir', "sommutsfiltering_strelka2.sh")))
}

if (opt$'mode' == "mutect2"){
    mutect2_cmd <- paste("Rscript", file.path(opt$'script-dir', "FilterMutect2sommuts.r"), "--anno-file", opt$'anno-file',
        "--output-dir", opt$'output-dir', "--total-rd-threshold-snv", opt$'total-rd-threshold-snv', "--total-rd-threshold-indel", opt$'total-rd-threshold-indel',
        "--tumor-rd-threshold-snv", opt$'tumor-rd-threshold-snv', "--tumor-rd-threshold-indel", opt$'tumor-rd-threshold-indel',
        "--normal-rd-threshold-snv", opt$'normal-rd-threshold-snv', "--normal-rd-threshold-indel", opt$'tumor-rd-threshold-indel',
        "--mq-threshold-snv", opt$'mq-threshold-snv', "--mq-threshold-indel", opt$'mq-threshold-indel', 
        "--tumor-afalt-snv", opt$'tumor-afalt-snv', "--tumor-afalt-indel", opt$'tumor-afalt-indel',
        "--normal-afalt-snv", opt$'normal-afalt-snv', "--normal-afalt-indel", opt$'normal-afalt-indel',
        "--tumor-rdalt-snv", opt$'tumor-rdalt-snv', "--tumor-rdalt-indel", opt$'tumor-rdalt-indel',
        "--normal-rdalt-snv", opt$'normal-rdalt-snv', "--normal-rdalt-indel", opt$'normal-rdalt-indel',
        "--picard-path", opt$'picard-path', 
        "--output-vcf", opt$'output-vcf', sep = " ")
    cmd.out <- paste(cmd.out, mutect2_cmd)

    ### write the script to shell file
    cat(cmd.out, file = file.path(opt$'output-dir', "sommutsfiltering_mutect2.sh"))

    ### execute the job using sbatch 
    # system(paste0("sbatch ", file.path(opt$'output-dir', "sommutsfiltering_mutect2.sh")))
}

if (opt$'mode' == "muse"){
    muse_cmd <- paste("Rscript", file.path(opt$'script-dir', "FilterMusesommuts.r"), "--anno-file", opt$'anno-file',
        "--output-dir", opt$'output-dir', "--total-rd-threshold-snv", opt$'total-rd-threshold-snv',
        "--tumor-rd-threshold-snv", opt$'tumor-rd-threshold-snv',
        "--normal-rd-threshold-snv", opt$'normal-rd-threshold-snv',
        "--tumor-afalt-snv", opt$'tumor-afalt-snv',
        "--normal-afalt-snv", opt$'normal-afalt-snv',
        "--tumor-rdalt-snv", opt$'tumor-rdalt-snv',
        "--normal-rdalt-snv", opt$'normal-rdalt-snv',
        "--picard-path", opt$'picard-path', 
        "--output-vcf", opt$'output-vcf', sep = " ")
    cmd.out <- paste(cmd.out, muse_cmd)

    ### write the script to shell file
    cat(cmd.out, file = file.path(opt$'output-dir', "sommutsfiltering_muse.sh"))

    ### execute the job using sbatch 
    # system(paste0("sbatch ", file.path(opt$'output-dir', "sommutsfiltering_muse.sh")))
}

if (opt$'mode' == "varscan"){
    varscan_cmd <- paste("Rscript", file.path(opt$'script-dir', "FilterVarscansommuts.r"), "--anno-file", opt$'anno-file',
        "--output-dir", opt$'output-dir', "--total-rd-threshold-snv", opt$'total-rd-threshold-snv', "--total-rd-threshold-indel", opt$'total-rd-threshold-indel',
        "--tumor-rd-threshold-snv", opt$'tumor-rd-threshold-snv', "--tumor-rd-threshold-indel", opt$'tumor-rd-threshold-indel',
        "--normal-rd-threshold-snv", opt$'normal-rd-threshold-snv', "--normal-rd-threshold-indel", opt$'tumor-rd-threshold-indel',
        "--tumor-afalt-snv", opt$'tumor-afalt-snv', "--tumor-afalt-indel", opt$'tumor-afalt-indel',
        "--normal-afalt-snv", opt$'normal-afalt-snv', "--normal-afalt-indel", opt$'normal-afalt-indel',
        "--tumor-rdalt-snv", opt$'tumor-rdalt-snv', "--tumor-rdalt-indel", opt$'tumor-rdalt-indel',
        "--normal-rdalt-snv", opt$'normal-rdalt-snv', "--normal-rdalt-indel", opt$'normal-rdalt-indel', 
        "--output-vcf", opt$'output-vcf', sep = " ")
    cmd.out <- paste(cmd.out, varscan_cmd)

    ### write the script to shell file
    cat(cmd.out, file = file.path(opt$'output-dir', "sommutsfiltering_varscan.sh"))

    ### execute the job using sbatch 
    # system(paste0("sbatch ", file.path(opt$'output-dir', "sommutsfiltering_varscan.sh")))
}

if (opt$'mode' == "lofreq"){
    lofreq_cmd <- paste("Rscript", file.path(opt$'script-dir', "FilterLofreqsommuts.r"), "--anno-file", opt$'anno-file',
        "--output-dir", opt$'output-dir', 
        "--tumor-rd-threshold-snv", opt$'tumor-rd-threshold-snv', "--tumor-rd-threshold-indel", opt$'tumor-rd-threshold-indel',
        "--tumor-afalt-snv", opt$'tumor-afalt-snv', "--tumor-afalt-indel", opt$'tumor-afalt-indel',
        "--tumor-rdalt-snv", opt$'tumor-rdalt-snv', "--tumor-rdalt-indel", opt$'tumor-rdalt-indel',
        "--output-vcf", opt$'output-vcf', sep = " ")
    cmd.out <- paste(cmd.out, lofreq_cmd)

    ### write the script to shell file
    cat(cmd.out, file = file.path(opt$'output-dir', "sommutsfiltering_lofreq.sh"))

    ### execute the job using sbatch 
    # system(paste0("sbatch ", file.path(opt$'output-dir', "sommutsfiltering_lofreq.sh")))
}