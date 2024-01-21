############################################### FilterMutect2sommuts.r
######################## This script performs graphical analysis of Mutect2 output vcf files from matched tumor-normal analysis
######################## and returns diagnostic plots and filtered variants in vcf format
################ Inputs:
################ @param anno-file: full path to variant annotation file
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
################ @param picard-path: full path to picard executable file
################ @param output-vcf: return vcf file with variants passing filters 
################ Outputs: diagnostic plots and filtered variants in vcf format (if requested)

### install required packages
require("optparse")

### set the arguments
option_list = list(
    make_option(c("-z", "--anno-file"), type="character", default=NULL, 
              help="full path to variant annotation file", metavar="character"),
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
    make_option(c("-c", "--picard-path"), type="character", default=NULL,
              help="full path to picard executable file", metavar="character"),
    make_option(c("-d", "--output-vcf"), type="numeric", default=0, 
              help="return vcf file with variants passing filters", metavar="numeric")); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

### load required libraries
require(data.table)
require(ggplot2)

### define output directory  
output.dir <- opt$'output-dir'
if (!dir.exists(output.dir)) {
    cat("Creating output directory\n")
    dir.create(output.dir)
}

### make folder to output filtering vcfs with SNVs
output.dir <- file.path(output.dir, "Filtered_mutect2vcfs")
if (!dir.exists(output.dir)){
    dir.create(output.dir)
}

### make subfolder for SNV vcfs
output.snvdir <- file.path(output.dir, "snv")
if (!dir.exists(output.snvdir)){
    dir.create(output.snvdir)
}

### make subfolder for INDEL vcfs
output.indeldir <- file.path(output.dir, "indel")
if (!dir.exists(output.indeldir)){
    dir.create(output.indeldir)
}

### if to output vcfs passing somatic filters
if (opt$'output-vcf' == 1){
    vcf.output.dir <- file.path(output.dir, "somatic_vcfs")
    if (!file.exists(vcf.output.dir)){
        dir.create(vcf.output.dir)
    }
}

### load annotation file
ANNO_FILE <- read.table(opt$'anno-file', sep="\t", header=TRUE)
ANNO_FILE <- cbind(ANNO_FILE[, c(1,2,3)],  ANNO_FILE$mutect2_vcf)
colnames(ANNO_FILE) <- c("subject_id", "tumor_id", "normal_id", "vcf")

### replace blank regions with NA
ANNO_FILE[ANNO_FILE == ""] <- NA

### list paths to sample-specific (i.e. tumor-normal pair) vcf files
vcf.paths <- c()
for (k in 1:nrow(ANNO_FILE)){
    if (!is.na(ANNO_FILE$vcf[k]) & file.exists(ANNO_FILE$vcf[k])){
        vcf.paths <- append(vcf.paths, ANNO_FILE$vcf[k])
    }
}
cat("Found total", length(vcf.paths), "existing vcf files\n")


### function that splits each vcf file into SNVs and INDELs
splitvcfs <- function(vcf.path, picard.path){
    system(paste(picard.path, " SplitVcfs I=", vcf.path, " SNP_OUTPUT=", file.path(dirname(vcf.path), "snvs.vcf.gz"),
        " INDEL_OUTPUT=", file.path(dirname(vcf.path), "indels.vcf.gz"), " CREATE_INDEX=true", " STRICT=false ",
        " VERBOSITY=ERROR QUIET=true"))
}

### split vcf file into SNVs and INDELs
vcf.snvpaths <- c()
vcf.indelpaths <- c()
cat("Splitting vcf files into SNVs and INDELs\n")
for (f in 1:length(vcf.paths)){
    splitvcfs(vcf.path=vcf.paths[f], picard.path=opt$'picard-path')
    vcf.snvpaths <- append(vcf.snvpaths, file.path(dirname(vcf.paths[f]), "snvs.vcf.gz"))
    vcf.indelpaths <- append(vcf.indelpaths, file.path(dirname(vcf.paths[f]), "indels.vcf.gz"))
}

####################### SNV filtering
#### set working directory
setwd(output.snvdir)

snv.list <- list()

start_time <- Sys.time()
for (i in 1:length(vcf.snvpaths)){
    cat("Processing", vcf.snvpaths[i], "\n")

    tumor.id <- ANNO_FILE$tumor_id[i]
    normal.id <- ANNO_FILE$normal_id[i]

    ### read vcf file
    vcf <- fread(vcf.snvpaths[i], data.table=FALSE, skip="#CHROM")

    ### Combine all SNVs into single dataframe with following columns:
    ### chromosome, position, reference base, alternative base, total read depth in normal sample excl filtered reads,
    ### total read depth in tumor sample excl filtered reads, read depth of reference base in normal sample
    ### read depth of reference base in tumor sample, read depth of alternative base in normal sample, 
    ### read depth of alternative base in tumor sample,
    ### allele frequency of alternative base in normal sample, allele frequency of alternative base in tumor sample,
    ### FILTER, average median mapping quality of reads in normal and tumor samples

    snv.df <- data.frame(chr = vcf$'#CHROM', pos = vcf$POS, ref = vcf$REF, alt = vcf$ALT, filter = vcf$FILTER)

    ### add average median mapping quality 
    MQ <- (as.numeric(system(paste("bcftools query -f '%INFO/MMQ\n'", vcf.snvpaths[i], "| awk -F, '{print $1}'"), intern = TRUE)) +
        as.numeric(system(paste("bcftools query -f '%INFO/MMQ\n'", vcf.snvpaths[i], "| awk -F, '{print $2}'"), intern = TRUE)))/2
    
    snv.df <- cbind(snv.df, mq = as.numeric(MQ))

    ### remove MQ variable
    rm(MQ)

    ### estimate number of reads supporting reference and alternative base in tumor and normal samples
    norm_ref_depth <- as.numeric(system(paste("bcftools query -s ", normal.id, " -f '[%AD]\n'", vcf.snvpaths[i], "| awk -F, '{print $1}'"), intern = TRUE))
    norm_alt_depth <-  as.numeric(system(paste("bcftools query -s ", normal.id, " -f '[%AD]\n'", vcf.snvpaths[i], "| awk -F, '{print $2}'"), intern = TRUE))
    tumor_ref_depth <-  as.numeric(system(paste("bcftools query -s ", tumor.id, " -f '[%AD]\n'", vcf.snvpaths[i], "| awk -F, '{print $1}'"), intern = TRUE))
    tumor_alt_depth <- as.numeric(system(paste("bcftools query -s ", tumor.id, " -f '[%AD]\n'", vcf.snvpaths[i], "| awk -F, '{print $2}'"), intern = TRUE))

    ### add depths and allele frequencies
    snv.df <- cbind(snv.df, totalnormal = norm_ref_depth+norm_alt_depth, normalrefdepth = norm_ref_depth, 
        normalaltdepth = norm_alt_depth, totaltumor = tumor_ref_depth + tumor_alt_depth, 
        tumorrefdepth = tumor_ref_depth, tumoraltdepth = tumor_alt_depth, afaltnormal = norm_alt_depth/(norm_alt_depth + norm_ref_depth),
        afalttumor = tumor_alt_depth/(tumor_alt_depth + tumor_ref_depth))
    snv.df$totaldepth <- snv.df$totalnormal+snv.df$totaltumor

     ### exclude multiallelic variants
    snv.df <- snv.df[which(!grepl("multiallelic", snv.df$filter)) ,]

    ### set id for sample
    id <- ANNO_FILE$subject_id[i]

    plot.dir <- file.path(output.snvdir, "plots")
    if (!file.exists(plot.dir)){
        dir.create(plot.dir)
    }
    setwd(plot.dir)

    ### make plots: 
    #### log10 distribution of totaldepth with cutoff at 40
    hist_plot <- hist(log10(snv.df$totaldepth), xlim = c(0, log10(max(snv.df$totaldepth)+1)), col = "skyblue", main = id, xlab = "Read depth (log10)", ylab = "Frequency")
    png(paste0("Histogram.allDP.mutect2.snvs.", id, ".png"), res=200, unit="in", height=8, width=11)
    plot(hist_plot, col = "skyblue", main = id, xlab = "Read depth (log10)", ylab = "Frequency")
    abline(v = log10(opt$'total-rd-threshold-snv'), col = "red", lty = 2)
    dev.off()

    #### average mapping quality distribution with cutoff at 30
    hist_plot <- hist(snv.df$mq, xlim = c(0, max(snv.df$mq)), col = "skyblue", main = id, xlab = "Mapping quality", ylab = "Frequency")
    png(paste0("Histogram.MQ.mutect2.snvs.", id, ".png"), res=200, unit="in", height=8, width=11)
    plot(hist_plot, col = "skyblue", main = id, xlab = "Mapping quality", ylab = "Frequency")
    abline(v = opt$'mq-threshold-snv', col = "red", lty = 2)
    dev.off()

    #### log10 distribution of totaltumor with cutoff at 25
    hist_plot <- hist(log10(snv.df$totaltumor), xlim = c(0, log10(max(snv.df$totaltumor)+1)), col = "skyblue", main = id, xlab = "Read depth in tumor (log10)", ylab = "Frequency")
    png(paste0("Histogram.TumorDP.mutect2.snvs.", id, ".png"), res=200, unit="in", height=8, width=11)
    plot(hist_plot, col = "skyblue", main = id, xlab = "Read depth in tumor (log10)", ylab = "Frequency")
    abline(v = log10(opt$'tumor-rd-threshold-snv'), col = "red", lty = 2)
    dev.off()

    #### log10 distribution of totalnormal with cutoff at 15
    hist_plot <- hist(log10(snv.df$totalnormal), xlim = c(0, log10(max(snv.df$totalnormal)+1)), col = "skyblue", main = id, xlab = "Read depth in normal (log10)", ylab = "Frequency")
    png(paste0("Histogram.NormalDP.mutect2.snvs.", id, ".png"), res=200, unit="in", height=8, width=11)
    plot(hist_plot, col = "skyblue", main = id, xlab = "Read depth in normal (log10)", ylab = "Frequency")
    abline(v = log10(opt$'normal-rd-threshold-snv'), col = "red", lty = 2)
    dev.off()

    snv.sub.df <- snv.df[snv.df$totaltumor > opt$'tumor-rd-threshold-snv' & snv.df$totalnormal > opt$'normal-rd-threshold-snv'
        & snv.df$mq > opt$'mq-threshold-snv',]
    cat("Number of variants left after read depth filtering ", nrow(snv.sub.df), "\n")

    ### continue with the plots
    #### afalttumor vs afaltnormal before read depth filter
    png(paste0("Scatterplot.AFAltTumorvsAFAltNormal.mutect2.snvs.beforeRDfilter.", id, ".png"), res=200, unit="in", height=8, width=11)
    print({plot1 <- ggplot(snv.df, aes(x=afaltnormal, y=afalttumor)) +
        geom_point() +
        xlab("AFAlt_in_normal") +
        ylab("AFAlt_in_tumor") +
        xlim(0,1) +
        ylim(0,1) +
        geom_vline(xintercept = as.numeric(opt$'normal-afalt-snv'), linetype = "dashed", color = "red") +
        geom_hline(yintercept = as.numeric(opt$'tumor-afalt-snv'), linetype = "dashed", color = "red")
        plot1})
    dev.off()

    #### afalttumor vs afaltnormal after read depth filter
    png(paste0("Scatterplot.AFAltTumorvsAFAltNormal.mutect2.snvs.afterRDfilter.", id, ".png"), res=200, unit="in", height=8, width=11)
    print({plot1 <- ggplot(snv.sub.df, aes(x=afaltnormal, y=afalttumor)) +
        geom_point() +
        xlab("AFAlt_in_normal") +
        ylab("AFAlt_in_tumor") +
        xlim(0,1) +
        ylim(0,1) +
        geom_vline(xintercept = as.numeric(opt$'normal-afalt-snv'), linetype = "dashed", color = "red") +
        geom_hline(yintercept = as.numeric(opt$'tumor-afalt-snv'), linetype = "dashed", color = "red")
        plot1})
    dev.off()

    #### afalttumor vs tumoraltdepth before read depth filter
    png(paste0("Scatterplot.AFAltTumorvsDPAltTumor.mutect2.snvs.beforeRDfilter.", id, ".png"), res=200, unit="in", height=8, width=11)
    print({plot1 <- ggplot(snv.df, aes(x=log10(tumoraltdepth+1), y=afalttumor)) +
        geom_point() +
        xlab("DPAltTumor (log10)") +
        ylab("AFAlt_in_tumor") +
        ylim(0,1) +
        geom_vline(xintercept = log10(opt$'tumor-rdalt-snv'), linetype = "dashed", color = "red") +
        geom_hline(yintercept = as.numeric(opt$'tumor-afalt-snv'), linetype = "dashed", color = "red")
        plot1})
    dev.off()
    
    #### afalttumor vs tumoraltdepth after read depth filter
    png(paste0("Scatterplot.AFAltTumorvsDPAltTumor.mutect2.snvs.afterRDfilter.", id, ".png"), res=200, unit="in", height=8, width=11)
    print({plot1 <- ggplot(snv.sub.df, aes(x=log10(tumoraltdepth+1), y=afalttumor)) +
        geom_point() +
        xlab("DPAltTumor (log10)") +
        ylab("AFAlt_in_tumor") +
        ylim(0,1) +
        geom_vline(xintercept = log10(opt$'tumor-rdalt-snv'), linetype = "dashed", color = "red") +
        geom_hline(yintercept = as.numeric(opt$'tumor-afalt-snv'), linetype = "dashed", color = "red")
        plot1})
    dev.off()

    #### afalttumor vs normalaltdepth before read depth filter
    png(paste0("Scatterplot.AFAltTumorvsDPAltNormal.mutect2.snvs.beforeRDfilter.", id, ".png"), res=200, unit="in", height=8, width=11)
    print({plot1 <- ggplot(snv.df, aes(x=log10(normalaltdepth+1), y=afalttumor)) +
        geom_point() +
        xlab("DPAltNormal (log10)") +
        ylab("AFAlt_in_tumor") +
        ylim(0,1) +
        geom_vline(xintercept = log10(opt$'normal-rdalt-snv'), linetype = "dashed", color = "red") +
        geom_hline(yintercept = as.numeric(opt$'tumor-afalt-snv'), linetype = "dashed", color = "red")
        plot1})
    dev.off()
    
    #### afalttumor vs normalaltdepth after read depth filter
    png(paste0("Scatterplot.AFAltTumorvsDPAltNormal.mutect2.snvs.afterRDfilter.", id, ".png"), res=200, unit="in", height=8, width=11)
    print({plot1 <- ggplot(snv.sub.df, aes(x=log10(normalaltdepth+1), y=afalttumor)) +
        geom_point() +
        xlab("DPAltNormal (log10)") +
        ylab("AFAlt_in_tumor") +
        ylim(0,1) +
        geom_vline(xintercept = log10(opt$'normal-rdalt-snv'), linetype = "dashed", color = "red") +
        geom_hline(yintercept = as.numeric(opt$'tumor-afalt-snv'), linetype = "dashed", color = "red")
        plot1})
    dev.off()

    snv.sub.df <- snv.sub.df[snv.sub.df$afaltnormal < opt$'normal-afalt-snv' & snv.sub.df$afalttumor > opt$'tumor-afalt-snv' 
        & snv.sub.df$tumoraltdepth >= opt$'tumor-rdalt-snv' & snv.sub.df$normalaltdepth < opt$'normal-rdalt-snv',]
    cat("Number of variants left after read depth filtering ", nrow(snv.sub.df), "\n")

    snv.list[[i]] <- data.frame(sample = id, before = nrow(snv.df), after = nrow(snv.sub.df))

    if (opt$'output-vcf' == 1){

        output.file <- file.path(vcf.output.dir, paste(id, "_hc_snv_variants.vcf"))
        vcf <- vcf[match(paste(snv.sub.df[,1], snv.sub.df[,2], snv.sub.df[,3], snv.sub.df[,4], sep=""), 
            paste(vcf[,1], vcf[,2], vcf[,4], vcf[,5], sep="")),]
        ### add vcf header
        system(paste("grep '^#'", vcf.snvpaths[i], ">", output.file, sep=" "))
        ### add rest of variants
        writeLines(vcf, con=output.file)
    }

    #### remove intermediate SNV file
    system(paste("rm", vcf.snvpaths[i], sep=" "))

    cat("Finished processing sample ", id, "\n")
}

### save counts of SNV variants
write.table(do.call(rbind, snv.list), file = "Nmut@snvs.mutect2.txt", sep = "\t", row.names = FALSE)

end_time <- Sys.time()
cat("Elapsed minutes for processing SNVs", end_time - start_time)


####################### INDEL filtering
#### set working directory
setwd(output.indeldir)

indel.list <- list()

start_time <- Sys.time()
for (i in 1:length(vcf.indelpaths)){
    cat("Processing", vcf.indelpaths[i], "\n")

    tumor.id <- ANNO_FILE$tumor_id[i]
    normal.id <- ANNO_FILE$normal_id[i]

    ### read vcf file
    vcf <- fread(vcf.indelpaths[i], data.table=FALSE, skip="#CHROM")

    ### Combine all INDELs into single dataframe with following columns:
    ### chromosome, position, reference base, alternative base, total read depth in normal sample excl filtered reads,
    ### total read depth in tumor sample excl filtered reads, read depth of reference base in normal sample
    ### read depth of reference base in tumor sample, read depth of alternative base in normal sample, 
    ### read depth of alternative base in tumor sample,
    ### allele frequency of alternative base in normal sample, allele frequency of alternative base in tumor sample,
    ### FILTER, average median mapping quality of reads in normal and tumor samples

    indel.df <- data.frame(chr = vcf$'#CHROM', pos = vcf$POS, ref = vcf$REF, alt = vcf$ALT, filter = vcf$FILTER)

    ### add average median mapping quality 
    MQ <- (as.numeric(system(paste("bcftools query -f '%INFO/MMQ\n'", vcf.indelpaths[i], "| awk -F, '{print $1}'"), intern = TRUE)) +
        as.numeric(system(paste("bcftools query -f '%INFO/MMQ\n'", vcf.indelpaths[i], "| awk -F, '{print $2}'"), intern = TRUE)))/2
    
    indel.df <- cbind(indel.df, mq = as.numeric(MQ))

    ### remove MQ variable
    rm(MQ)

    ### estimate number of reads supporting reference and alternative base in tumor and normal samples
    norm_ref_depth <- as.numeric(system(paste("bcftools query -s ", normal.id, " -f '[%AD]\n'", vcf.indelpaths[i], "| awk -F, '{print $1}'"), intern = TRUE))
    norm_alt_depth <-  as.numeric(system(paste("bcftools query -s ", normal.id, " -f '[%AD]\n'", vcf.indelpaths[i], "| awk -F, '{print $2}'"), intern = TRUE))
    tumor_ref_depth <-  as.numeric(system(paste("bcftools query -s ", tumor.id, " -f '[%AD]\n'", vcf.indelpaths[i], "| awk -F, '{print $1}'"), intern = TRUE))
    tumor_alt_depth <- as.numeric(system(paste("bcftools query -s ", tumor.id, " -f '[%AD]\n'", vcf.indelpaths[i], "| awk -F, '{print $2}'"), intern = TRUE))

    ### add depths and allele frequencies
    indel.df <- cbind(indel.df, totalnormal = norm_ref_depth+norm_alt_depth, normalrefdepth = norm_ref_depth, 
        normalaltdepth = norm_alt_depth, totaltumor = tumor_ref_depth + tumor_alt_depth, 
        tumorrefdepth = tumor_ref_depth, tumoraltdepth = tumor_alt_depth, afaltnormal = norm_alt_depth/(norm_alt_depth + norm_ref_depth),
        afalttumor = tumor_alt_depth/(tumor_alt_depth + tumor_ref_depth))
    indel.df$totaldepth <- indel.df$totalnormal+indel.df$totaltumor

     ### exclude multiallelic variants
    indel.df <- indel.df[which(!grepl("multiallelic", indel.df$filter)) ,]

    ### set id for sample
    id <- ANNO_FILE$subject_id[i]

    plot.dir <- file.path(output.indeldir, "plots")
    if (!file.exists(plot.dir)){
        dir.create(plot.dir)
    }
    setwd(plot.dir)

    ### make plots: 
    #### log10 distribution of totaldepth with cutoff at 40
    hist_plot <- hist(log10(indel.df$totaldepth), xlim = c(0, log10(max(indel.df$totaldepth)+1)), col = "skyblue", main = id, xlab = "Read depth (log10)", ylab = "Frequency")
    png(paste0("Histogram.allDP.mutect2.indels.", id, ".png"), res=200, unit="in", height=8, width=11)
    plot(hist_plot, col = "skyblue", main = id, xlab = "Read depth (log10)", ylab = "Frequency")
    abline(v = log10(opt$'total-rd-threshold-indel'), col = "red", lty = 2)
    dev.off()

    #### average mapping quality distribution with cutoff at 30
    hist_plot <- hist(indel.df$mq, xlim = c(0, max(indel.df$mq)), col = "skyblue", main = id, xlab = "Mapping quality", ylab = "Frequency")
    png(paste0("Histogram.MQ.mutect2.indels.", id, ".png"), res=200, unit="in", height=8, width=11)
    plot(hist_plot, col = "skyblue", main = id, xlab = "Mapping quality", ylab = "Frequency")
    abline(v = opt$'mq-threshold-indel', col = "red", lty = 2)
    dev.off()

    #### log10 distribution of totaltumor with cutoff at 25
    hist_plot <- hist(log10(indel.df$totaltumor), xlim = c(0, log10(max(indel.df$totaltumor)+1)), col = "skyblue", main = id, xlab = "Read depth in tumor (log10)", ylab = "Frequency")
    png(paste0("Histogram.TumorDP.mutect2.indels.", id, ".png"), res=200, unit="in", height=8, width=11)
    plot(hist_plot, col = "skyblue", main = id, xlab = "Read depth in tumor (log10)", ylab = "Frequency")
    abline(v = log10(opt$'tumor-rd-threshold-indel'), col = "red", lty = 2)
    dev.off()

    #### log10 distribution of totalnormal with cutoff at 15
    hist_plot <- hist(log10(indel.df$totalnormal), xlim = c(0, log10(max(indel.df$totalnormal)+1)), col = "skyblue", main = id, xlab = "Read depth in normal (log10)", ylab = "Frequency")
    png(paste0("Histogram.NormalDP.mutect2.indels.", id, ".png"), res=200, unit="in", height=8, width=11)
    plot(hist_plot, col = "skyblue", main = id, xlab = "Read depth in normal (log10)", ylab = "Frequency")
    abline(v = log10(opt$'normal-rd-threshold-indel'), col = "red", lty = 2)
    dev.off()

    indel.sub.df <- indel.df[indel.df$totaltumor > opt$'tumor-rd-threshold-indel' & indel.df$totalnormal > opt$'normal-rd-threshold-indel'
        & indel.df$mq > opt$'mq-threshold-indel',]

    cat("Number of variants left after read depth filtering ", nrow(indel.sub.df), "\n")

    ### continue with the plots
    #### afalttumor vs afaltnormal before read depth filter
    png(paste0("Scatterplot.AFAltTumorvsAFAltNormal.mutect2.indels.beforeRDfilter.", id, ".png"), res=200, unit="in", height=8, width=11)
    print({plot1 <- ggplot(indel.df, aes(x=afaltnormal, y=afalttumor)) +
        geom_point() +
        xlab("AFAlt_in_normal") +
        ylab("AFAlt_in_tumor") +
        xlim(0,1) +
        ylim(0,1) +
        geom_vline(xintercept = as.numeric(opt$'normal-afalt-indel'), linetype = "dashed", color = "red") +
        geom_hline(yintercept = as.numeric(opt$'tumor-afalt-indel'), linetype = "dashed", color = "red")
        plot1})
    dev.off()

    #### afalttumor vs afaltnormal after read depth filter
    png(paste0("Scatterplot.AFAltTumorvsAFAltNormal.mutect2.indels.afterRDfilter.", id, ".png"), res=200, unit="in", height=8, width=11)
    print({plot1 <- ggplot(indel.sub.df, aes(x=afaltnormal, y=afalttumor)) +
        geom_point() +
        xlab("AFAlt_in_normal") +
        ylab("AFAlt_in_tumor") +
        xlim(0,1) +
        ylim(0,1) +
        geom_vline(xintercept = as.numeric(opt$'normal-afalt-indel'), linetype = "dashed", color = "red") +
        geom_hline(yintercept = as.numeric(opt$'tumor-afalt-indel'), linetype = "dashed", color = "red")
        plot1})
    dev.off()

    #### afalttumor vs tumoraltdepth before read depth filter
    png(paste0("Scatterplot.AFAltTumorvsDPAltTumor.mutect2.indels.beforeRDfilter.", id, ".png"), res=200, unit="in", height=8, width=11)
    print({plot1 <- ggplot(indel.df, aes(x=log10(tumoraltdepth+1), y=afalttumor)) +
        geom_point() +
        xlab("DPAltTumor (log10)") +
        ylab("AFAlt_in_tumor") +
        ylim(0,1) +
        geom_vline(xintercept = log10(opt$'tumor-rdalt-indel'), linetype = "dashed", color = "red") +
        geom_hline(yintercept = as.numeric(opt$'tumor-afalt-indel'), linetype = "dashed", color = "red")
        plot1})
    dev.off()
    
    #### afalttumor vs tumoraltdepth after read depth filter
    png(paste0("Scatterplot.AFAltTumorvsDPAltTumor.mutect2.indels.afterRDfilter.", id, ".png"), res=200, unit="in", height=8, width=11)
    print({plot1 <- ggplot(indel.sub.df, aes(x=log10(tumoraltdepth+1), y=afalttumor)) +
        geom_point() +
        xlab("DPAltTumor (log10)") +
        ylab("AFAlt_in_tumor") +
        ylim(0,1) +
        geom_vline(xintercept = log10(opt$'tumor-rdalt-indel'), linetype = "dashed", color = "red") +
        geom_hline(yintercept = as.numeric(opt$'tumor-afalt-indel'), linetype = "dashed", color = "red")
        plot1})
    dev.off()

    #### afalttumor vs normalaltdepth before read depth filter
    png(paste0("Scatterplot.AFAltTumorvsDPAltNormal.mutect2.indels.beforeRDfilter.", id, ".png"), res=200, unit="in", height=8, width=11)
    print({plot1 <- ggplot(indel.df, aes(x=log10(normalaltdepth+1), y=afalttumor)) +
        geom_point() +
        xlab("DPAltNormal (log10)") +
        ylab("AFAlt_in_tumor") +
        ylim(0,1) +
        geom_vline(xintercept = log10(opt$'normal-rdalt-indel'), linetype = "dashed", color = "red") +
        geom_hline(yintercept = as.numeric(opt$'tumor-afalt-indel'), linetype = "dashed", color = "red")
        plot1})
    dev.off()
    
    #### afalttumor vs normalaltdepth after read depth filter
    png(paste0("Scatterplot.AFAltTumorvsDPAltNormal.mutect2.indels.afterRDfilter.", id, ".png"), res=200, unit="in", height=8, width=11)
    print({plot1 <- ggplot(indel.sub.df, aes(x=log10(normalaltdepth+1), y=afalttumor)) +
        geom_point() +
        xlab("DPAltNormal (log10)") +
        ylab("AFAlt_in_tumor") +
        ylim(0,1) +
        geom_vline(xintercept = log10(opt$'normal-rdalt-indel'), linetype = "dashed", color = "red") +
        geom_hline(yintercept = as.numeric(opt$'tumor-afalt-indel'), linetype = "dashed", color = "red")
        plot1})
    dev.off()

    indel.sub.df <- indel.sub.df[indel.sub.df$afaltnormal < opt$'normal-afalt-indel' & indel.sub.df$afalttumor >= as.numeric(opt$'tumor-afalt-indel') 
        & indel.sub.df$tumoraltdepth >= opt$'tumor-rdalt-indel' & indel.sub.df$normalaltdepth < opt$'normal-rdalt-indel',] 
    cat("Number of variants left after read depth filtering ", nrow(indel.sub.df), "\n")

    indel.list[[i]] <- data.frame(sample = id, before = nrow(indel.df), after = nrow(indel.sub.df))

    if (opt$'output-vcf' == 1){

        output.file <- file.path(vcf.output.dir, paste(id, "_hc_indel_variants.vcf"))
        vcf <- vcf[match(paste(indel.sub.df[,1], indel.sub.df[,2], indel.sub.df[,3], indel.sub.df[,4], sep=""), 
            paste(vcf[,1], vcf[,2], vcf[,4], vcf[,5], sep="")),]
        ### add vcf header
        system(paste("grep '^#'", vcf.indelpaths[i], ">", output.file, sep=" "))
        ### add rest of variants
        writeLines(vcf, con=output.file)
    }

    #### remove intermediate SNV file
    system(paste("rm", vcf.indelpaths[i], sep=" "))

    cat("Finished processing sample ", id, "\n")
}

### save counts of INDEL variants
write.table(do.call(rbind, indel.list), file = "Nmut@indels.mutect2.txt", sep = "\t", row.names = FALSE)

end_time <- Sys.time()
cat("Elapsed minutes for processing INDELs", end_time - start_time)