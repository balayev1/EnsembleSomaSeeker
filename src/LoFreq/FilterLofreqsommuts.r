############################################### FilterLofreqsommuts.r
######################## This script performs graphical analysis of Strelka2 output vcf files from matched tumor-normal analysis
######################## and returns diagnostic plots and filtered variants in vcf format
################ Inputs:
################ @param anno-file: full path to variant annotation file
################ @param output-dir: output directory
################ @param tumor-rd-threshold-snv: tumor read depth threshold for SNVs
################ @param tumor-rd-threshold-indel: tumor read depth threshold for INDELs
################ @param tumor-afalt-snv: allele frequency threshold of alternate allele in tumor for SNVs
################ @param tumor-afalt-indel: allele frequency threshold of alternate allele in tumor for INDELs
################ @param tumor-rdalt-snv: read depth threshold of alternate allele in tumor for SNVs
################ @param tumor-rdalt-indel: read depth threshold of alternate allele in tumor for INDELs
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
    make_option(c("-a", "--tumor-rd-threshold-snv"), type="numeric", default=25,
              help="tumor read depth threshold for SNVs", metavar="numeric"),
    make_option(c("-A", "--tumor-rd-threshold-indel"), type="numeric", default=25, 
              help="tumor read depth threshold for INDELs", metavar="numeric"),
    make_option(c("-t", "--tumor-afalt-snv"), type="numeric", default=0.05,
              help="allele frequency threshold of alternate allele in tumor for SNVs", metavar="numeric"),
    make_option(c("-T", "--tumor-afalt-indel"), type="numeric", default=0.05,
              help="allele frequency threshold of alternate allele in tumor for INDELs", metavar="numeric"),
    make_option(c("-p", "--tumor-rdalt-snv"), type="numeric", default=5,
              help="read depth threshold of alternate allele in tumor for SNVs", metavar="numeric"),
    make_option(c("-P", "--tumor-rdalt-indel"), type="numeric", default=5,
            help="read depth threshold of alternate allele in tumor for INDELs", metavar="numeric"),
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

### build output structure
### make folder to output filtering vcfs with SNVs
output.dir <- file.path(output.dir, "Filtered_lofreqvcfs")
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
ANNO_FILE <- cbind(ANNO_FILE[, c(1,2,3)],  ANNO_FILE$lofreq_snv_vcf, ANNO_FILE$lofreq_indel_vcf)
colnames(ANNO_FILE) <- c("subject_id", "tumor_id", "normal_id", "snv", "indel")

### replace blank regions with NA
ANNO_FILE[ANNO_FILE == ""] <- NA

####################### SNV filtering
#### set working directory
setwd(output.snvdir)

### list paths to sample-specific (i.e. tumor-normal pair) vcf SNV files
vcf.snvpaths <- c()
for (k in 1:nrow(ANNO_FILE)){
    if (!is.na(ANNO_FILE$snv[k]) & file.exists(ANNO_FILE$snv[k])){
        vcf.snvpaths <- append(vcf.snvpaths, ANNO_FILE$snv[k])
    }
}
cat("Found total", length(vcf.snvpaths), "existing SNV vcf files\n")

snv.list <- list()
start_time <- Sys.time()
for (i in 1:length(vcf.snvpaths)){
    cat("Processing", vcf.snvpaths[i], "\n")

    ### read vcf file
    vcf <- fread(vcf.snvpaths[i], data.table=FALSE, skip="#CHROM")


    ### Combine all SNVs into single dataframe with following columns:
    ### chromosome, position, reference base, alternative base,
    ### total read depth in tumor sample excl filtered reads,
    ### read depth of reference base in tumor sample,
    ### read depth of alternative base in tumor sample,
    ### allele frequency of alternative base in tumor sample,
    ### FILTER status

    snv.df <- data.frame(chr = vcf$'#CHROM', pos = vcf$POS, ref = vcf$REF, alt = vcf$ALT, filter = vcf$FILTER)
    
    ### add number of reads supporting each base in tumor samples to dataframe
    snv.df$ref_tumor <- as.numeric(system(paste("bcftools query -f '%DP4\n'", vcf.snvpaths[i], "| awk -F, '{print $1}'"), intern = TRUE)) +
        as.numeric(system(paste("bcftools query -f '%DP4\n'", vcf.snvpaths[i], "| awk -F, '{print $2}'"), intern = TRUE))
    snv.df$alt_tumor <- as.numeric(system(paste("bcftools query -f '%DP4\n'", vcf.snvpaths[i], "| awk -F, '{print $3}'"), intern = TRUE)) +
        as.numeric(system(paste("bcftools query -f '%DP4\n'", vcf.snvpaths[i], "| awk -F, '{print $4}'"), intern = TRUE))

    ### add depths and allele frequencies
    snv.df <- cbind(snv.df, totaltumor = snv.df$ref_tumor + snv.df$alt_tumor, 
        tumorrefdepth = snv.df$ref_tumor, tumoraltdepth = snv.df$alt_tumor,
        afalttumor = snv.df$alt_tumor/(snv.df$alt_tumor + snv.df$ref_tumor))

    ### set id for sample
    id <- ANNO_FILE$subject_id[i]

    plot.dir <- file.path(output.snvdir, "plots")
    if (!file.exists(plot.dir)){
        dir.create(plot.dir)
    }
    setwd(plot.dir)

    ### make plots: 
    #### log10 distribution of totaltumor with cutoff at 25
    hist_plot <- hist(log10(snv.df$totaltumor), xlim = c(0, log10(max(snv.df$totaltumor)+1)), col = "skyblue", main = id, xlab = "Read depth in tumor (log10)", ylab = "Frequency")
    png(paste0("Histogram.TumorDP.lofreq.snvs.", id, ".png"), res=200, unit="in", height=8, width=11)
    plot(hist_plot, col = "skyblue", main = id, xlab = "Read depth in tumor (log10)", ylab = "Frequency")
    abline(v = log10(opt$'tumor-rd-threshold-snv'), col = "red", lty = 2)
    dev.off()

    snv.sub.df <- snv.df[snv.df$totaltumor > opt$'tumor-rd-threshold-snv',]
    cat("Number of variants left after read depth filtering ", nrow(snv.sub.df), "\n")

    ### continue with the plots
    #### afalttumor vs tumoraltdepth before read depth filter
    png(paste0("Scatterplot.AFAltTumorvsDPAltTumor.lofreq.snvs.beforeRDfilter.", id, ".png"), res=200, unit="in", height=8, width=11)
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
    png(paste0("Scatterplot.AFAltTumorvsDPAltTumor.lofreq.snvs.afterRDfilter.", id, ".png"), res=200, unit="in", height=8, width=11)
    print({plot1 <- ggplot(snv.sub.df, aes(x=log10(tumoraltdepth+1), y=afalttumor)) +
        geom_point() +
        xlab("DPAltTumor (log10)") +
        ylab("AFAlt_in_tumor") +
        ylim(0,1) +
        geom_vline(xintercept = log10(opt$'tumor-rdalt-snv'), linetype = "dashed", color = "red") +
        geom_hline(yintercept = as.numeric(opt$'tumor-afalt-snv'), linetype = "dashed", color = "red")
        plot1})
    dev.off()

    snv.sub.df <- snv.sub.df[snv.sub.df$afalttumor > opt$'tumor-afalt-snv' & snv.sub.df$tumoraltdepth >= opt$'tumor-rdalt-snv',]
    cat("Number of variants left after alternate allele stats filtering ", nrow(snv.sub.df), "\n")

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

    cat("Finished processing sample ", id, "\n")
}

### save counts of SNV variants
write.table(do.call(rbind, snv.list), file = "Nmut@snvs.lofreq.txt", sep = "\t", row.names = FALSE)

end_time <- Sys.time()
cat("Elapsed minutes for processing SNVs", end_time - start_time)



####################### INDEL filtering
#### set working directory
setwd(output.indeldir)

### list paths to sample-specific (i.e. tumor-normal pair) vcf INDEL files
vcf.indelpaths <- c()
for (k in 1:nrow(ANNO_FILE)){
    if (!is.na(ANNO_FILE$indel[k]) & file.exists(ANNO_FILE$indel[k])){
        vcf.indelpaths <- append(vcf.indelpaths, ANNO_FILE$indel[k])
    }
}
cat("Found total", length(vcf.indelpaths), "existing INDEL vcf files\n")


indel.list <- list()
start_time <- Sys.time()
for (i in 1:length(vcf.indelpaths)){
    cat("Processing", vcf.indelpaths[i], "\n")

    ### read vcf file
    vcf <- fread(vcf.indelpaths[i], data.table=FALSE, skip="#CHROM")


    ### Combine all INDELs into single dataframe with following columns:
    ### chromosome, position, reference base, alternative base,
    ### total read depth in tumor sample excl filtered reads,
    ### read depth of reference base in tumor sample,
    ### read depth of alternative base in tumor sample,
    ### allele frequency of alternative base in tumor sample,
    ### FILTER status

    indel.df <- data.frame(chr = vcf$'#CHROM', pos = vcf$POS, ref = vcf$REF, alt = vcf$ALT, filter = vcf$FILTER)

    ### add number of reads supporting each base in tumor samples to dataframe
    indel.df$ref_tumor <- as.numeric(system(paste("bcftools query -f '%DP4\n'", vcf.indelpaths[i], "| awk -F, '{print $1}'"), intern = TRUE)) +
        as.numeric(system(paste("bcftools query -f '%DP4\n'", vcf.indelpaths[i], "| awk -F, '{print $2}'"), intern = TRUE))
    indel.df$alt_tumor <- as.numeric(system(paste("bcftools query -f '%DP4\n'", vcf.indelpaths[i], "| awk -F, '{print $3}'"), intern = TRUE)) +
        as.numeric(system(paste("bcftools query -f '%DP4\n'", vcf.indelpaths[i], "| awk -F, '{print $4}'"), intern = TRUE))

    ### add depths and allele frequencies
    indel.df <- cbind(indel.df, totaltumor = indel.df$ref_tumor + indel.df$alt_tumor, 
        tumorrefdepth = indel.df$ref_tumor, tumoraltdepth = indel.df$alt_tumor,
        afalttumor = indel.df$alt_tumor/(indel.df$alt_tumor + indel.df$ref_tumor))

    ### set id for sample
    id <- ANNO_FILE$subject_id[i]

    plot.dir <- file.path(output.indeldir, "plots")
    if (!file.exists(plot.dir)){
        dir.create(plot.dir)
    }
    setwd(plot.dir)

    ### make plots: 
    #### log10 distribution of totaltumor with cutoff at 25
    hist_plot <- hist(log10(indel.df$totaltumor), xlim = c(0, log10(max(indel.df$totaltumor)+1)), col = "skyblue", main = id, xlab = "Read depth in tumor (log10)", ylab = "Frequency")
    png(paste0("Histogram.TumorDP.lofreq.indels.", id, ".png"), res=200, unit="in", height=8, width=11)
    plot(hist_plot, col = "skyblue", main = id, xlab = "Read depth in tumor (log10)", ylab = "Frequency")
    abline(v = log10(opt$'tumor-rd-threshold-indel'), col = "red", lty = 2)
    dev.off()

    indel.sub.df <- indel.df[indel.df$totaltumor > opt$'tumor-rd-threshold-indel',]
    cat("Number of variants left after read depth filtering ", nrow(indel.sub.df), "\n")

    ### continue with the plots
    #### afalttumor vs tumoraltdepth before read depth filter
    png(paste0("Scatterplot.AFAltTumorvsDPAltTumor.lofreq.indels.beforeRDfilter.", id, ".png"), res=200, unit="in", height=8, width=11)
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
    png(paste0("Scatterplot.AFAltTumorvsDPAltTumor.lofreq.indels.afterRDfilter.", id, ".png"), res=200, unit="in", height=8, width=11)
    print({plot1 <- ggplot(indel.sub.df, aes(x=log10(tumoraltdepth+1), y=afalttumor)) +
        geom_point() +
        xlab("DPAltTumor (log10)") +
        ylab("AFAlt_in_tumor") +
        ylim(0,1) +
        geom_vline(xintercept = log10(opt$'tumor-rdalt-indel'), linetype = "dashed", color = "red") +
        geom_hline(yintercept = as.numeric(opt$'tumor-afalt-indel'), linetype = "dashed", color = "red")
        plot1})
    dev.off()

    indel.sub.df <- indel.sub.df[indel.sub.df$afalttumor > opt$'tumor-afalt-indel'
        & indel.sub.df$tumoraltdepth >= opt$'tumor-rdalt-indel', ]
    cat("Number of variants left after alternate allele stats filtering ", nrow(indel.sub.df), "\n")

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
    cat("Finished processing sample ", id, "\n")
}

### save counts of INDEL variants
write.table(do.call(rbind, indel.list), file = "Nmut@indels.lofreq.txt", sep = "\t", row.names = FALSE)

end_time <- Sys.time()
cat("Elapsed minutes for processing INDELs", end_time - start_time)