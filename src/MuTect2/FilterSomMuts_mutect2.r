######################## Filtering of Mutect2 output vcf files
######################## vcf files are the output from running Varscan command in Mutect2SomMuts_runscript.sh file

### load required libraries
require(data.table)
require(ggplot2)

### define parent directory  
parent.dir <- "/home/aventeic/balay011/WGS_skChordomas_Bai/SomaticSNV_calling/mutect2"
if (!dir.exists(parent.dir)) {
    stop("Parent directory with VCF files does not exist. Please check\n")
}

### build output structure
### make folder to output filtering vcfs with SNVs
output.dir <- file.path(parent.dir, "Filtered_vcfs")
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

### list paths to sample-specific (i.e. tumor-normal pair) vcf files
vcf.paths <- as.vector(list.files(parent.dir, pattern = "*_filtered.vcf.gz$", full.names=TRUE, recursive=TRUE))
length(vcf.paths)
# [1] 80

### function that splits each vcf file into SNVs and INDELs
picard.path <- "~/.conda/envs/gatk_env/bin/picard"
splitvcfs <- function(vcf.path){
    vcf.snvpath <- gsub("_filtered.vcf.gz", "_snv.vcf.gz", vcf.path)
    vcf.indelpath <- gsub("_filtered.vcf.gz", "_indel.vcf.gz", vcf.path)
    system(paste(picard.path, "SplitVcfs I=", vcf.path, "SNP_OUTPUT=", vcf.snvpath,
        "INDEL_OUTPUT=", vcf.indelpath, "CREATE_INDEX=true", "STRICT=false", sep = " "))
}

### load annotation file
ANNO_FILE <- read.table("/home/aventeic/balay011/scripts/WGS_sommuts_anno.txt", sep="\t")
colnames(ANNO_FILE) <- c("run_id", "is_tumor", "vcf_path")

####################### SNV filtering
#### set working directory
setwd(output.snvdir)

snv.list <- list()

start_time <- Sys.time()
for (i in 1:length(vcf.paths)){
    cat("Processing", vcf.paths[i], "for SNVs\n")

    tumor.id <- ANNO_FILE$run_id[ANNO_FILE$is_tumor == 'Yes' & ANNO_FILE$vcf_path == basename(vcf.paths[i])]
    normal.id <- ANNO_FILE$run_id[ANNO_FILE$is_tumor == 'No' & ANNO_FILE$vcf_path == basename(vcf.paths[i])]

    ### split variants into SNVs and INDELs from vcf file
    splitvcfs(vcf.paths[i])

    vcf.snvpath <- gsub("_filtered.vcf.gz", "_snv.vcf.gz", vcf.paths[i])
    
    ### read vcf file
    vcf <- fread(vcf.snvpath, data.table=FALSE, skip="#CHROM")

    ### Combine all SNVs into single dataframe with following columns:
    ### chromosome, position, reference base, alternative base, total read depth in normal sample excl filtered reads,
    ### total read depth in tumor sample excl filtered reads, read depth of reference base in normal sample
    ### read depth of reference base in tumor sample, read depth of alternative base in normal sample, 
    ### read depth of alternative base in tumor sample,
    ### allele frequency of alternative base in normal sample, allele frequency of alternative base in tumor sample,
    ### FILTER, average median mapping quality of reads in normal and tumor samples

    snv.df <- data.frame(chr = vcf$'#CHROM', pos = vcf$POS, ref = vcf$REF, alt = vcf$ALT, filter = vcf$FILTER)

    ### add average median mapping quality 
    MQ <- (as.numeric(system(paste("bcftools query -f '%INFO/MMQ\n'", vcf.snvpath, "| awk -F, '{print $1}'"), intern = TRUE)) +
        as.numeric(system(paste("bcftools query -f '%INFO/MMQ\n'", vcf.snvpath, "| awk -F, '{print $2}'"), intern = TRUE)))/2
    
    snv.df <- cbind(snv.df, mq = as.numeric(MQ))

    ### remove MQ variable
    rm(MQ)

    ### estimate number of reads supporting reference and alternative base in tumor and normal samples
    norm_ref_depth <- as.numeric(system(paste("bcftools query -s ", normal.id, " -f '[%AD]\n'", vcf.snvpath, "| awk -F, '{print $1}'"), intern = TRUE))
    norm_alt_depth <-  as.numeric(system(paste("bcftools query -s ", normal.id, " -f '[%AD]\n'", vcf.snvpath, "| awk -F, '{print $2}'"), intern = TRUE))
    tumor_ref_depth <-  as.numeric(system(paste("bcftools query -s ", tumor.id, " -f '[%AD]\n'", vcf.snvpath, "| awk -F, '{print $1}'"), intern = TRUE))
    tumor_alt_depth <- as.numeric(system(paste("bcftools query -s ", tumor.id, " -f '[%AD]\n'", vcf.snvpath, "| awk -F, '{print $2}'"), intern = TRUE))

    ### add depths and allele frequencies
    snv.df <- cbind(snv.df, totalnormal = norm_ref_depth+norm_alt_depth, normalrefdepth = norm_ref_depth, 
        normalaltdepth = norm_alt_depth, totaltumor = tumor_ref_depth + tumor_alt_depth, 
        tumorrefdepth = tumor_ref_depth, tumoraltdepth = tumor_alt_depth, afaltnormal = norm_alt_depth/(norm_alt_depth + norm_ref_depth),
        afalttumor = tumor_alt_depth/(tumor_alt_depth + tumor_ref_depth))
    snv.df$totaldepth <- snv.df$totalnormal+snv.df$totaltumor

     ### exclude multiallelic variants
    snv.df <- snv.df[which(!grepl("multiallelic", snv.df$filter)) ,]

    ### set id for sample
    id <- tail(strsplit(dirname(vcf.snvpath), split="/")[[1]],1)[1]

    ### make plots: 
    #### log10 distribution of totaldepth with cutoff at 40
    hist_plot <- hist(log10(snv.df$totaldepth), xlim = c(0, log10(max(snv.df$totaldepth)+1)), col = "skyblue", main = id, xlab = "Read depth (log10)", ylab = "Frequency")
    png(paste0("Histogram.allDP.mutect2.snvs.", id, ".png"), res=200, unit="in", height=8, width=11)
    plot(hist_plot, col = "skyblue", main = id, xlab = "Read depth (log10)", ylab = "Frequency")
    abline(v = log10(40), col = "red", lty = 2)
    dev.off()

    #### average mapping quality distribution with cutoff at 30
    hist_plot <- hist(snv.df$mq, xlim = c(0, max(snv.df$mq)), col = "skyblue", main = id, xlab = "Mapping quality", ylab = "Frequency")
    png(paste0("Histogram.MQ.mutect2.snvs.", id, ".png"), res=200, unit="in", height=8, width=11)
    plot(hist_plot, col = "skyblue", main = id, xlab = "Mapping quality", ylab = "Frequency")
    abline(v = 30, col = "red", lty = 2)
    dev.off()

    #### log10 distribution of totaltumor with cutoff at 25
    hist_plot <- hist(log10(snv.df$totaltumor), xlim = c(0, log10(max(snv.df$totaltumor)+1)), col = "skyblue", main = id, xlab = "Read depth in tumor (log10)", ylab = "Frequency")
    png(paste0("Histogram.TumorDP.mutect2.snvs.", id, ".png"), res=200, unit="in", height=8, width=11)
    plot(hist_plot, col = "skyblue", main = id, xlab = "Read depth in tumor (log10)", ylab = "Frequency")
    abline(v = log10(25), col = "red", lty = 2)
    dev.off()

    #### log10 distribution of totalnormal with cutoff at 15
    hist_plot <- hist(log10(snv.df$totalnormal), xlim = c(0, log10(max(snv.df$totalnormal)+1)), col = "skyblue", main = id, xlab = "Read depth in normal (log10)", ylab = "Frequency")
    png(paste0("Histogram.NormalDP.mutect2.snvs.", id, ".png"), res=200, unit="in", height=8, width=11)
    plot(hist_plot, col = "skyblue", main = id, xlab = "Read depth in normal (log10)", ylab = "Frequency")
    abline(v = log10(15), col = "red", lty = 2)
    dev.off()

    snv.sub.df <- snv.df[snv.df$totaltumor > 25 & snv.df$totalnormal > 15 & snv.df$mq > 30,]
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
        geom_vline(xintercept = 0.02, linetype = "dashed", color = "red") +
        geom_hline(yintercept = 0.05, linetype = "dashed", color = "red")
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
        geom_vline(xintercept = 0.02, linetype = "dashed", color = "red") +
        geom_hline(yintercept = 0.05, linetype = "dashed", color = "red")
        plot1})
    dev.off()

    #### afalttumor vs tumoraltdepth before read depth filter
    png(paste0("Scatterplot.AFAltTumorvsDPAltTumor.mutect2.snvs.beforeRDfilter.", id, ".png"), res=200, unit="in", height=8, width=11)
    print({plot1 <- ggplot(snv.df, aes(x=log10(tumoraltdepth+1), y=afalttumor)) +
        geom_point() +
        xlab("DPAltTumor (log10)") +
        ylab("AFAlt_in_tumor") +
        ylim(0,1) +
        geom_vline(xintercept = log10(5), linetype = "dashed", color = "red") +
        geom_hline(yintercept = 0.05, linetype = "dashed", color = "red")
        plot1})
    dev.off()
    
    #### afalttumor vs tumoraltdepth after read depth filter
    png(paste0("Scatterplot.AFAltTumorvsDPAltTumor.mutect2.snvs.afterRDfilter.", id, ".png"), res=200, unit="in", height=8, width=11)
    print({plot1 <- ggplot(snv.sub.df, aes(x=log10(tumoraltdepth+1), y=afalttumor)) +
        geom_point() +
        xlab("DPAltTumor (log10)") +
        ylab("AFAlt_in_tumor") +
        ylim(0,1) +
        geom_vline(xintercept = log10(5), linetype = "dashed", color = "red") +
        geom_hline(yintercept = 0.05, linetype = "dashed", color = "red")
        plot1})
    dev.off()

    #### afalttumor vs normalaltdepth before read depth filter
    png(paste0("Scatterplot.AFAltTumorvsDPAltNormal.mutect2.snvs.beforeRDfilter.", id, ".png"), res=200, unit="in", height=8, width=11)
    print({plot1 <- ggplot(snv.df, aes(x=log10(normalaltdepth+1), y=afalttumor)) +
        geom_point() +
        xlab("DPAltNormal (log10)") +
        ylab("AFAlt_in_tumor") +
        ylim(0,1) +
        geom_vline(xintercept = log10(2), linetype = "dashed", color = "red") +
        geom_hline(yintercept = 0.05, linetype = "dashed", color = "red")
        plot1})
    dev.off()
    
    #### afalttumor vs normalaltdepth after read depth filter
    png(paste0("Scatterplot.AFAltTumorvsDPAltNormal.mutect2.snvs.afterRDfilter.", id, ".png"), res=200, unit="in", height=8, width=11)
    print({plot1 <- ggplot(snv.sub.df, aes(x=log10(normalaltdepth+1), y=afalttumor)) +
        geom_point() +
        xlab("DPAltNormal (log10)") +
        ylab("AFAlt_in_tumor") +
        ylim(0,1) +
        geom_vline(xintercept = log10(2), linetype = "dashed", color = "red") +
        geom_hline(yintercept = 0.05, linetype = "dashed", color = "red")
        plot1})
    dev.off()

    snv.sub.df <- snv.sub.df[snv.sub.df$afaltnormal < 0.02 & snv.sub.df$afalttumor > 0.05 & snv.sub.df$tumoraltdepth >= 5 & snv.sub.df$normalaltdepth < 2,]
    cat("Number of variants left after read depth filtering ", nrow(snv.sub.df), "\n")

    snv.list[[i]] <- data.frame(sample = id, before = nrow(snv.df), after = nrow(snv.sub.df))

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
for (i in 1:length(vcf.paths)){
    cat("Processing", vcf.paths[i], "for INDELs\n")

    tumor.id <- ANNO_FILE$run_id[ANNO_FILE$is_tumor == 'Yes' & ANNO_FILE$vcf_path == basename(vcf.paths[i])]
    normal.id <- ANNO_FILE$run_id[ANNO_FILE$is_tumor == 'No' & ANNO_FILE$vcf_path == basename(vcf.paths[i])]

    vcf.indelpath <- gsub("_filtered.vcf.gz", "_indel.vcf.gz", vcf.paths[i])
    
    ### read vcf file
    vcf <- fread(vcf.indelpath, data.table=FALSE, skip="#CHROM")

    ### Combine all INDELs into single dataframe with following columns:
    ### chromosome, position, reference base, alternative base, total read depth in normal sample excl filtered reads,
    ### total read depth in tumor sample excl filtered reads, read depth of reference base in normal sample
    ### read depth of reference base in tumor sample, read depth of alternative base in normal sample, 
    ### read depth of alternative base in tumor sample,
    ### allele frequency of alternative base in normal sample, allele frequency of alternative base in tumor sample,
    ### FILTER, average median mapping quality of reads in normal and tumor samples

    indel.df <- data.frame(chr = vcf$'#CHROM', pos = vcf$POS, ref = vcf$REF, alt = vcf$ALT, filter = vcf$FILTER)

    ### add average median mapping quality 
    MQ <- (as.numeric(system(paste("bcftools query -f '%INFO/MMQ\n'", vcf.indelpath, "| awk -F, '{print $1}'"), intern = TRUE)) +
        as.numeric(system(paste("bcftools query -f '%INFO/MMQ\n'", vcf.indelpath, "| awk -F, '{print $2}'"), intern = TRUE)))/2
    
    indel.df <- cbind(indel.df, mq = as.numeric(MQ))

    ### remove MQ variable
    rm(MQ)

    ### estimate number of reads supporting reference and alternative base in tumor and normal samples
    norm_ref_depth <- as.numeric(system(paste("bcftools query -s ", normal.id, " -f '[%AD]\n'", vcf.indelpath, "| awk -F, '{print $1}'"), intern = TRUE))
    norm_alt_depth <-  as.numeric(system(paste("bcftools query -s ", normal.id, " -f '[%AD]\n'", vcf.indelpath, "| awk -F, '{print $2}'"), intern = TRUE))
    tumor_ref_depth <-  as.numeric(system(paste("bcftools query -s ", tumor.id, " -f '[%AD]\n'", vcf.indelpath, "| awk -F, '{print $1}'"), intern = TRUE))
    tumor_alt_depth <- as.numeric(system(paste("bcftools query -s ", tumor.id, " -f '[%AD]\n'", vcf.indelpath, "| awk -F, '{print $2}'"), intern = TRUE))

    ### add depths and allele frequencies
    indel.df <- cbind(indel.df, totalnormal = norm_ref_depth+norm_alt_depth, normalrefdepth = norm_ref_depth, 
        normalaltdepth = norm_alt_depth, totaltumor = tumor_ref_depth + tumor_alt_depth, 
        tumorrefdepth = tumor_ref_depth, tumoraltdepth = tumor_alt_depth, afaltnormal = norm_alt_depth/(norm_alt_depth + norm_ref_depth),
        afalttumor = tumor_alt_depth/(tumor_alt_depth + tumor_ref_depth))
    indel.df$totaldepth <- indel.df$totalnormal+indel.df$totaltumor

     ### exclude multiallelic variants
    indel.df <- indel.df[which(!grepl("multiallelic", indel.df$filter)) ,]

    ### set id for sample
    id <- tail(strsplit(dirname(vcf.indelpath), split="/")[[1]],1)[1]

    ### make plots: 
    #### log10 distribution of totaldepth with cutoff at 40
    hist_plot <- hist(log10(indel.df$totaldepth), xlim = c(0, log10(max(indel.df$totaldepth)+1)), col = "skyblue", main = id, xlab = "Read depth (log10)", ylab = "Frequency")
    png(paste0("Histogram.allDP.mutect2.indels.", id, ".png"), res=200, unit="in", height=8, width=11)
    plot(hist_plot, col = "skyblue", main = id, xlab = "Read depth (log10)", ylab = "Frequency")
    abline(v = log10(40), col = "red", lty = 2)
    dev.off()

    #### average mapping quality distribution with cutoff at 30
    hist_plot <- hist(indel.df$mq, xlim = c(0, max(indel.df$mq)), col = "skyblue", main = id, xlab = "Mapping quality", ylab = "Frequency")
    png(paste0("Histogram.MQ.mutect2.indels.", id, ".png"), res=200, unit="in", height=8, width=11)
    plot(hist_plot, col = "skyblue", main = id, xlab = "Mapping quality", ylab = "Frequency")
    abline(v = 50, col = "red", lty = 2)
    dev.off()

    #### log10 distribution of totaltumor with cutoff at 25
    hist_plot <- hist(log10(indel.df$totaltumor), xlim = c(0, log10(max(indel.df$totaltumor)+1)), col = "skyblue", main = id, xlab = "Read depth in tumor (log10)", ylab = "Frequency")
    png(paste0("Histogram.TumorDP.mutect2.indels.", id, ".png"), res=200, unit="in", height=8, width=11)
    plot(hist_plot, col = "skyblue", main = id, xlab = "Read depth in tumor (log10)", ylab = "Frequency")
    abline(v = log10(25), col = "red", lty = 2)
    dev.off()

    #### log10 distribution of totalnormal with cutoff at 15
    hist_plot <- hist(log10(indel.df$totalnormal), xlim = c(0, log10(max(indel.df$totalnormal)+1)), col = "skyblue", main = id, xlab = "Read depth in normal (log10)", ylab = "Frequency")
    png(paste0("Histogram.NormalDP.mutect2.indels.", id, ".png"), res=200, unit="in", height=8, width=11)
    plot(hist_plot, col = "skyblue", main = id, xlab = "Read depth in normal (log10)", ylab = "Frequency")
    abline(v = log10(15), col = "red", lty = 2)
    dev.off()

    indel.sub.df <- indel.df[indel.df$totaltumor > 25 & indel.df$totalnormal > 15 & indel.df$mq > 50,]
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
        geom_vline(xintercept = 0.02, linetype = "dashed", color = "red") +
        geom_hline(yintercept = 0.05, linetype = "dashed", color = "red")
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
        geom_vline(xintercept = 0.02, linetype = "dashed", color = "red") +
        geom_hline(yintercept = 0.05, linetype = "dashed", color = "red")
        plot1})
    dev.off()

    #### afalttumor vs tumoraltdepth before read depth filter
    png(paste0("Scatterplot.AFAltTumorvsDPAltTumor.mutect2.indels.beforeRDfilter.", id, ".png"), res=200, unit="in", height=8, width=11)
    print({plot1 <- ggplot(indel.df, aes(x=log10(tumoraltdepth+1), y=afalttumor)) +
        geom_point() +
        xlab("DPAltTumor (log10)") +
        ylab("AFAlt_in_tumor") +
        ylim(0,1) +
        geom_vline(xintercept = log10(5), linetype = "dashed", color = "red") +
        geom_hline(yintercept = 0.05, linetype = "dashed", color = "red")
        plot1})
    dev.off()
    
    #### afalttumor vs tumoraltdepth after read depth filter
    png(paste0("Scatterplot.AFAltTumorvsDPAltTumor.mutect2.indels.afterRDfilter.", id, ".png"), res=200, unit="in", height=8, width=11)
    print({plot1 <- ggplot(indel.sub.df, aes(x=log10(tumoraltdepth+1), y=afalttumor)) +
        geom_point() +
        xlab("DPAltTumor (log10)") +
        ylab("AFAlt_in_tumor") +
        ylim(0,1) +
        geom_vline(xintercept = log10(5), linetype = "dashed", color = "red") +
        geom_hline(yintercept = 0.05, linetype = "dashed", color = "red")
        plot1})
    dev.off()

    #### afalttumor vs normalaltdepth before read depth filter
    png(paste0("Scatterplot.AFAltTumorvsDPAltNormal.mutect2.indels.beforeRDfilter.", id, ".png"), res=200, unit="in", height=8, width=11)
    print({plot1 <- ggplot(indel.df, aes(x=log10(normalaltdepth+1), y=afalttumor)) +
        geom_point() +
        xlab("DPAltNormal (log10)") +
        ylab("AFAlt_in_tumor") +
        ylim(0,1) +
        geom_vline(xintercept = log10(2), linetype = "dashed", color = "red") +
        geom_hline(yintercept = 0.05, linetype = "dashed", color = "red")
        plot1})
    dev.off()
    
    #### afalttumor vs normalaltdepth after read depth filter
    png(paste0("Scatterplot.AFAltTumorvsDPAltNormal.mutect2.indels.afterRDfilter.", id, ".png"), res=200, unit="in", height=8, width=11)
    print({plot1 <- ggplot(indel.sub.df, aes(x=log10(normalaltdepth+1), y=afalttumor)) +
        geom_point() +
        xlab("DPAltNormal (log10)") +
        ylab("AFAlt_in_tumor") +
        ylim(0,1) +
        geom_vline(xintercept = log10(2), linetype = "dashed", color = "red") +
        geom_hline(yintercept = 0.05, linetype = "dashed", color = "red")
        plot1})
    dev.off()

    indel.sub.df <- indel.sub.df[indel.sub.df$afaltnormal < 0.02 & indel.sub.df$afalttumor > 0.05 & indel.sub.df$tumoraltdepth >= 5 & indel.sub.df$normalaltdepth < 2,]
    cat("Number of variants left after read depth filtering ", nrow(indel.sub.df), "\n")

    indel.list[[i]] <- data.frame(sample = id, before = nrow(indel.df), after = nrow(indel.sub.df))
    cat("Finished processing sample ", id, "\n")
}

### save counts of INDEL variants
write.table(do.call(rbind, indel.list), file = "Nmut@indel.mutect2.txt", sep = "\t", row.names = FALSE)

end_time <- Sys.time()
cat("Elapsed minutes for processing INDELs", end_time - start_time)