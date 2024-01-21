txt <- read.table("/scratch.global/balay011/WGS_skChordomas_Bai/scripts/WGS_chordomaBai_anno.txt", header=FALSE, row.names=1)

## number of subjects
n <- length(unique(txt$V3))

## change to common directory if caller was executed
mutect2 <- "/home/aventeic/balay011/WGS_skChordomas_Bai/SomaticSNV_calling/mutect2/samples"
varscan_snv <- "/home/aventeic/balay011/WGS_skChordomas_Bai/SomaticSNV_calling/varscan/samples"
varscan_indel <- "/home/aventeic/balay011/WGS_skChordomas_Bai/SomaticSNV_calling/varscan/samples"
jsm <- NULL
somaticsniper <- NULL
vardict_vcf <- NULL
muse <- "/home/aventeic/balay011/WGS_skChordomas_Bai/SomaticSNV_calling/muse"
lofreq_snv <- "/home/aventeic/balay011/WGS_skChordomas_Bai/SomaticSNV_calling/lofreq/samples"
lofreq_indel <- "/home/aventeic/balay011/WGS_skChordomas_Bai/SomaticSNV_calling/lofreq/samples"
scalpel_vcf <- NULL
strelka_snv <- "/home/aventeic/balay011/WGS_skChordomas_Bai/SomaticSNV_calling/strelka2/samples"
strelka_indel <- "/home/aventeic/balay011/WGS_skChordomas_Bai/SomaticSNV_calling/strelka2/samples"
arbitrary_snvs <- "/home/aventeic/balay011/WGS_skChordomas_Bai/SomaticSNV_calling/varscan/samples"
arbitrary_indels <- NULL

subj <- c(); tumor_bam <- c(); normal_bam <- c(); mutect2_vcf <- rep('NULL', n); varscan_snv_vcf <- rep('NULL', n); varscan_indel_vcf <- rep('NULL', n)
jsm_vcf <- rep('NULL', n); somaticsniper_vcf <- rep('NULL', n); vardict_vcf <- rep('NULL', n); muse_vcf <- rep('NULL', n); lofreq_snv_vcf <- rep('NULL', n)
lofreq_indel_vcf <- rep('NULL', n); scalpel_vcf <- rep('NULL', n); strelka_snv_vcf <- rep('NULL', n); strelka_indel_vcf <- rep('NULL', n)
arbitrary_snvs_vcf <- rep('NULL', n); arbitrary_indels_vcf <- rep('NULL', n)


for (i in 1:length(unique(txt$V3))){
    txt_sub <- txt[txt$V3 == unique(txt$V3)[i],]
    subj <- append(subj, unique(txt_sub$V3))
    tumor_bam <- append(tumor_bam, txt_sub$V4[txt_sub$V2 == "Yes"])
    normal_bam <- append(normal_bam, txt_sub$V4[txt_sub$V2 == "No"])

    if (!is.null(mutect2)){
        mutect2_vcf[i] <- file.path(mutect2, unique(txt_sub$V3), paste0(unique(txt_sub$V3), "_filtered.vcf.gz"))
    }

    if (!is.null(varscan_snv)){
        varscan_snv_vcf[i] <- file.path(varscan_snv, unique(txt_sub$V3), paste0(unique(txt_sub$V3), ".snp.vcf"))
    }

    if (!is.null(varscan_indel)){
        varscan_indel_vcf[i] <- file.path(varscan_indel, unique(txt_sub$V3), paste0(unique(txt_sub$V3), ".indel.vcf"))
    }

    if (!is.null(muse)){
        muse_vcf[i] <- file.path(muse, paste0(unique(txt_sub$V3), ".MuSE.vcf"))
    }

    if (!is.null(lofreq_snv)){
        lofreq_snv_vcf[i] <- file.path(lofreq_snv, unique(txt_sub$V3), "somatic_final.snvs.vcf.gz")
    }

    if (!is.null(lofreq_indel)){
        lofreq_indel_vcf[i] <- file.path(lofreq_indel, unique(txt_sub$V3), "somatic_final.indels.vcf.gz")
    }

    if (!is.null(strelka_snv)){
        strelka_snv_vcf[i] <- file.path(strelka_snv, unique(txt_sub$V3), "results/variants/somatic.snvs.vcf.gz")
    }

    if (!is.null(strelka_indel)){
        strelka_indel_vcf[i] <- file.path(strelka_indel, unique(txt_sub$V3), "results/variants/somatic.indels.vcf.gz")
    }
}

## Make single dataframe
mut_df <- data.frame(subj = subj, tumor_bam = tumor_bam, normal_bam = normal_bam, mutect2_vcf = mutect2_vcf, varscan_snv_vcf = varscan_snv_vcf, 
    varscan_indel_vcf = varscan_indel_vcf, jsm_vcf = jsm_vcf, somaticsniper_vcf = somaticsniper_vcf, vardict_vcf = vardict_vcf, muse_vcf = muse_vcf,
    lofreq_snv_vcf = lofreq_snv_vcf, lofreq_indel_vcf = lofreq_indel_vcf, scalpel_vcf = scalpel_vcf, strelka_snv_vcf = strelka_snv_vcf, 
    strelka_indel_vcf = strelka_indel_vcf, arbitrary_snvs_vcf = arbitrary_snvs_vcf, arbitrary_indels_vcf = arbitrary_indels_vcf)

## Save the dataframe
write.table(mut_df, file = "/home/aventeic/balay011/scripts/WGS_chordomaBai_somaticseq.txt", sep="\t", row.names=FALSE)