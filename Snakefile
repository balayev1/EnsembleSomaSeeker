"""
Snakefile to identify somatic SNVs and INDELs in single subject/patient in tumor-normal pair mode
"""

# Validation of config file
# from snakemake.utils import validate

# validate(config, schema="schemas/config_schema.yaml")

# MUTECT2 RULE
rule mutect2:
    conda:
        "envs/mutect2_ensemblesomaseeker.yaml"
    threads:
        config['rule_cores'].get('mutect2', 2)
    input:
        ref = config["reference_genome"],
        tumor_bam = config["sample_data"]["tumor_bam_path"],
        normal_bam = config["sample_data"]["normal_bam_path"],
    output:
        unfiltered_vcf = f"{config['outputdir']}/mutect2/{config['subject_id']}_unfiltered.vcf.gz",
        f1r2 = f"{config['outputdir']}/mutect2/{config['subject_id']}.f1r2.tar.gz",
    params:
        # Move directory logic here
        out_dir = f"{config['outputdir']}/mutect2",
        tmpdir = config.get("tmpdir"),
        tumor_id = config["sample_data"]["tumor_sample_id"],
        normal_id = config["sample_data"]["normal_sample_id"],
        optional_args = (
            (f" --intervals {config.get('intervals_bed')} --interval-padding 0" 
             if config.get('intervals_bed') else "") +
            (f" --germline-resource {config['mutect2_params'].get('germline_resource')}" 
             if config['mutect2_params'].get('germline_resource') else "") +
            (f" --panel-of-normals {config['mutect2_params'].get('panel_of_normals')}" 
             if config['mutect2_params'].get('panel_of_normals') else "") +
            (f" --tmp-dir {config.get('tmpdir')}"
             if config.get('tmpdir') else "")
        )
    shell:
        """        
        mkdir -p {params.out_dir}

        gatk Mutect2 \
            -R {input.ref} \
            -I {input.tumor_bam} -tumor {params.tumor_id} \
            -I {input.normal_bam} -normal {params.normal_id} \
            -O {output.unfiltered_vcf} \
            --f1r2-tar-gz {output.f1r2} \
            --native-pair-hmm-threads {threads} \
            --af-of-alleles-not-in-resource 0.0000025 \
            {params.optional_args}
        """

# FILTER MUTECT2 RULE
rule filter_mutect2:
    conda:
        "envs/mutect2_ensemblesomaseeker.yaml"
    input:
        ref = config["reference_genome"],
        tumor_bam = config["sample_data"]["tumor_bam_path"],
        normal_bam = config["sample_data"]["normal_bam_path"],
        unfiltered_vcf = f"{config['outputdir']}/mutect2/{config['subject_id']}_unfiltered.vcf.gz",
        biallelic_resource = config['filter_mutect2_params']['biallelic_resource'],
        f1r2 = f"{config['outputdir']}/mutect2/{config['subject_id']}.f1r2.tar.gz",
    output:
        tumor_pileup = temp(f"{config['outputdir']}/mutect2_filtered/{config['subject_id']}_tumor.pileups.table"),
        normal_pileup = temp(f"{config['outputdir']}/mutect2_filtered/{config['subject_id']}_normal.pileups.table"),
        contamination_table = temp(f"{config['outputdir']}/mutect2_filtered/{config['subject_id']}_contamination.table"),
        segmentation_table = temp(f"{config['outputdir']}/mutect2_filtered/{config['subject_id']}_segments.table"),
        read_orientation_table = temp(f"{config['outputdir']}/mutect2_filtered/{config['subject_id']}.artifact_prior_tables.tar.gz"),
        filtered_vcf = f"{config['outputdir']}/mutect2_filtered/{config['subject_id']}_filtered.vcf.gz",
    params:
        # Define the directory here to avoid the quote error in shell
        out_dir = f"{config['outputdir']}/mutect2_filtered"
    shell:
        """        
        mkdir -p {params.out_dir}

        gatk GetPileupSummaries \
            -I {input.tumor_bam} \
            -L {input.biallelic_resource} \
            -V {input.biallelic_resource} \
            -O {output.tumor_pileup}
        
        gatk GetPileupSummaries \
            -I {input.normal_bam} \
            -L {input.biallelic_resource} \
            -V {input.biallelic_resource} \
            -O {output.normal_pileup}
        
        gatk CalculateContamination \
            -I {output.tumor_pileup} \
            --matched-normal {output.normal_pileup} \
            -O {output.contamination_table} \
            --tumor-segmentation {output.segmentation_table}
        
        gatk LearnReadOrientationModel \
            -I {input.f1r2} \
            -O {output.read_orientation_table}
        
        gatk FilterMutectCalls \
            -R {input.ref} \
            -V {input.unfiltered_vcf} \
            --orientation-bias-artifact-priors {output.read_orientation_table} \
            --contamination-table {output.contamination_table} \
            --tumor-segmentation {output.segmentation_table} \
            -O {output.filtered_vcf}
        """

# PREPARE BED FOR STRELKA2 RULE 
rule prepare_strelka_bed:
    conda:
        "envs/bcftools_ensemblesomaseeker.yaml"
    input:
        bed = config.get("intervals_bed", "")
    output:
        gz = f"{config['outputdir']}/strelka2/intervals.bed.gz",
        tbi = f"{config['outputdir']}/strelka2/intervals.bed.gz.tbi"
    shell:
        """
        mkdir -p {config['outputdir']}/strelka2

        if [[ "{input.bed}" == *.gz ]]; then
            cp {input.bed} {output.gz}
        else
            bgzip -c {input.bed} > {output.gz}
        fi
        
        tabix -f -p bed {output.gz}
        """

# STRELKA2 RULE
rule strelka2:
    conda:
        "envs/strelka2_ensemblesomaseeker.yaml"
    threads:
        config['rule_cores'].get('strelka2', 2)
    input:
        ref = config["reference_genome"],
        tumor_bam = config["sample_data"]["tumor_bam_path"],
        normal_bam = config["sample_data"]["normal_bam_path"],
        prepared_bed = (f"{config['outputdir']}/strelka2/intervals.bed.gz" 
                        if config.get("intervals_bed") else []),
    output:
        snvs = f"{config['outputdir']}/strelka2/{config['subject_id']}/results/variants/somatic.snvs.vcf.gz",
        indels = f"{config['outputdir']}/strelka2/{config['subject_id']}/results/variants/somatic.indels.vcf.gz",
        run_dir = directory(f"{config['outputdir']}/strelka2/{config['subject_id']}"),
    params:
        optional_args = (
            (f" --callRegions {config['outputdir']}/strelka2/intervals.bed.gz" 
                        if config.get("intervals_bed") else "")
        )
    shell:
        """        
        rm -rf {output.run_dir}

        mkdir -p {output.run_dir}

        python2 configureStrelkaSomaticWorkflow.py \
            --tumorBam {input.tumor_bam} \
            --normalBam {input.normal_bam} \
            --referenceFasta {input.ref} \
            --runDir {output.run_dir} \
            {params.optional_args}

        python2 {output.run_dir}/runWorkflow.py \
            --mode local \
            --jobs {threads} \
            --memGb 64
        """

# VARSCAN2 RULE
rule varscan2:
    conda:
        "envs/varscan2_ensemblesomaseeker.yaml"
    threads:
        config['rule_cores'].get('varscan2', 2)
    input:
        ref = config["reference_genome"],
        tumor_bam = config["sample_data"]["tumor_bam_path"],
        normal_bam = config["sample_data"]["normal_bam_path"],
    output:
        snvs = f"{config['outputdir']}/varscan2/{config['subject_id']}.snp.vcf",
        indels = f"{config['outputdir']}/varscan2/{config['subject_id']}.indel.vcf",
    params:
        # Move directory and prefix here
        out_dir = f"{config['outputdir']}/varscan2",
        prefix = f"{config['outputdir']}/varscan2/{config['subject_id']}",
        optional_args = (f" -l {config.get('intervals_bed')}" if config.get('intervals_bed') else "")
    shell:
        """
        mkdir -p {params.out_dir}

        samtools mpileup \
            --no-BAQ \
            --min-MQ 1 \
            -f {input.ref} \
            {params.optional_args} \
            {input.normal_bam} \
            {input.tumor_bam} | \
        varscan somatic \
            - \
            {params.prefix} \
            --mpileup 1 \
            --min-coverage 1 \
            --min-var-freq 0.001 \
            --output-vcf 1
        """

# MUSE RULE
checkpoint muse:
    conda:
        "envs/muse_ensemblesomaseeker.yaml"
    threads:
        config['rule_cores'].get('muse', 2)
    input:
        ref = config["reference_genome"],
        tumor_bam = config["sample_data"]["tumor_bam_path"],
        normal_bam = config["sample_data"]["normal_bam_path"],
    output:
        target_vcf = f"{config['outputdir']}/muse/{config['subject_id']}.MuSE.vcf"
    params:
        out_dir = f"{config['outputdir']}/muse",
        prefix = f"{config['outputdir']}/muse/{config['subject_id']}",
        dbsnp = (f" -D {config.get('dbsnp_resource')}" 
                 if config.get('dbsnp_resource') else "")
    shell:
        """
        mkdir -p {params.out_dir}

        (
            {params.muse_bin} call \
                -f {input.ref} \
                -O {params.prefix} \
                -n {threads} \
                {input.tumor_bam} \
                {input.normal_bam} && \
            
            {params.muse_bin} sump \
                -I {params.prefix}.MuSE.txt \
                -O {output.target_vcf} \
                -G \
                -n {threads} \
                {params.dbsnp}
        ) || (echo "MuSE failed, creating empty recovery file" && touch {output.target_vcf})
        """

# LOFREQ RULE
rule lofreq:
    conda:
        "envs/lofreq_ensemblesomaseeker.yaml"
    threads:
        config['rule_cores'].get('lofreq', 2)
    input: 
        ref = config["reference_genome"],
        tumor_bam = config["sample_data"]["tumor_bam_path"],
        normal_bam = config["sample_data"]["normal_bam_path"],
    output:
        realign_tumor = temp(f"{config['outputdir']}/lofreq/{config['subject_id']}_tumor_BD_BI.bam"),
        realign_normal = temp(f"{config['outputdir']}/lofreq/{config['subject_id']}_normal_BD_BI.bam"),
        realign_tumor_bai = temp(f"{config['outputdir']}/lofreq/{config['subject_id']}_tumor_BD_BI.bam.bai"),
        realign_normal_bai = temp(f"{config['outputdir']}/lofreq/{config['subject_id']}_normal_BD_BI.bam.bai"),
        snvs = f"{config['outputdir']}/lofreq/{config['subject_id']}_somatic_final.snvs.vcf.gz",
        indels = f"{config['outputdir']}/lofreq/{config['subject_id']}_somatic_final.indels.vcf.gz",
    params:
        # Move directory logic here
        out_dir = f"{config['outputdir']}/lofreq",
        prefix = f"{config['outputdir']}/lofreq/{config['subject_id']}_",
        optional_args = (
            (f" --bed {config.get('intervals_bed')}" 
            if config.get('intervals_bed') else "") +
            (f" --dbsnp {config.get('dbsnp_resource')}" 
            if config.get('dbsnp_resource') else "")
        )
    shell:
        """        
        mkdir -p {params.out_dir}

        lofreq indelqual \
            --dindel \
            --ref {input.ref} \
            --out {output.realign_normal} \
            {input.normal_bam}
        
        lofreq indelqual \
            --dindel \
            --ref {input.ref} \
            --out {output.realign_tumor} \
            {input.tumor_bam}
        
        samtools index -o {output.realign_tumor_bai} {output.realign_tumor}

        samtools index -o {output.realign_normal_bai} {output.realign_normal}

        lofreq somatic \
            --normal {output.realign_normal} \
            --tumor {output.realign_tumor} \
            --outprefix {params.prefix} \
            --ref {input.ref} \
            --threads {threads} \
            --call-indels \
            --min-cov 7 \
            {params.optional_args}
        """

# Define input VCFs for merging in SomaticSeq
def get_somaticseq_inputs(wildcards):
    muse_vcf = checkpoints.muse.get(**wildcards).output.target_vcf
    
    inputs = {
        "mutect2_vcf": f"{config['outputdir']}/mutect2_filtered/{wildcards.subject_id}_filtered.vcf.gz",
        "strelka2_snvs": f"{config['outputdir']}/strelka2/{wildcards.subject_id}/results/variants/somatic.snvs.vcf.gz",
        "strelka2_indels": f"{config['outputdir']}/strelka2/{wildcards.subject_id}/results/variants/somatic.indels.vcf.gz",
        "varscan2_snvs": f"{config['outputdir']}/varscan2/{wildcards.subject_id}.snp.vcf",
        "varscan2_indels": f"{config['outputdir']}/varscan2/{wildcards.subject_id}.indel.vcf",
        "lofreq_snvs": f"{config['outputdir']}/lofreq/{wildcards.subject_id}_somatic_final.snvs.vcf.gz",
        "lofreq_indels": f"{config['outputdir']}/lofreq/{wildcards.subject_id}_somatic_final.indels.vcf.gz",
    }
    
    import os
    if os.path.exists(muse_vcf) and os.path.getsize(muse_vcf) > 0:
        inputs["muse_vcf"] = muse_vcf
    else:
        inputs["muse_vcf"] = None 
        
    return inputs

# SOMATICSEQ RULE
rule somaticseq:
    conda:
        "envs/somaticseq_ensemblesomaseeker.yaml"
    threads:
        config['rule_cores'].get('somaticseq', 2)
    input: 
        unpack(get_somaticseq_inputs),
        ref = config["reference_genome"],
        tumor_bam = config["sample_data"]["tumor_bam_path"],
        normal_bam = config["sample_data"]["normal_bam_path"],
    output:
        snvs_vcf = f"{config['outputdir']}/somaticseq/{config['subject_id']}/Consensus.sSNV.vcf",
        indels_vcf = f"{config['outputdir']}/somaticseq/{config['subject_id']}/Consensus.sINDEL.vcf",
        snvs_tsv = f"{config['outputdir']}/somaticseq/{config['subject_id']}/Ensemble.sSNV.tsv",
        indels_tsv = f"{config['outputdir']}/somaticseq/{config['subject_id']}/Ensemble.sINDEL.tsv",
    params:
        prefix = f"{config['outputdir']}/somaticseq/{config['subject_id']}",
        muse_arg = lambda wildcards, input: f"--muse-vcf {input.muse_vcf}" if hasattr(input, 'muse_vcf') and input.muse_vcf else ""
        optional_args = (
            (f" --inclusion-region {config.get('intervals_bed')}" 
            if config.get('intervals_bed') else "") +
            (f" --dbsnp-vcf {config.get('dbsnp_resource')}" 
            if config.get('dbsnp_resource') else "") +
            (f" --cosmic-vcf {config.get('cosmic_resource')}" 
            if config.get('cosmic_resource') else "")
        )
    shell:
        """
        mkdir -p {config['outputdir']}/somaticseq/{config['subject_id']}

        somaticseq_parallel.py paired \
            --tumor-bam-file {input.tumor_bam} \
            --normal-bam-file {input.normal_bam} \
            --mutect2-vcf {input.mutect2_vcf} \
            --varscan-snv {input.varscan2_snvs} \
            --varscan-indel {input.varscan2_indels} \
            {params.muse_arg} \
            --lofreq-snv {input.lofreq_snvs} \
            --lofreq-indel {input.lofreq_indels} \
            --strelka-snv {input.strelka2_snvs} \
            --strelka-indel {input.strelka2_indels} \
            --output-directory {params.prefix} \
            --genome-reference {input.ref} \
            --threads {threads} \
            --pass-threshold 0.5 \
            --lowqual-threshold 0.1 \
            {params.optional_args}
        """