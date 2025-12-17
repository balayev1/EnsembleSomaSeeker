"""
Snakefile to identify somatic SNVs and INDELs in single subject/patient in tumor-normal pair mode
"""

# Validation of config file
from snakemake.utils import validate

validate(config, schema="schemas/config_schema.yaml")

# MUTECT2 RULE
rule mutect2:
    conda:
        "envs/gatk_mutect2.yaml"
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
        tmpdir = config.get("tmpdir"),
        tumor_id = config["sample_data"]["tumor_sample_id"],
        normal_id = config["sample_data"]["normal_sample_id"],

        optional_args = (
            (f" --intervals {config['mutect2_params'].get('intervals_bed')} --interval-padding 0" 
             if config['mutect2_params'].get('intervals_bed') else "") +
            (f" --germline-resource {config['mutect2_params'].get('germline_resource')}" 
             if config['mutect2_params'].get('germline_resource') else "") +
            (f" --panel-of-normals {config['mutect2_params'].get('panel_of_normals')}" 
             if config['mutect2_params'].get('panel_of_normals') else "") +
            (f" --tmp-dir {config.get('tmpdir')}"
             if config.get('tmpdir') else "")
        )
    shell:
        """        
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
        "envs/gatk_mutect2.yaml"
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
    shell:
        """        
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