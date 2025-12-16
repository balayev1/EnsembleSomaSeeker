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
        config['mutect2_cores']
    input:
        ref = config["reference_genome"],
        tumor_bam = config["sample_data"]["tumor_bam_path"],
        normal_bam = config["sample_data"]["normal_bam_path"],
    output:
        vcf = f"{config['outputdir']}/mutect2/{config['subject_id']}_unfiltered.vcf.gz",
        f1r2 = f"{config['outputdir']}/mutect2/{config['subject_id']}.f1r2.tar.gz",
    params:
        tmpdir = config.get("tmpdir"),
        tumor_id = config["sample_data"]["tumor_sample_id"],
        normal_id = config["sample_data"]["normal_sample_id"],
        optional_args = (
            (f" --intervals {config['mutect2_params']['intervals_bed']}" 
             if config['mutect2_params'].get('intervals_bed') else "") +
            (f" --germline-resource {config['mutect2_params']['germline_resource']}" 
             if config['mutect2_params'].get('germline_resource') else "") +
            (f" --panel-of-normals {config['mutect2_params']['panel_of_normals']}" 
             if config['mutect2_params'].get('panel_of_normals') else "") +
            (" --mitochondria-mode" 
             if config['mutect2_params'].get('mitochondria_mode', 'false') == 'true' else "")
        )
    shell:
        """
        JAVA_OPTS="-Djava.io.tmpdir={params.tmpdir}"
        
        gatk Mutect2 \
            -R {input.ref} \
            -I {input.tumor_bam} -tumor {params.tumor_id} \
            -I {input.normal_bam} -normal {params.normal_id} \
            -O {output.vcf} \
            --stats {output.stats} \
            --native-pair-hmm-threads {threads} \

            {params.optional_args}
        """
    