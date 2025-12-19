# EnsembleSomaSeeker

EnsembleSomaSeeker is a somatic mutation calling pipeline that uses BAM/SAM/CRAM files containing reads in tumor-normal pair mode, and identifies somatic mutations using two-step approach. In step I, it takes the reads from tumor-normal sample pair and finds somatic single-nucleotide variants (SNV) and short insertions/deletions (INDEL) using five different mutation-calling software: Mutect2 from GATK package, Strelka2, Varscan2, MuSE and Lofreq. Somatic mutations are then collapsed together using SomaticSeq package. In step 2, EnsembleSomaSeeker selects mutations called by at least two different software and passes them through several hard filtering techniques resulting in the final VCF file per tumor-normal pair. Pipeline can be run on any whole-genome (WGS) sequencing tumor-normal sample pair in BAM/SAM/CRAM format and is flexible in terms of parameter choice selection for each mutation calling software or filtering options.

## Setup
Clone repository:
```bash 
git clone git@github.com:balayev1/EnsembleSomaSeeker.git 
```
Change path to EnsembleSomaSeeker repository:
```bash
cd EnsembleSomaSeeker
```

## Usage
EnsembleSomaSeeker can now be executed using simple command in your terminal:
```bash
bash EnsembleSomaSeeker.sh \
    -t [path/to/tumor.bam] \
    -n [path/to/normal.bam] \
    -T [tumor_id] \
    -N [normal_id] \
    -R [path/to/reference_fasta] \
    -o [output_dir] \
    -b [path/to/biallelic.vcf]
```

**Warning:** You must have dictionary and index files of reference fasta in `.dict` and `.fai` formats within same directory of your fasta file. VCF file with biallelic SNPs is also required for filtering part of Mutect2. Please check `data` folder to see contents of each file. 