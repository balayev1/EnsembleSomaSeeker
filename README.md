# EnsembleSomaSeeker

EnsembleSomaSeeker is a somatic mutation calling pipeline that uses BAM/SAM/CRAM files containing reads in tumor-normal pair mode, and identifies somatic mutations using two-step approach. In step I, it takes the reads from tumor-normal sample pair and finds somatic single-nucleotide variants (SNV) and short insertions/deletions (INDEL) using five different mutation-calling software: Mutect2 from GATK package, Strelka2, Varscan2, MuSE and Lofreq. Somatic mutations are then collapsed together using SomaticSeq package. In step 2, EnsembleSomaSeeker selects mutations called by at least two different software and passes them through several hard filtering techniques resulting in the final VCF file per tumor-normal pair. Pipeline can be run on any whole-genome (WGS) sequencing tumor-normal sample pair in BAM/SAM/CRAM format and is flexible in terms of parameter choice selection for each mutation calling software or filtering options.

## Setup
Clone repository:
``` 
git clone git@github.com:balayev1/EnsembleSomaSeeker.git 
```
Change path to EnsembleSomaSeeker repository:
```
cd EnsembleSomaSeeker
```
To implement the test run, you have to download reference genome fasta [fasta file](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/GRCh38.primary_assembly.genome.fa.gz) and decompress it (3GB in size).
```
# Download
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/GRCh38.primary_assembly.genome.fa.gz -P data/

# Decompress the file
gzip -d data/GRCh38.primary_assembly.genome.fa.gz 
```

Now, we are ready for test run:
``` 
bash EnsembleSomaSeeker.sh \
    -t test/SRR14097759_sorted_MD_BQ_chr1_1000000.bam \
    -n test/SRR14097760_sorted_MD_BQ_chr1_1000000.bam \
    -R data/GRCh38.primary_assembly.genome.fa \
    -o test \
    -b data/common_all_biallelic.vcf.gz
```