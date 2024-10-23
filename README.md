# SV homology
## Overview

This repository contains code to find and annotate homologous flanks of structural variants (SV). Currently only deletions and insertion are supported. This tool was used to detect the homology landscape of [1019 samples from the 1000 Genomes Project sequenced with Nanopore](https://www.biorxiv.org/content/10.1101/2024.04.18.590093v1), but can be applied to other SV callset.

### System requirement
- blast (tested on ?)
- Python (tested on ?)
- DELLY (tested on ?)
- pysam (tested on ?)
- cyvcf2 (tested on ?)
- RepeatMasker (tested on ?)

### File requirement
- The vcf file is read in using cyvcf2 and the tool relies on following fields next to mandatory vcf fields:
  - INFO field - SVTYPE - DEL/INS
  - INFO field - SVLEN - Length of the SV in base pairs
- Genome reference as fasta file
- [RepeatMasker](https://github.com/Dfam-consortium/RepeatMasker) bed file for the reference genome

### Usage
Detecting homologous DNA streches flanking SV
```bash
python homology_detection.py \
  --vcf_file path/to/file.vcf \ #Path to the VCF file to be processed
  --top_workdir path/to/workdir \ #Path to the working directory
  --reference_fasta path/to/fasta.fa \ #Path to the reference FASTA file
  --blastn_bin path/to/blastn \ #Path to the binary of blastn
  --num_parallel X \ #Number of parallel processes
  --max_homology_length X #Maximium length of size of search window
```
