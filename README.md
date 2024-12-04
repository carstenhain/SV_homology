# SV homology
## Overview

This repository contains code to find and annotate homologous flanks of structural variants (SV). Currently only deletions and insertion are supported. This tool was used to detect the homology landscape of [1019 samples from the 1000 Genomes Project sequenced with Nanopore](https://www.biorxiv.org/content/10.1101/2024.04.18.590093v1), but can be applied to other SV callset.

### System requirement
- blast (tested on 2.12.0+)
- Python (tested on 3.7.12)
- DELLY (tested on 0.9.1)
- pysam (tested on 0.21.0)
- cyvcf2 (tested on 0.30.28)
- RepeatMasker (tested on 4.1.2)

### File requirement
- The VCF file is read in using cyvcf2 and the tool relies on following fields next to mandatory VCF fields:
  - INFO field - SVTYPE - DEL/INS
  - INFO field - SVLEN - Length of the SV in base pairs
  - The VCF file must be gzipped and end on .vcf.gz
- Genome reference as fasta file
- [RepeatMasker](https://github.com/Dfam-consortium/RepeatMasker) bed file for the reference genome

### Usage
Detecting homologous DNA streches flanking SV
```bash
python homology_detection.py \
  --vcf_file path/to/file.vcf.gz \ #Path to the gzipped VCF file to be processed
  --top_workdir path/to/workdir \ #Path to the working directory
  --reference_fasta path/to/fasta.fa \ #Path to the reference FASTA file
  --blastn_bin path/to/blastn \ #Path to the binary of blastn
  --num_parallel X \ #Number of parallel processes
  --max_homology_length X #Maximium length of size of search window
```
The outout of this tool contains an annotated VCF file and different bed and fasta file indicating the homologus region 
```bash
*.homology.vcf.gz        # Annotated VCF file
*.del.homology.bed       # BED file indicating the positions of homologous flanks of deletions on the reference genome
*.ins_ref.fa             # "Pseudoassemblies" containing the sequence of an insertion flanked by locus specific reference sequence
*.ins_ref.bed            # BED file indicating the position of the inserted sequence on *.ins_ref.fa 
*.ins_ref.homology.bed    # BED file indicating the positions of homologous flanks of insertions on *.ins_ref.fa
```
