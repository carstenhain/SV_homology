# SV homology
## Overview

This projects finds homologous regions at the breakpoints of structural variants (SV). Those homologous sequences hint to a non-allelic homologous recombination (NAHR) mechanism of formation. Currently the tool words for deletions, inversions and duplications. Support for insertions was present in a earlier version and will be added soon.
This tool (in a earlier version but with the same logic) was used to detect the homology landscape of [1019 samples from the 1000 Genomes Project sequenced with Nanopore](https://www.nature.com/articles/s41586-025-09290-7), but can be applied to other SV calls to test for homologous flanks and a potential NAHR mechanism.

### System requirement
- blast (tested on 2.14.1+)
- Python (tested on 3.9.18)
- pysam (tested on 0.21.0)
- cyvcf2 (tested on 0.30.28)
- nextflow (tested on  24.04.4.5917)
later additions will need 
- DELLY (earlier version tested on 0.9.1)
- RepeatMasker (earlier version tested on 4.1.2)

### File requirement
- The VCF file is read in using cyvcf2 and the tool relies on following fields next to mandatory VCF fields:
  - INFO field - SVTYPE - DEL/DUP/INV
  - INFO field - SVLEN - Length of the SV in base pairs
- Genome reference as fasta file

### Usage
Detecting homologous DNA streches flanking SV
```bash
nextflow run SV_homology/pipeline/homology_detection.nf \
--vcf_file_path path/to/file.vcf \ ## Path to VCF file containing SV calls
--outdir path/to/outdir \ ## path to output folder
--final_name name \ ## name of the callset
--blastn_path path/to/blastn \ ## Path to blastn binary
--reference path/to/fasta.fa \ ## Path to the reference FASTA file
--max_homology_length L \ ## Maximum length of the search window 
--num_parallel = 10 ## Number of parallel processes
```
This tool outputs an annotated VCF file and different bed file indicating the homologus region into 'outdir'
```bash
*.homology.vcf                # Annotated VCF file
final_homology.bed            # BED file containing the homologous flanks of SV on the reference genome
final_candidate_homology.bed  # BED file containing candidate homologous flanks of the SV that did not pass all filtering criteria, they might be useful in cases of inaccurate SV breakpoint calling
```
