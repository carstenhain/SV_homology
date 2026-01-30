#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.vcf_file_path = "/home/hain/EMBL/inversions/homology/final_inv_calls.sort.vcf"
params.outdir = "/home/hain/EMBL/inversions/homology/"
params.final_name = "final_inv_calls_homology.vcf"

params.blastn_path = "/home/hain/miniconda3/envs/jupyter_wave/bin/blastn"
params.reference = "/home/hain/EMBL/References/chm13v2.0.fa"
params.max_homology_length = 10000
params.num_parallel = 10

/*
 *nextflow run SV_homology/pipeline/homology_detection.nf \
 * --vcf_file_path /home/hain/EMBL/inversions/homology/final_inv_calls.vcf \
 * --outdir /home/hain/EMBL/inversions/homology/ \
 * --final_name final_inv_calls_homology.vcf \
 * --num_parallel 10
 */


process split_vcf {

    input:
    path vcf_file
    val num_parts
    path split_script

    output:
    path "splits/*.vcf", emit: split_vcfs

    script:
    """
    mkdir -p splits

    python ${split_script} \
        --vcf_file ${vcf_file} \
        --workdir splits \
        --num_splits ${num_parts}
    """
}


process detect_homology {

    input:
    path vcf_part
    path homology_script

    output:
    path "modified/*.homology.vcf.gz", emit: homology_vcf
    path "modified/*.homology.vcf.gz.csi", emit: homology_vcf_index
    path "modified/*.homology.bed", emit: homology_bed
    path "modified/*.candidate_homology.bed", emit: candidate_homology_bed

    script:
    """
    mkdir -p modified

    python ${homology_script} \
        --vcf_file ${vcf_part} \
        --workdir modified \
        --reference_fasta ${params.reference} \
        --blastn_bin ${params.blastn_path} \
        --max_homology_length ${params.max_homology_length}

    # sort, compresses and index the homology VCF
    for vcf in modified/*.homology.vcf; do
        bcftools sort \$vcf -O z -o \$vcf.gz
        bcftools index \$vcf.gz
    done
    """
}


process concat_outputs {

    publishDir params.outdir, mode: 'copy'

    input:
    path homology_vcfs
    path homology_vcf_indices
    path homology_beds
    path candidate_beds

    output:
    path params.final_name
    path "final_homology.bed"
    path "final_candidate_homology.bed"

    script:
    """
    # concatenate VCFs
    bcftools concat -a ${homology_vcfs.sort()} -O z -o ${params.final_name}

    # concatenate BED files
    cat ${homology_beds.sort()} > final_homology.bed
    cat ${candidate_beds.sort()} > final_candidate_homology.bed
    """
}

workflow {

    if (!params.vcf_file_path) error "Missing required parameter: vcf_file_path"
    if (!params.outdir) error "Missing required parameter: outdir"
    if (!params.blastn_path) error "Missing required parameter: blastn_path"
    if (!params.reference) error "Missing required parameter: reference"

    def repoRoot = projectDir.parent
    def split_vcf_script = "${repoRoot}/detection/split_vcf.py"
    def homology_detection_script = "${repoRoot}/detection/homology_detection.py"

    vcf_ch = channel.fromPath(params.vcf_file_path, checkIfExists: true)
    split_script_ch = channel.fromPath(split_vcf_script, checkIfExists: true)
    //homology_script_ch = channel.fromPath(homology_detection_script, checkIfExists: true)

    split_vcfs_ch = split_vcf(vcf_ch, params.num_parallel, split_script_ch)

    homology_results = detect_homology(
        split_vcfs_ch.flatten(),
        channel.value(homology_detection_script)
    )

    concat_outputs(
        homology_results.homology_vcf.collect(),
        homology_results.homology_vcf_index.collect(),
        homology_results.homology_bed.collect(),
        homology_results.candidate_homology_bed.collect()
    )
}