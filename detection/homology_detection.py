import multiprocessing
import os
import cyvcf2 # type: ignore
import datetime
from multiprocessing import Pool
import pysam # type: ignore
import subprocess
import numpy as np # type: ignore
from Bio.Seq import Seq # type: ignore
import argparse
import pandas as pd # type: ignore
from tqdm import tqdm # type: ignore

def make_reference_ins(record:cyvcf2.Variant, reference:pysam.FastaFile) -> tuple:
    """
    Build reference allele for an insertion
    Inserts inserted sequence between two flanks from the reference genome

    Args:
        record (cyvcf2.VCFRecord): insertion record, contains inserted sequence and insertion position
        reference (pysam.FastaFile): reference genome sequence

    Returns:
        tuple: sequence of reference allele with insertion, start and end position of the insertion
    """
    ### make reference with inserted sequence
    # from left to right 
    # reference genome: 2x SV length upstream of insertion
    # inserted sequence
    # reference genome: 2x SV length downstream of insertion
    ins_ref = reference.fetch(record.CHROM, record.POS - 2 * record.INFO["SVLEN"], record.POS - 1).rstrip() + \
        record.ALT[0].rstrip() + \
            reference.fetch(record.CHROM, record.POS, record.POS + 2 * record.INFO["SVLEN"]).rstrip()
    return ins_ref, 2 * record.INFO["SVLEN"], 3 * record.INFO["SVLEN"]

def make_intervals(sv_length:int, number_uneven_intervals:int) -> list:
    """
    Makes different sized intervals to use for the homology search
    Use only centered intervals for now --> add uncentered intervals later

    Args:
        sv_length (int): Length of the structural variant
        number_uneven_intervals (int): Number of uneven intervals to generate

    Returns:
        list: List of intervals for homology search
    """
    intervals = []
    # add centered intervals up to half the SV length
    for i in [10, 15, 20, 25, 30, 35, 40, 45, 50, 75, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 1000, 2000, 3000, 4000, 5000, 10000, 50000, 100000, int(sv_length / 2) - 1]:
        if i < sv_length / 2 and i > 0:
            intervals.append([-i, i])
    
    # add uncentered intervals with size = SV length
    for k in range(0, number_uneven_intervals + 1):
        intervals.append([-int(sv_length * k / number_uneven_intervals), int(sv_length * (1 - k / number_uneven_intervals))])

    return intervals



def parse_blast_line(line:str) -> dict:
    """
    Parse blastn output line into dictionary

    Args:
        line (str): blastn output line

    Returns:
        dict: dictionary containing hit length, pid and positions of query and subject
    """    

    ### load hit length, pid and posiiton of query and subject into variables
    hit_length = int(line.split("\t")[3])
    hit_pid = float(line.split("\t")[2])

    hit_q_start = int(line.split("\t")[6])
    hit_q_end = int(line.split("\t")[7])

    hit_s_start = int(line.split("\t")[8])
    hit_s_end = int(line.split("\t")[9])

    ### return as dictionary
    return {
        "hit_length": hit_length,
        "hit_pid": hit_pid,
        "hit_q_start": hit_q_start,
        "hit_q_end": hit_q_end,
        "hit_s_start": hit_s_start,
        "hit_s_end": hit_s_end
    }

def parse_blast_output(blast_output:list, interval:list, breakpoint_position:int) -> pd.DataFrame:
    """
    Parse blast output lines into dataframe

    Args:
        blast_output (str): Lines of blastn output in format 6

    Returns:
        pd.DataFrame: BLA
    """

    parsed_lines = []

    ### iterate over blast hits and parse
    for line in blast_output:

        ### skip empty lines
        if len(line.strip()) == 0:
            continue
            
        parsed_lines.append(parse_blast_line(line))
    
    #make dataframe
    blast_df = pd.DataFrame(parsed_lines)

    ### return empty DataFrame if no hits
    if blast_df.shape[0] == 0:
        return blast_df

    ### add information about interval and breakpoint position
    blast_df["interval_left"] = interval[0]
    blast_df["interval_right"] = interval[1]
    blast_df["bp_position"] = breakpoint_position

    ### add direction information
    blast_df["hit_q_dir"] = np.where(blast_df["hit_q_end"] >= blast_df["hit_q_start"], "+", "-")
    blast_df["hit_s_dir"] = np.where(blast_df["hit_s_end"] >= blast_df["hit_s_start"], "+", "-")

    ### add distance to breakpoint information
    blast_df["hit_q_bp_distance"] = np.where(
        blast_df["bp_position"].between(
            blast_df[["hit_q_start", "hit_q_end"]].min(axis=1), # type: ignore
            blast_df[["hit_q_start", "hit_q_end"]].max(axis=1) # type: ignore
        ),
        0,
        np.amin([
            (blast_df["bp_position"] - blast_df["hit_q_start"]).abs(),
            (blast_df["bp_position"] - blast_df["hit_q_end"]).abs()
        ], axis=0)
    )
    blast_df["hit_s_bp_distance"] = np.where(
        blast_df["bp_position"].between(
            blast_df[["hit_s_start", "hit_s_end"]].min(axis=1), # type: ignore
            blast_df[["hit_s_start", "hit_s_end"]].max(axis=1) # type: ignore
        ),
        0,
        np.amin([
            (blast_df["bp_position"] - blast_df["hit_s_start"]).abs(),
            (blast_df["bp_position"] - blast_df["hit_s_end"]).abs()
        ], axis=0)
    )

    ### add distance between breakpoint and left/right end of homology segment
    blast_df["hist_q_start_bp_dist"] = blast_df["bp_position"] - blast_df[["hit_q_start", "hit_q_end"]].min(axis=1) # type: ignore
    blast_df["hist_q_end_bp_dist"] = blast_df[["hit_q_start", "hit_q_end"]].max(axis=1) - blast_df["bp_position"] # type: ignore
    blast_df["hist_s_start_bp_dist"] = blast_df["bp_position"] - blast_df[["hit_s_start", "hit_s_end"]].min(axis=1) # type: ignore
    blast_df["hist_s_end_bp_dist"] = blast_df[["hit_s_start", "hit_s_end"]].max(axis=1) - blast_df["bp_position"] # type: ignore

    ### add hit delta - compare position of breakpoint in hit between s and q side
    col_hit_delta_left = []
    col_hit_delta_right = []
    for idx, row in blast_df.iterrows():
        col_hit_delta_left.append(np.abs(row["hist_q_start_bp_dist"] - row["hist_s_start_bp_dist"]))
        col_hit_delta_right.append(np.abs(row["hist_q_end_bp_dist"] - row["hist_s_end_bp_dist"]))
        """
        if row["hit_q_dir"] == row["hit_s_dir"]:
            col_hit_delta_left.append(np.abs(row["hist_q_start_bp_dist"] - row["hist_s_start_bp_dist"]))
            col_hit_delta_right.append(np.abs(row["hist_q_end_bp_dist"] - row["hist_s_end_bp_dist"]))
        else:
            col_hit_delta_left.append(np.abs(row["hist_q_start_bp_dist"] - row["hist_s_end_bp_dist"]))
            col_hit_delta_right.append(np.abs(row["hist_q_end_bp_dist"] - row["hist_s_start_bp_dist"]))
        """
    blast_df["delta_left"] = col_hit_delta_left
    blast_df["delta_right"] = col_hit_delta_right
    
    return blast_df



def blast(subject_seq:str, query_seq:str, 
          breakpoint_position:int,
          interval:list,
          workdir:str, 
          blastn_bin_path:str,
          direction:str,
          distance_to_breakpoint_cutoff:int,
          delta_breakpoints_cutoff:int) -> tuple:
    """
    This method blasts two sequences against each other to find candidate homologus regions. It calculates metrics and filters the candidates to exclude all
    where the homologous regions are too far away from the SV breakpoint, where the SV breakpoint is located very differently in the left and right homologus region,
    and the method checks for sensible directions of the homologous flanks (e.g. same direction for DEL and inverse direction for INV)

    Args:
        subject_seq (str): Sequence of the subject (e.g. window around left breakpoint)
        query_seq (str): Sequence of the query (e.g. window around right breakpoint)
        breakpoint_position (int): Position of the breakpoint within the sequences
        interval (list): Search window used for this iteration
        workdir (str): Working directory for temporary files
        blastn_bin_path (str): Path to the blastn binary
        direction (str): Filter for directionality of homologous flanks, possible values are "same", "reverse" and "any"
        distance_to_breakpoint_cutoff (int): Maximum allowed distance from breakpoint to homologous region
        delta_breakpoints_cutoff (int): Maximum allowed difference in breakpoint position between left and right homologous regions

    Returns:
        tuple: Tuple containing two DataFrames - one with filtered blast hits that passed all filters and one with hits that almost passed the filters
    """
    
    best_hit = {"length": 0, "pid": 0, "q_pos": [0, 0], "s_pos": [0, 0], "q_dir":".", "s_dir":"."}

    ### write subject and query sequence to file
    with open(os.path.join(workdir, "S.fa"), "w") as s_fo:
        s_fo.write(f">S\n{subject_seq}")
    with open(os.path.join(workdir, "Q.fa"), "w") as q_fo:
        q_fo.write(f">Q\n{query_seq}")

    ### run blastn, store output
    result = subprocess.run(
        [
            blastn_bin_path,
            "-subject",
            os.path.join(workdir, "S.fa"),
            "-query",
            os.path.join(workdir, "Q.fa"),
            "-word_size",
            "5",
            "-outfmt",
            "6",
            "-perc_identity",
            "80",
            "-dust",
            "no",
            "-soft_masking",
            "false"
        ], stdout=subprocess.PIPE)

    ### clean up written sequences
    os.system(f"rm {os.path.join(workdir, 'S.fa')}")
    os.system(f"rm {os.path.join(workdir, 'Q.fa')}")

    ### parse blast output
    blast_output = parse_blast_output (result.stdout.decode('ASCII').rstrip().split("\n"), interval=interval, breakpoint_position=breakpoint_position)
    
    ### return empty DataFrame if no hits
    if blast_output.shape[0] == 0:
        return (pd.DataFrame(), pd.DataFrame())

    ### filter hits
    # filter for hit to span the breakpoint
    blast_output["passes_spans_breakpoint_filter"] = (blast_output["hit_q_bp_distance"] <= distance_to_breakpoint_cutoff) & (blast_output["hit_s_bp_distance"] <= distance_to_breakpoint_cutoff)
    # filter for same direction - possible values are same, reverse, any
    if direction == "same":
        blast_output["passes_same_direction_filter"] = blast_output["hit_q_dir"] == blast_output["hit_s_dir"]
    elif direction == "reverse":
        blast_output["passes_same_direction_filter"] = blast_output["hit_q_dir"] != blast_output["hit_s_dir"]
    else:
        blast_output["passes_same_direction_filter"] = True
    # filter for relative position of breakpoint in homology segment
    blast_output["passes_relative_position_filter"] = (blast_output["delta_left"] <= delta_breakpoints_cutoff) & (blast_output["delta_right"] <= delta_breakpoints_cutoff)

    ### sort by hit length
    blast_output = blast_output.sort_values(by="hit_length", ascending=False)

    ### check how many filters are passed
    blast_output["NUM_PASSED"] = blast_output[["passes_spans_breakpoint_filter", "passes_same_direction_filter", "passes_relative_position_filter"]].sum(axis=1) # type: ignore

    ### return DataFrame of passing hits and DataFrame on almost passing hits
    return blast_output[blast_output["NUM_PASSED"] == 3].copy().reset_index(drop=True), blast_output[blast_output["NUM_PASSED"] >= 2].copy().reset_index(drop=True)  # type: ignore

def write_hit_as_bed (hit:pd.Series, record:cyvcf2.Variant, bed_outfile):
    """
    Write blast hit as bed intervals

    Args:
        hit (dict): Pandas series containing blast hit information
        record (cyvcf2.Variant): SV record
        bed_outfile (_type_): Opened BED file for outputting homologous regions
    """

    ### calculate left hit
    left_start = record.POS + hit["interval_left"] + hit['hit_s_start']
    left_end = record.POS + hit["interval_left"] + hit['hit_s_end']
    left_direction = "+"
    if left_end < left_start:
        left_direction = "-"
    
    ### write left homology
    bed_outfile.write("\t".join(
        [
            record.CHROM,
            str(np.amin([left_start, left_end])), 
            str(np.amax([left_start, left_end])),
            f"{record.ID}_LEFT",
            "1",
            left_direction
        ]
    ) + "\n")

    ### calculate right hit
    right_start = record.POS + np.abs(record.INFO['SVLEN']) + hit["interval_left"] + hit['hit_q_start']
    right_end = record.POS + np.abs(record.INFO['SVLEN']) + hit["interval_left"] + hit['hit_q_end']
    right_direction = "+"
    if right_end < right_start:
        right_direction = "-"
    
    ### write right homology
    bed_outfile.write("\t".join(
        [
            record.CHROM,
            str(np.amin([right_start, right_end])), 
            str(np.amax([right_start, right_end])),
            f"{record.ID}_RIGHT",
            "1",
            right_direction
        ]
    ) + "\n")


def detect_homology(record:cyvcf2.Variant, reference_file_path:str, workdir:str, blast_bin_path:str, max_search_window:int, bed_outfile, failing_hits=None) -> cyvcf2.Variant:
    """
    Detects homology for a single structural variant

    Args:
        record (cyvcf2.Variant): SV record to process
        reference_file_path (str): Reference genome FASTA file path
        workdir (str): Working directory, results will be written there
        blast_bin_path (str): Path to blastn binary
        max_search_window (int): Maximum search window size
        bed_outfile: Opened BED file for outputting homologous regions

    Raises:
        ValueError: This function only supports INS, DEL, INV and DUP SV types for now, for other SV types it raises an error

    Returns:
        cyvcf2.Variant: SV record with annotated homology information
    """    
    
    
    ### build search windows
    # set maximum length either the smaller of (a) supplied max_search_window or (b) SV length
    intervals = make_intervals(np.amin([np.abs(record.INFO["SVLEN"]), max_search_window]), 10)

    ### load reference genome
    reference = pysam.FastaFile(reference_file_path)

    ### build reference for INS
    if record.INFO["SVTYPE"] == "INS":
        ins_ref, ins_start, ins_end = make_reference_ins(record, reference)
        """
        # write insertion reference sequence and position of the insertion to file
        ins_fasta_fo.write(f">{record.ID.replace('>', '_')}\n{i_r}\n")
        ins_sv_bed_fo.write(f"{record.ID.replace('>', '_')}\t{i_s}\t{i_e}\n")
        """

    ### organize directionality
    directionality = "NA"
    if record.INFO["SVTYPE"] in ["INS", "DEL", "DUP"]:
        directionality = "same"
    elif record.INFO["SVTYPE"] == "INV":
        directionality = "reverse"

    ### prepare to store all hits
    hit_dfs = []
    candidate_dfs = []

    ### iterate over each interval and detect homology
    for interval in intervals:

        ### select subject and query sequence at SV breakpoints
        # select from reference genome for DEL, DUP and INV
        # select from insertion reference for INS
        q_seq = "N"
        s_seq = "N"
        if record.INFO["SVTYPE"] == "INS":
             s_seq = ins_ref[ins_start + interval[0]: ins_start + interval[1]]
             q_seq = ins_ref[ins_end + interval[0]: ins_end + interval[1]]
        elif record.INFO["SVTYPE"] in ["DEL", "INV", "DUP"]:
            s_seq = reference.fetch(record.CHROM, record.POS + interval[0], record.POS + interval[1])
            q_seq = reference.fetch(record.CHROM, record.POS + np.abs(record.INFO["SVLEN"]) + interval[0], record.POS + np.abs(record.INFO["SVLEN"]) + interval[1])
        else:
            raise ValueError("SVTYPE not supported")

        ### run blast to get best hit for this interval
        passing_hits, candidate_hits = blast(
            subject_seq=s_seq, 
            query_seq=q_seq, 
            breakpoint_position=np.abs(interval[0]),
            interval=interval,
            workdir=workdir, 
            blastn_bin_path=blast_bin_path, 
            direction=directionality,
            distance_to_breakpoint_cutoff=10,
            delta_breakpoints_cutoff=np.amin([20, np.abs(record.INFO["SVLEN"]) / 10])
        )

        ### store passing hits
        if passing_hits.shape[0] > 0:
            hit_dfs.append(passing_hits)

        ### store candidate hits
        if candidate_hits.shape[0] > 0:
            candidate_dfs.append(candidate_hits)

    ### concatenate all hits and sort by hit length
    try:
        all_hits_df = pd.concat(hit_dfs).sort_values(by="hit_length", ascending=False).reset_index(drop=True)
    except ValueError:
        all_hits_df = pd.DataFrame()
    
    ### concatenate candidate hits and sort by hit length
    try:
        candidate_hits_df = pd.concat(candidate_dfs).sort_values(by="hit_length", ascending=False).reset_index(drop=True)
    except ValueError:
        candidate_hits_df = pd.DataFrame()

    
    candidate_hits_df.to_csv(os.path.join(workdir, f"{record.ID}_candidates.tsv"), sep="\t", index=False)

    ### select the best hit and write
    if all_hits_df.shape[0] > 0:
        best_hit = all_hits_df.iloc[0]

        ### annotate VCF record
        record.INFO["HOMLEN"] = int(best_hit["hit_length"])
        record.INFO["HOMPID"] = float(best_hit["hit_pid"])

        write_hit_as_bed (hit=best_hit, record=record, bed_outfile=bed_outfile)
    
    ### if an output file for the candidate hits is provide write them there
    if failing_hits is not None and candidate_hits_df.shape[0] > 0:
        for idx, hit in candidate_hits_df[candidate_hits_df["hit_length"] >= 100].iterrows():
            write_hit_as_bed (hit=hit, record=record, bed_outfile=failing_hits)

    return record


if __name__ == "__main__":

    ### command line parameter
    parser = argparse.ArgumentParser(description="Detects homologous DNA flanks of SVs")
    parser.add_argument("--vcf_file", type=str, help="Path to the VCF file to be processed", required=True)
    parser.add_argument("--workdir", type=str, help="Path to the working directory", required=True)
    parser.add_argument("--reference_fasta", type=str, help="Path to the reference FASTA file", required=True)
    parser.add_argument("--blastn_bin", type=str, help="Path to the binary of blastn", required=True)
    parser.add_argument("--max_homology_length", type=int, help="Maximium length of size of search window", required=True)
    parser.add_argument("--job_name", type=str, help="Job name for the run", required=False)
    args = parser.parse_args()

    ### read parameter
    vcf_file_path = args.vcf_file
    workdir = args.workdir
    reference_full_path = args.reference_fasta
    blastn_bin = args.blastn_bin
    max_search_window = args.max_homology_length

    ### set job name, use either command line parameter or infer
    if args.job_name:
        job_name = args.job_name
    else:
        job_name = os.path.basename(vcf_file_path).replace(".vcf", "").replace(".vcf.gz", "")
    
    ### open vcf
    vcf_file = cyvcf2.VCF(vcf_file_path)
    vcf_file.add_info_to_header({"ID": "HOMLEN", "Description": "Homology length", "Type": "Integer", "Number": "1"})
    vcf_file.add_info_to_header({"ID": "HOMPID", "Description": "Percent identity of homologous regions", "Type": "Float", "Number": "1"})

    ### prepare annotated vcf
    annot_vcf_writer = cyvcf2.Writer(os.path.join(workdir, f"{job_name}.homology.vcf"), vcf_file)
    
    ### open bed file for string homologus regions
    hom_bed_fo = open(os.path.join(workdir, f"{job_name}.homology.bed"), "w")

    ### open bed file for candidate hits
    candidate_bed_fo = open(os.path.join(workdir, f"{job_name}.candidate_homology.bed"), "w")

    ### iterate over variants in VCF and detect homology
    for record in tqdm(vcf_file):

        ### skip INS for now
        if record.INFO["SVTYPE"] == "INS":
            continue

        ### detect homology
        ann_record = detect_homology(
            record,
            reference_file_path=reference_full_path,
            workdir=workdir,
            blast_bin_path=blastn_bin,
            max_search_window=max_search_window,
            bed_outfile = hom_bed_fo,
            failing_hits = candidate_bed_fo
        )

        ### write annotated record to VCF
        annot_vcf_writer.write_record(ann_record)

    ### close all
    hom_bed_fo.close()
    candidate_bed_fo.close()
    vcf_file.close()
    annot_vcf_writer.close()



#python homology_detection.py --vcf_file ../test/test.vcf --workdir ../../workdir/ --reference_fasta ../../../References/chm13v2.0.fa --blastn_bin ../../../../miniconda3/envs/jupyter_wave/bin/blastn --max_homology_length 10000 --job_name tmp