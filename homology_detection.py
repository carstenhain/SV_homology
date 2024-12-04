import multiprocessing
import os
import cyvcf2
import datetime
from multiprocessing import Pool
import pysam
import subprocess
import numpy as np
from Bio.Seq import Seq
import microhomology
import argparse

################# get hit metrics #################
### including direction, distance to breakpoint and distance from breakpoint to left and right end of the homology segment
def get_hit_info(breakpoint_pos, hit_start, hit_end):
    info = {"sorted_coord": [np.amin([hit_start, hit_end]), np.amax([hit_start, hit_end])], "direction": "NA",
            "distance_to_breakpoint": 1000, "distance_hit_breakpoint": [0, 0]}

    if hit_end < hit_start:
        info["direction"] = "-"
    else:
        info["direction"] = "+"

    if breakpoint_pos >= info["sorted_coord"][0] and breakpoint_pos <= info["sorted_coord"][1]:
        info["distance_to_breakpoint"] = 0
    else:
        info["distance_to_breakpoint"] = np.amin(
            [np.abs(info["sorted_coord"][0] - breakpoint_pos), np.abs(info["sorted_coord"][1] - breakpoint_pos)])

    info["distance_hit_breakpoint"] = [breakpoint_pos - info["sorted_coord"][0],
                                       info["sorted_coord"][1] - breakpoint_pos]

    return info

############### build INS reference ###############
### build reference sequence for INS and return sequence and start, end position of the insertion
### TODO: change this depending on normalised or non-normalized vcf!!!
def make_reference_ins(record, reference):
    ins_ref = reference.fetch(record.CHROM, record.POS - 2 * record.INFO["SVLEN"], record.POS - 1).rstrip() + record.ALT[0].rstrip() + reference.fetch(record.CHROM, record.POS, record.POS + 2 * record.INFO["SVLEN"]).rstrip()
    return ins_ref, 2 * record.INFO["SVLEN"], 3 * record.INFO["SVLEN"]

################# build intervals #################
### based on SV length
def make_intervals(sv_length, number_uneven_intervals):
    intervals = []
    # add centered intervals up to half the SV length
    for i in [10, 15, 20, 25, 30, 35, 40, 45, 50, 75, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 1000, 2000, 3000, 4000, 5000, 10000, 50000, 100000, int(sv_length / 2) - 1]:
        if i < sv_length / 2 and i > 0:
            intervals.append([-i, i])
    # add uncentered intervals
    for k in range(0, number_uneven_intervals + 1):
        intervals.append(
            [-int(sv_length * k / number_uneven_intervals), int(sv_length * (1 - k / number_uneven_intervals))])
    return intervals

############### run and parse blast ###############
### blast subject against query, filter hits and return longest hit
def blast(subject_seq, query_seq, breakpoint_position, workdir, blast_cmd,
          same_direction,
          distance_to_breakpoint_cutoff,
          delta_breakpoints_cutoff):
    best_hit = {"length": 0, "pid": 0, "q_pos": [0, 0], "s_pos": [0, 0], "q_dir":".", "s_dir":"."}

    # write subject and query to file
    s_fo = open(f"{workdir}S.fa", "w")
    s_fo.write(f">S\n{subject_seq}")
    s_fo.close()

    q_fo = open(f"{workdir}Q.fa", "w")
    q_fo.write(f">Q\n{query_seq}")
    q_fo.close()

    # blast call
    result = subprocess.run(
        [
            blast_cmd,
            "-subject",
            f"{workdir}S.fa",
            "-query",
            f"{workdir}Q.fa",
            "-word_size",
            "5",
            "-outfmt",
            "6",
            "-perc_identity",
            "80"#,
            #"-dust",
            #"no",
            #"-soft_masking",
            #"false"
        ], stdout=subprocess.PIPE)

    # iterate over blast results
    for idx, line in enumerate(result.stdout.decode('ASCII').rstrip().split("\n")):

        if len(line.rstrip()) > 0:

            # TODO parse - change to one parsing method combine or expand get_hit_info
            hit_length = int(line.split("\t")[3])
            hit_pid = float(line.split("\t")[2])
            hit_q_pos = [int(line.split("\t")[6]), int(line.split("\t")[7])]
            hit_s_pos = [int(line.split("\t")[8]), int(line.split("\t")[9])]

            hit_q_info = get_hit_info(breakpoint_position, hit_q_pos[0], hit_q_pos[1])
            hit_s_info = get_hit_info(breakpoint_position, hit_s_pos[0], hit_s_pos[1])

            # difference of distance to left and distance to right end between the segment in query and subject --> used for restricting the results to thos where the breakpoint is at the relativly same position inside the homology strech
            if same_direction:
                hit_delta = \
                    np.abs(np.subtract(hit_q_info["distance_hit_breakpoint"], hit_s_info["distance_hit_breakpoint"]))
            else:
                hit_delta = [
                    np.abs(hit_q_info["distance_hit_breakpoint"][0] - hit_s_info["distance_hit_breakpoint"][1]),
                    np.abs(hit_q_info["distance_hit_breakpoint"][1] - hit_s_info["distance_hit_breakpoint"][0])]

            # filter for spans_breakpoint, same_direction and same relative position
            passes_spans_breakpoint_filter = False
            passes_same_direction_filter = False
            passes_relative_position_filter = True

            if same_direction:
                if hit_q_info["direction"] == hit_s_info["direction"]:
                    passes_same_direction_filter = True
            if not same_direction:
                if hit_q_info["direction"] != hit_s_info["direction"]:
                    passes_same_direction_filter = True

            if hit_q_info["distance_to_breakpoint"] <= distance_to_breakpoint_cutoff and hit_s_info[
                "distance_to_breakpoint"] <= distance_to_breakpoint_cutoff:
                passes_spans_breakpoint_filter = True

            if np.amax(hit_delta) > delta_breakpoints_cutoff:
                passes_relative_position_filter = False

            # evaluation of filtering
            if passes_spans_breakpoint_filter and passes_same_direction_filter and passes_relative_position_filter:
                if hit_length > best_hit["length"]:
                    best_hit = {"length": hit_length, "pid": hit_pid, "q_pos": hit_q_info["sorted_coord"],
                    "s_pos": hit_s_info["sorted_coord"], "breakpoint_balance": hit_delta, "s_dir":hit_s_info["direction"], "q_dir":hit_q_info["direction"]}

    return best_hit

##### worker - applies homology detection to one vcf ######
def detect_homology_vcf(worker_name, vcf_file_path, reference_file_path, workdir, resultdir, blast_cmd, max_search_window):

    print(f"This is worker {worker_name}: Analysing {vcf_file_path} and storing intermediates in {workdir} and results in {resultdir}")

    tmp_fo = open(f"{resultdir}flanks.bed", "w")

    ################## load reference #################

    reference = pysam.FastaFile(reference_file_path)

    ############ in out vcf, add INFO fields ##########

    vcf = cyvcf2.VCF(vcf_file_path)
    vcf.add_info_to_header({"ID": "BLASTHOMLEN", "Description": "Homology length by BLAST", "Type": "Integer", "Number": "1"})
    vcf.add_info_to_header({"ID": "BLASTHOMPID", "Description": "Homology identity by BLAST", "Type": "Integer", "Number": "1"})
    vcf.add_info_to_header({"ID": "MICROHOMLEN", "Description": "Microhomology length by DELLY", "Type": "Integer", "Number": "1"})
    vcf.add_info_to_header({"ID": "HOMLEN", "Description": "Homology length", "Type": "Integer", "Number": "1"})
    annot_vcf = cyvcf2.Writer(f"{resultdir}{vcf_file_path.split('/')[-1]}", vcf)

    ################ open output files ################

    ins_fasta_fo = open(f"{workdir}ins_ref.fa", "w")
    ins_sv_bed_fo = open(f"{workdir}ins_ref.bed", "w")
    ins_homology_bed_fo = open(f"{workdir}ins_ref.homology.bed", "w")
    del_homology_bed_fo = open(f"{workdir}del.homology.bed", "w")

    ############### iterate over records ##############

    k = 0
    total_start = datetime.datetime.now()

    for record in vcf:

        k += 1
        if (k - 1) % 100 == 0 and k != 1:
            print(f"Worker {worker_name} analysed {k - 1} SV in {datetime.datetime.now() - total_start}")

        if record.INFO["SVTYPE"] in ["INS", "DEL", "INV"]:

            ############### gather microhomology ##############

            mh = 0
            flanks = ["", ""]

            if record.INFO["SVTYPE"] == "INS":
                mh = microhomology.get_ins_microhomology(record, reference, flanks, 100000, workdir)
            if record.INFO["SVTYPE"] == "DEL":
                mh = microhomology.get_del_microhomology(record, reference, flanks, 1000000000, workdir)
            if record.INFO["SVTYPE"] == "INV":
                mh = microhomology.get_inv_microhomology(record, reference, flanks, 100, workdir)

            ################## gather homology ################

            ### build blank best hist for this SV
            best_hit = {"length": 0, "pid": 0, "ref_q_pos": [0, 0], "ref_s_pos": [0, 0],
                        "breakpoint_balance": [0, 0], "q_dir": ".", "s_dir": "."}

            ### build search windows based on SV length
            intervals = make_intervals(np.amin([np.abs(record.INFO["SVLEN"]), max_search_window]), 10)

            ### build reference for INS
            if record.INFO["SVTYPE"] == "INS":
                i_r, i_s, i_e = make_reference_ins(record, reference)
                ins_fasta_fo.write(f">{record.ID.replace('>', '_')}\n{i_r}\n")
                ins_sv_bed_fo.write(f"{record.ID.replace('>', '_')}\t{i_s}\t{i_e}\n")

            #### perform blast search for each interval ######

            for interval in intervals:

                ### build subject and query sequence, differs for DEL and INS
                q_seq = ""
                s_seq = ""
                if record.INFO["SVTYPE"] == "INS":
                    try:
                        s_seq = i_r[i_s + interval[0]: i_s + interval[1]]
                        q_seq = i_r[i_e + interval[0]: i_e + interval[1]]
                    except TypeError:
                        print(record.ID)
                if record.INFO["SVTYPE"] == "DEL":
                    s_seq = reference.fetch(record.CHROM, record.POS + interval[0], record.POS + interval[1])
                    q_seq = reference.fetch(record.CHROM, record.POS + np.abs(record.INFO["SVLEN"]) + interval[0], record.POS + np.abs(record.INFO["SVLEN"]) + interval[1])
                if record.INFO["SVTYPE"] == "INV":
                    s_seq = reference.fetch(record.CHROM, record.POS + interval[0], record.POS + interval[1])
                    q_seq = reference.fetch(record.CHROM, record.POS + np.abs(record.INFO["SVLEN"]) + interval[0],
                                            record.POS + np.abs(record.INFO["SVLEN"]) + interval[1])

                ### get blast hits
                #   filtering parameters
                #     distance to cutoff 10
                #     delta is either 50 or 1/10 of the SV length, which restricts this further for small SV
                #     TODO add directionality filter parameter when adding inversions
                same_direction = True
                if record.INFO["SVTYPE"] == "INV":
                    same_direction = False

                tmp_hit = blast(s_seq, q_seq, -interval[0], workdir, blast_cmd,
                                same_direction,
                                distance_to_breakpoint_cutoff=10,
                                delta_breakpoints_cutoff=np.amin([20, np.abs(record.INFO["SVLEN"]) / 10])) # alter cutoff war 50, mal drüber nachdenken ob 10 nicht der deutlich bessere cutoff ist oder ob man den flexibel gestalten sollte... (also min 10, max 50 und dazwischen scalierend bis 2 % der sv length

                ### check if the hit is better than the best hit
                #   put hit into best_hit, which differs for DEL and INS
                if tmp_hit["length"] > best_hit["length"]:
                    if record.INFO["SVTYPE"] == "INS":
                        best_hit = {
                            "length": tmp_hit["length"],
                            "pid": tmp_hit["pid"],
                            "ref_q_pos": [i_e + interval[0] + tmp_hit['q_pos'][0],
                                          i_e + interval[0] + tmp_hit['q_pos'][1]],
                            "ref_s_pos": [i_s + interval[0] + tmp_hit['s_pos'][0],
                                          i_s + interval[0] + tmp_hit['s_pos'][1]],
                            "breakpoint_balance": tmp_hit['breakpoint_balance'],
                            "q_dir": tmp_hit["q_dir"],
                            "s_dir": tmp_hit["s_dir"]
                        }
                    if record.INFO["SVTYPE"] == "DEL":
                        best_hit = {
                            "length": tmp_hit["length"],
                            "pid": tmp_hit["pid"],
                            "ref_q_pos": [record.POS + record.INFO["SVLEN"] + interval[0] + tmp_hit['q_pos'][0],
                                          record.POS + record.INFO["SVLEN"] + interval[0] + tmp_hit['q_pos'][1]],
                            "ref_s_pos": [record.POS + interval[0] + tmp_hit['s_pos'][0],
                                          record.POS + interval[0] + tmp_hit['s_pos'][1]],
                            "breakpoint_balance": tmp_hit['breakpoint_balance'],
                            "q_dir": tmp_hit["q_dir"],
                            "s_dir": tmp_hit["s_dir"]
                        }
                    if record.INFO["SVTYPE"] == "INV":
                        best_hit = {
                            "length": tmp_hit["length"],
                            "pid": tmp_hit["pid"],
                            "ref_q_pos": [record.POS + record.INFO["SVLEN"] + interval[0] + tmp_hit['q_pos'][0],
                                          record.POS + record.INFO["SVLEN"] + interval[0] + tmp_hit['q_pos'][1]],
                            "ref_s_pos": [record.POS + interval[0] + tmp_hit['s_pos'][0],
                                          record.POS + interval[0] + tmp_hit['s_pos'][1]],
                            "breakpoint_balance": tmp_hit['breakpoint_balance'],
                            "q_dir": tmp_hit["q_dir"],
                            "s_dir": tmp_hit["s_dir"]
                        }

            ################# write output ###################

            ### write homology segments of reference (DEL) or INS-reference (INS)
            if record.INFO["SVTYPE"] == "INS":
                ins_homology_bed_fo.write(f"{record.ID.replace('>', '_')}\t{best_hit['ref_s_pos'][0]}\t{best_hit['ref_s_pos'][1]}\t{record.ID}\t{np.amax(best_hit['breakpoint_balance'])}\t{best_hit['s_dir']}\n")
                ins_homology_bed_fo.write(f"{record.ID.replace('>', '_')}\t{best_hit['ref_q_pos'][0]}\t{best_hit['ref_q_pos'][1]}\t{record.ID}\t{np.amax(best_hit['breakpoint_balance'])}\t{best_hit['q_dir']}\n")
            if record.INFO["SVTYPE"] in ["DEL", "INV"]:
                del_homology_bed_fo.write(f"{record.CHROM}\t{best_hit['ref_s_pos'][0]}\t{best_hit['ref_s_pos'][1]}\t{record.ID}\t{np.amax(best_hit['breakpoint_balance'])}\t{best_hit['s_dir']}\n")
                del_homology_bed_fo.write(f"{record.CHROM}\t{best_hit['ref_q_pos'][0]}\t{best_hit['ref_q_pos'][1]}\t{record.ID}\t{np.amax(best_hit['breakpoint_balance'])}\t{best_hit['q_dir']}\n")

            ### write in info fields
            record.INFO["BLASTHOMLEN"] = best_hit["length"]
            record.INFO["BLASTHOMPID"] = best_hit["pid"]
            record.INFO["MICROHOMLEN"] = mh
            if mh <= 50 and mh > best_hit["length"]:
                record.INFO["HOMLEN"] = mh
            else:
                record.INFO["HOMLEN"] = best_hit["length"]
            #print(best_hit)
            ### write record
            annot_vcf.write_record(record)

    ################ close output files ###############

    ins_fasta_fo.close()
    ins_sv_bed_fo.close()
    ins_homology_bed_fo.close()
    del_homology_bed_fo.close()
    annot_vcf.close()
    tmp_fo.close()





####### main - parameter parsing and wrapper ######
if __name__ == "__main__":

    ### command line parameter
    parser = argparse.ArgumentParser(description="Detects homologous DNA flanks of SVs")
    parser.add_argument("--vcf_file", type=str, help="Path to the VCF file to be processed")
    parser.add_argument("--top_workdir", type=str, help="Path to the working directory")
    parser.add_argument("--reference_fasta", type=str, help="Path to the reference FASTA file")
    parser.add_argument("--blastn_bin", type=str, help="Path to the binary of blastn")
    parser.add_argument("--num_parallel", type=int, help="Number of parallel processes")
    parser.add_argument("--max_homology_length", type=int, help="Maximium length of size of search window")
    args = parser.parse_args()

    ### run parameter
    vcf_file_full_path = args.vcf_file
    top_workdir = args.top_workdir
    reference_full_path = args.reference_fasta
    blast_cmd = args.blastn_bin
    num_parallelization = args.num_parallel
    max_search_window = args.max_homology_length

    ### arguments for parallel processes

    args_name = []
    args_vcf = []
    args_ref = []
    args_wd = []
    args_rd = []
    args_blast_cmd = []
    args_max_search_window = []

    #################### split vcf ####################
    print("Split VCF")

    ### count variants
    num_variants = 0
    for record in cyvcf2.VCF(vcf_file_full_path):
        num_variants += 1

    ### size of each split
    bin_size = int(num_variants / num_parallelization) + 1

    ### report
    print(f"Found {num_variants}, splitting into {num_parallelization} files with bin size {bin_size}")

    ### initialize vcf streams
    vcf = cyvcf2.VCF(vcf_file_full_path)
    vcf_writer = []
    for i in range(0, num_parallelization):
        vcf_writer.append(cyvcf2.Writer(f"{top_workdir}tmp_S{i + 1}.vcf.gz", vcf))
        args_vcf.append(f"{top_workdir}tmp_S{i + 1}.vcf.gz")

    ### write in the respective stream
    for i, record in enumerate(vcf):
        vcf_writer[int(i / bin_size)].write_record(record)

    ### close all streams
    for writer in vcf_writer:
        writer.close()

    ################# build workdirs ##################
    os.system(f"rm -r {top_workdir}results")
    os.system(f"mkdir {top_workdir}results")
    for i in range(0, num_parallelization):
        os.system(f"rm -r {top_workdir}wd_S{i + 1}")
        os.system(f"mkdir {top_workdir}wd_S{i + 1}")
        args_wd.append(f"{top_workdir}wd_S{i + 1}/")
        args_rd.append(f"{top_workdir}results/")
        args_name.append(f"Worker {i}")
        args_ref.append(reference_full_path)
        args_blast_cmd.append(blast_cmd)
        args_max_search_window.append(max_search_window)

    ############### launch subprocesses ###############
    print("Launch subprocesses")
    p = Pool(num_parallelization)
    return_codes = p.starmap(detect_homology_vcf, zip(args_name, args_vcf, args_ref, args_wd, args_rd, args_blast_cmd, args_max_search_window))

    ############### concatenate results ###############

    vcf_concat = []
    del_hom_bed_concat = []
    ins_bed_concat = []
    ins_hom_bed_concat = []
    ins_fa_concat = []

    for i in range(0, num_parallelization):
        vcf_concat.append(f"{top_workdir}results/tmp_S{i + 1}.vcf.gz")
        del_hom_bed_concat.append(f"{top_workdir}wd_S{i + 1}/del.homology.bed")
        ins_bed_concat.append(f"{top_workdir}wd_S{i + 1}/ins_ref.bed")
        ins_hom_bed_concat.append(f"{top_workdir}wd_S{i + 1}/ins_ref.homology.bed")
        ins_fa_concat.append(f"{top_workdir}wd_S{i + 1}/ins_ref.fa")

    os.system(f"bcftools concat {' '.join(vcf_concat)} -o {vcf_file_full_path.replace('.vcf.gz', '.homology.vcf.gz')}")
    os.system(f"cat {' '.join(del_hom_bed_concat)} > {vcf_file_full_path.replace('.vcf.gz', '.del.homology.bed')}")
    os.system(f"cat {' '.join(ins_bed_concat)} > {vcf_file_full_path.replace('.vcf.gz', '.ins_ref.bed')}")
    os.system(f"cat {' '.join(ins_hom_bed_concat)} > {vcf_file_full_path.replace('.vcf.gz', '.ins_ref.homology.bed')}")
    os.system(f"cat {' '.join(ins_fa_concat)} > {vcf_file_full_path.replace('.vcf.gz', '.ins_ref.fa')}")

    os.system(f"gzip -dfk {vcf_file_full_path.replace('.vcf.gz', '.homology.vcf.gz')}")

    ##################### clean up ####################
    """
    os.system(f"rm -r {top_workdir}results")
    for i in range(0, num_parallelization):
        os.system(f"rm -r {top_workdir}wd_S{i + 1}")
        os.system(f"rm {top_workdir}tmp_S{i + 1}.vcf.gz")
    """

