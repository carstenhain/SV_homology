import argparse
import os
import numpy as np
import cyvcf2
import pandas as pd

def annotate_vcf (vcf_file_path, out_vcf_file_path, repeatmasker_cmd, trf_cmd, rm2bed_cmd, sv_length_cutoff, workdir):
    ### extract all sv sequences - either inserted sequence or deleted part of reference

    if not workdir.endswith("/"):
        workdir += "/"

    out_fo = open(f"{workdir}sv_seq.fasta", "w")
    for record in cyvcf2.VCF(vcf_file_path):

        svlen = 10000000000
        if record.INFO["SVTYPE"] == "DEL":
            svlen = np.abs(record.INFO["END"] - record.POS)
        if record.INFO["SVTYPE"] == "INS":
            svlen = np.abs(record.INFO["INSLEN"])

        if svlen <= sv_length_cutoff:
            if record.INFO["SVTYPE"] == "DEL":
                out_fo.write(f">{record.ID}\n{record.REF}\n")
            if record.INFO["SVTYPE"] == "INS":
                out_fo.write(f">{record.ID}\n{record.ALT[0]}\n")
    out_fo.close()

    os.system(f"{repeatmasker_cmd} -pa 30 {workdir}sv_seq.fasta -dir {workdir}")
    os.system(f"python {rm2bed_cmd} {workdir}sv_seq.fasta.out -d {workdir}")

    bed_df = pd.read_csv(f"{workdir}sv_seq.fasta_rm.bed", sep="\t",
                         names=["CONTIG", "START", "END", "RM_DETAIL", "RM_LENGTH", "RM_DIRECTION", "RM_CLASS", "RM_GROUP", "NUMBER", "INDEX"])

    vcf = cyvcf2.VCF(vcf_file_path)
    vcf.add_info_to_header({"ID": "RM_ANNOT", "Description": "Type of RepeatMasker annotation", "Type": "String", "Number": "1"})
    writer = cyvcf2.Writer(out_vcf_file_path, vcf)

    k = 0
    for idx, record in enumerate(vcf):

        tmp_df = bed_df[bed_df["CONTIG"] == record.ID]

        rm_annot = "NA"

        if tmp_df.shape[0] > 0:

            svlen = len(record.ALT[0])
            sum_sine = np.sum(tmp_df[tmp_df["RM_CLASS"] == "SINE"]["RM_LENGTH"])
            sum_retro = np.sum(tmp_df[tmp_df["RM_CLASS"] == "Retroposon"]["RM_LENGTH"])
            sum_line = np.sum(tmp_df[tmp_df["RM_CLASS"] == "LINE"]["RM_LENGTH"])
            sum_satellite = np.sum(tmp_df[tmp_df["RM_CLASS"] == "Satellite"]["RM_LENGTH"])
            sum_simple_repeat = np.sum(tmp_df[tmp_df["RM_CLASS"] == "Simple_repeat"]["RM_LENGTH"])

            frac_sine = sum_sine / svlen
            frac_longest_sine = tmp_df[tmp_df["RM_CLASS"] == "SINE"]["RM_LENGTH"].max() / svlen
            num_sine = tmp_df[tmp_df["RM_CLASS"] == "SINE"].shape[0]
            frac_retro = sum_retro / svlen
            frac_longest_retro = tmp_df[tmp_df["RM_CLASS"] == "Retroposon"]["RM_LENGTH"].max() / svlen
            num_retro = tmp_df[tmp_df["RM_CLASS"] == "Retroposon"].shape[0]
            frac_line = sum_line / svlen
            frac_longest_line = tmp_df[tmp_df["RM_CLASS"] == "LINE"]["RM_LENGTH"].max() / svlen
            num_line = tmp_df[tmp_df["RM_CLASS"] == "LINE"].shape[0]
            frac_satellite = sum_satellite / svlen
            frac_longest_satellite = tmp_df[tmp_df["RM_CLASS"] == "Satellite"]["RM_LENGTH"].max() / svlen
            frac_simple_repeat = sum_simple_repeat / svlen
            frac_longest_simple_repeat = tmp_df[tmp_df["RM_CLASS"] == "Simple_repeat"]["RM_LENGTH"].max() / svlen

            if frac_sine > 0.8:
                if frac_longest_sine > 0.8:
                    rm_annot = "SINGLE_SINE"
                else:
                    if num_sine == 2:
                        rm_annot = "DOUBLE_SINE"
                    else:
                        rm_annot = "MULTI_SINE"
            if frac_retro > 0.8:
                if frac_longest_retro > 0.8:
                    rm_annot = "SINGLE_RETRO"
                else:
                    if num_retro == 2:
                        rm_annot = "DOUBLE_RETRO"
                    else:
                        rm_annot = "MULTI_RETRO"
            if frac_line > 0.8:
                if frac_longest_line > 0.8:
                    rm_annot = "SINGLE_LINE"
                else:
                    if num_line == 2:
                        rm_annot = "DOUBLE_LINE"
                    else:
                        rm_annot = "MULTI_LINE"
            if frac_satellite > 0.8:
                rm_annot = "SATELLITE"
            if frac_simple_repeat > 0.8:
                rm_annot = "SIMPLE_REPEAT"

        record.INFO["RM_ANNOT"] = rm_annot
        writer.write_record(record)

    writer.close()

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Annotates SV with RepeatMasker and Tandem Repeat Finder annotations")

    parser.add_argument("--vcf_file", type=str, help="Path to the VCF file containg SV calls")
    parser.add_argument("--annotated_vcf_file", type=str, help="Path to the annotated VCF file (output)")
    parser.add_argument("--repeatmasker_cmd", type=str, help="Path to binary of RepeatMasker")
    parser.add_argument("--trf_cmd", type=str, help="Path to the binary of tandem repeat finder")
    parser.add_argument("--rm2bed_cmd", type=str, help="Path to RM2Bed.py")
    parser.add_argument("--SV_length_cutoff", type=int, help="Maximum SV length to attempt annotation")
    parser.add_argument("--workdir", type=str, help="Path to working directory")

    args = parser.parse_args()

    annotate_vcf(args.vcf_file, args.annotated_vcf_file, args.repeatmasker_cmd, args.trf_cmd, args.rm2bed_cmd, args.SV_length_cutoff, args.workdir)

    #/g/easybuild/x86_64/Rocky/8/haswell/software/TRF/4.09.1-GCCcore-11.2.0/bin/trf
    #/g/easybuild/x86_64/Rocky/8/haswell/software/RepeatMasker/4.1.2-p1-foss-2021b/RepeatMasker
