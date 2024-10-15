import argparse
import os
import cyvcf2
import numpy as np
import pandas as pd
from collections import Counter


def annotate (
        vcf_file,
        homology_bed_full_path,
        insertion_homology_bed_full_path,
        insertion_alleles_fasta_repeatmasker_bed_full_path,
        repeatmasker_bed_full_path,
        sd_bed_full_path,
        bedtools_cmd):

    ### annotate deletions with RepeatMasker

    # intersect reference homologies with the RepeatMasker track
    annotated_homology_bed_full_path = homology_bed_full_path.replace(".bed", ".RM_INTERSECT.bed")
    #os.system(f"{bedtools_cmd} intersect -wo -a {homology_bed_full_path} -b {repeatmasker_bed_full_path} > {annotated_homology_bed_full_path}")

    # read intersection into dataframe
    custom_header = ["HOM_CHROM", "HOM_START", "HOM_END", "SV_ID", "SCORE", "HOM_DIRECTION", "RM_CHROM", "RM_START",
                     "RM_END", "RM_DETAIL", "?", "RM_DIRECTION", "RM_GROUP", "RM_CLASS", "RM_SCORE_1", "RM_SCORE_2",
                     "OVERLAP"]
    bed_df = pd.read_csv(annotated_homology_bed_full_path, sep="\t", header=None, names=custom_header)
    frac_hom_overlap = []
    frac_rm_overlap = []
    for idx, row in bed_df.iterrows():
        if (row["HOM_END"] - row["HOM_START"]) > 0:
            frac_hom_overlap.append(row["OVERLAP"] / (row["HOM_END"] - row["HOM_START"]))
        else:
            frac_hom_overlap.append(0)
        frac_rm_overlap.append(row["OVERLAP"] / (row["RM_END"] - row["RM_START"]))
    bed_df.insert(17, "FRACTION HOMOLOGY OVERLAP", frac_hom_overlap)
    bed_df.insert(18, "FRACTION RM OVERLAP", frac_rm_overlap)

    # remove very partial intersection
    fil_bed_df = bed_df[(bed_df["FRACTION HOMOLOGY OVERLAP"] > 0.85) & (bed_df["FRACTION RM OVERLAP"] > 0.1)]

    # perform classification
    classifications_RM_deletions = {}

    for idx, sv in enumerate(set(fil_bed_df["SV_ID"])):
        tmp = fil_bed_df[fil_bed_df["SV_ID"] == sv]
        start_positions = list(set(tmp["HOM_START"]))
        if len(start_positions) == 2:
            rm_overlap_start = tmp[tmp["HOM_START"] == start_positions[0]]
            rm_overlap_end = tmp[tmp["HOM_START"] == start_positions[1]]
            if np.amax([rm_overlap_start["FRACTION RM OVERLAP"].max(), rm_overlap_end["FRACTION RM OVERLAP"].max()]) > 0.85:
                best_rm_overlap_start = rm_overlap_start[
                    rm_overlap_start["OVERLAP"] == rm_overlap_start["OVERLAP"].max()]
                best_rm_overlap_end = rm_overlap_end[rm_overlap_end["OVERLAP"] == rm_overlap_end["OVERLAP"].max()]
                distance = np.amin(np.abs(
                    [best_rm_overlap_end["RM_START"].iloc[0] - best_rm_overlap_start["RM_END"].iloc[0],
                     best_rm_overlap_start["RM_START"].iloc[0] - best_rm_overlap_end["RM_END"].iloc[0]]))
                if best_rm_overlap_start["RM_CLASS"].iloc[0] == best_rm_overlap_end["RM_CLASS"].iloc[
                    0] and distance > 20:
                    classifications_RM_deletions[best_rm_overlap_start["SV_ID"].iloc[0]] = {
                        "CLASS": best_rm_overlap_start["RM_CLASS"].iloc[0],
                        "GROUP": best_rm_overlap_start["RM_GROUP"].iloc[0],
                        "START_ELEMENT": best_rm_overlap_start["RM_DETAIL"].iloc[0],
                        "END_ELEMENT": best_rm_overlap_end["RM_DETAIL"].iloc[0]
                    }

    print(f"Finished DEL annotation with RepeatMasker elements, found {len(classifications_RM_deletions)} annotated variants")
    classes = []
    for key in classifications_RM_deletions:
        classes.append(classifications_RM_deletions[key]["CLASS"])
    class_counts_dict = dict(Counter(classes))
    for key in class_counts_dict:
        print(f"{key}\t{class_counts_dict[key]}")

    ### annotate insertions with RepeatMasker
    insertion_annotated_homology_bed_full_path = insertion_homology_bed_full_path.replace(".bed", ".RM_INTERSECT.bed")
    #os.system(f"{bedtools_cmd} intersect -wo -a {insertion_homology_bed_full_path} -b {insertion_alleles_fasta_repeatmasker_bed_full_path} > {insertion_annotated_homology_bed_full_path}")

    custom_header = ["HOM_CHROM", "HOM_START", "HOM_END", "SV_ID", "SCORE", "HOM_DIRECTION", "RM_CHROM", "RM_START",
                     "RM_END", "RM_DETAIL", "?", "RM_DIRECTION", "RM_GROUP", "RM_CLASS", "RM_SCORE_1", "RM_SCORE_2",
                     "OVERLAP"]
    bed_df = pd.read_csv(insertion_annotated_homology_bed_full_path, sep="\t", header=None, names=custom_header)
    frac_hom_overlap = []
    frac_rm_overlap = []
    for idx, row in bed_df.iterrows():
        if (row["HOM_END"] - row["HOM_START"]) > 0:
            frac_hom_overlap.append(row["OVERLAP"] / (row["HOM_END"] - row["HOM_START"]))
        else:
            frac_hom_overlap.append(0)
        frac_rm_overlap.append(row["OVERLAP"] / (row["RM_END"] - row["RM_START"]))
    bed_df.insert(17, "FRACTION HOMOLOGY OVERLAP", frac_hom_overlap)
    bed_df.insert(18, "FRACTION RM OVERLAP", frac_rm_overlap)

    fil_bed_df = bed_df[(bed_df["FRACTION HOMOLOGY OVERLAP"] > 0.85) & (bed_df["FRACTION RM OVERLAP"] > 0.1)]

    classifications_RM_insertions = {}

    for idx, sv in enumerate(set(fil_bed_df["SV_ID"])):
        tmp = fil_bed_df[fil_bed_df["SV_ID"] == sv]
        start_positions = list(set(tmp["HOM_START"]))
        if len(start_positions) == 2:
            rm_overlap_start = tmp[tmp["HOM_START"] == start_positions[0]]
            rm_overlap_end = tmp[tmp["HOM_START"] == start_positions[1]]
            if np.amax([rm_overlap_start["FRACTION RM OVERLAP"].max(),
                        rm_overlap_end["FRACTION RM OVERLAP"].max()]) > 0.85:
                best_rm_overlap_start = rm_overlap_start[
                    rm_overlap_start["OVERLAP"] == rm_overlap_start["OVERLAP"].max()]
                best_rm_overlap_end = rm_overlap_end[rm_overlap_end["OVERLAP"] == rm_overlap_end["OVERLAP"].max()]
                distance = np.amin(np.abs(
                    [best_rm_overlap_end["RM_START"].iloc[0] - best_rm_overlap_start["RM_END"].iloc[0],
                     best_rm_overlap_start["RM_START"].iloc[0] - best_rm_overlap_end["RM_END"].iloc[0]]))
                if best_rm_overlap_start["RM_CLASS"].iloc[0] == best_rm_overlap_end["RM_CLASS"].iloc[
                    0] and distance > 20:
                    classifications_RM_insertions[best_rm_overlap_start["SV_ID"].iloc[0]] = {
                        "CLASS": best_rm_overlap_start["RM_CLASS"].iloc[0],
                        "GROUP": best_rm_overlap_start["RM_GROUP"].iloc[0],
                        "START_ELEMENT": best_rm_overlap_start["RM_DETAIL"].iloc[0],
                        "END_ELEMENT": best_rm_overlap_end["RM_DETAIL"].iloc[0]
                    }
    print(f"Finished INS annotation with RepeatMasker elements, found {len(classifications_RM_insertions)} annotated variants")
    classes = []
    for key in classifications_RM_insertions:
        classes.append(classifications_RM_insertions[key]["CLASS"])
    class_counts_dict = dict(Counter(classes))
    for key in class_counts_dict:
        print(f"{key}\t{class_counts_dict[key]}")


    ### annotate deletions with SD

    # intersect reference homologies with the RepeatMasker track
    sd_homology_bed_full_path = homology_bed_full_path.replace(".bed", ".SD_INTERSECT.bed")
    #os.system(f"{bedtools_cmd} intersect -wo -a {homology_bed_full_path} -b {sd_bed_full_path} > {sd_homology_bed_full_path}")

    # read intersection into dataframe
    custom_header = ["HOM_CHROM", "HOM_START", "HOM_END", "SV_ID", "SCORE", "HOM_DIRECTION", "SD_CHROM", "SD_START",
                     "SD_END", "SD_PARTNER_CHROM", "SD_PARTNER_START", "SD_PARTNER_END", "OVERLAP"]
    bed_df = pd.read_csv(sd_homology_bed_full_path, sep="\t", header=None, names=custom_header)
    frac_hom_overlap = []
    frac_sd_overlap = []
    sd_id = []
    for idx, row in bed_df.iterrows():
        if (row["HOM_END"] - row["HOM_START"]) > 0:
            frac_hom_overlap.append(row["OVERLAP"] / (row["HOM_END"] - row["HOM_START"]))
        else:
            frac_hom_overlap.append(0)
        frac_sd_overlap.append(row["OVERLAP"] / (row["SD_END"] - row["SD_START"]))

        if row["SD_START"] < row["SD_PARTNER_START"]:
            sd_id.append(f"{row['SD_CHROM']}:{row['SD_START']}-{row['SD_END']}_{row['SD_PARTNER_CHROM']}:{row['SD_PARTNER_START']}-{row['SD_PARTNER_END']}")
        else:
            sd_id.append(f"{row['SD_PARTNER_CHROM']}:{row['SD_PARTNER_START']}-{row['SD_PARTNER_END']}_{row['SD_CHROM']}:{row['SD_START']}-{row['SD_END']}")

    bed_df.insert(13, "FRACTION_HOMOLOGY_OVERLAP", frac_hom_overlap)
    bed_df.insert(14, "FRACTION_SD_OVERLAP", frac_sd_overlap)
    bed_df.insert(15, "SD_ID", sd_id)

    # remove very partial intersection
    fil_bed_df = bed_df[(bed_df["FRACTION_HOMOLOGY_OVERLAP"] > 0.85) & (bed_df["OVERLAP"] >= 1000)]

    classifications_SD_deletions = {}

    for idx, sv in enumerate(set(fil_bed_df["SV_ID"])):
        sd_idents = []
        tmp = fil_bed_df[fil_bed_df["SV_ID"] == sv]
        start_positions = list(set(tmp["HOM_START"]))
        if len(start_positions) == 2:
            sd_overlap_start = tmp[tmp["HOM_START"] == start_positions[0]]
            sd_overlap_end = tmp[tmp["HOM_START"] == start_positions[1]]
            for sd_idx, sd in sd_overlap_start.iterrows():
                if sd_overlap_end[sd_overlap_end["SD_ID"] == sd['SD_ID']].shape[0] > 0:
                    sd_idents.append(sd['SD_ID'])

            classifications_SD_deletions[sv] = sd_idents

    print(f"Finished DEL annotation with SD elements, found {len(classifications_SD_deletions)} annotated variants")

    ### write classifications into VCF file

    ###########################################################################
    ########### COMBINE THIS FOR DEL RM, DEL SD, INS RM
    ###########################################################################


    vcf = cyvcf2.VCF(vcf_file)
    vcf.add_info_to_header({"ID": "RM_GROUP", "Description": "Group of repeat element mediating this deletion", "Type": "String", "Number": "1"})
    vcf.add_info_to_header({"ID": "RM_CLASS", "Description": "Class of repeat element mediating this deletion", "Type": "String", "Number": "1"})
    vcf.add_info_to_header({"ID": "RM_START", "Description": "Detail of repeat element at SV start", "Type": "String", "Number": "1"})
    vcf.add_info_to_header({"ID": "RM_END", "Description": "Detail of repeat element at SV end", "Type": "String", "Number": "1"})
    vcf.add_info_to_header({"ID": "SD", "Description": "Identifiers of SD likely mediating this SV, NA if not mediated by SDs", "Type": "String", "Number": "1"})
    writer = cyvcf2.Writer(vcf_file.replace(".vcf.gz", ".hom_annot.vcf.gz"), vcf)

    for record in vcf:
        if record.ID in classifications_RM_deletions:
            record.INFO["RM_GROUP"] = classifications_RM_deletions[record.ID]["GROUP"]
            record.INFO["RM_CLASS"] = classifications_RM_deletions[record.ID]["CLASS"]
            record.INFO["RM_START"] = classifications_RM_deletions[record.ID]["START_ELEMENT"]
            record.INFO["RM_END"] = classifications_RM_deletions[record.ID]["END_ELEMENT"]
        else:
            record.INFO["RM_GROUP"] = "NA"
            record.INFO["RM_CLASS"] = "NA"
            record.INFO["RM_START"] = "NA"
            record.INFO["RM_END"] = "NA"

        if record.ID in classifications_SD_deletions:
            record.INFO["SD"] = ",".join(classifications_SD_deletions[record.ID])
        else:
            record.INFO["SD"] = "NA"

        writer.write_record(record)

    writer.close()


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Annotates homologous DNA flanks of SVs using RepeatMasker and Segmental Duplications")

    parser.add_argument("--vcf_file", type=str, help="Path to the VCF file processed with homology_detection.py")
    parser.add_argument("--reference_homology_bed", type=str, help="Path to the BED file with the SV homology flanks on the reference (DEL, INV)")
    parser.add_argument("--insertion_homology_bed", type=str, help="Path to the BED file with the SV homology flanks on the insertion alleles")
    parser.add_argument("--insertion_fasta_repeatmasker_bed", type=str, help="Path to the RepeatMasker annotations (in bed format) for the FASTA file containg the insertion alleles")
    parser.add_argument("--repeatmasker_bed", type=str, help="Path to the BED file containing RepeatMasker annotations for the reference genome (can be gzipped)")
    parser.add_argument("--sd_bed", type=str, help="Path to the BED file containing Segmental Duplications for the reference genome (can be gzipped)")
    parser.add_argument("--bedtools_bin", type=str, help="Path to the binary of BEDtools")

    args = parser.parse_args()

    annotate(args.vcf_file, args.reference_homology_bed, args.insertion_homology_bed, args.insertion_fasta_repeatmasker_bed, args.repeatmasker_bed, args.sd_bed, args.bedtools_bin)