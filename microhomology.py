import os
import cyvcf2
import pysam
import subprocess
import numpy as np
from Bio.Seq import Seq
import random

########## get microhomology using delly ##########
### different functions for DEL, INS and INV
def get_del_microhomology(cyvcf2_record, reference, flanks, del_cutoff, workdir):
    start = cyvcf2_record.POS
    end = cyvcf2_record.POS + len(cyvcf2_record.REF)
    bases = ["A", "C", "G", "T"]
    padding = [[0, 100], [50, 50], [100, 0]]

    # build reference
    if end - start > del_cutoff:
        del_seq = reference.fetch(cyvcf2_record.CHROM, start, start + int(del_cutoff / 2)) + reference.fetch(cyvcf2_record.CHROM, end - int(del_cutoff / 2), end)
    else:
        del_seq = reference.fetch(cyvcf2_record.CHROM, start, end)
    ref_fo = open(f"{workdir}ref.fa", "w")
    ref_fo.write(">ref\n" + flanks[0] + reference.fetch(cyvcf2_record.CHROM, start - 2000, start) + del_seq + reference.fetch(cyvcf2_record.CHROM, end, end + 2000) + flanks[1] + "\n")
    ref_fo.close()

    # build reads
    fl_read = flanks[0] + reference.fetch(cyvcf2_record.CHROM, start - 2000, start) + reference.fetch(cyvcf2_record.CHROM, end, end + 2000) + flanks[1]
    reads = open(f"{workdir}alt.fa", "w")
    for idx, read in enumerate(["A", "B", "C"]):
        mut_positions = []
        mut_bases = []
        if len(fl_read) - 1 < 0:
            print("Read with sub 0 length")
            print(cyvcf2_record.ID)
        for k in range(0, 50):
            mut_positions.append(random.randrange(0, 1500 + len(flanks[0])))
            mut_bases.append(bases[random.randrange(0, 4)])
            try:
                mut_positions.append(random.randrange(len(flanks[0]) + 2500, len(fl_read) - 1))
            except ValueError:
                print("ValueError")
                print(cyvcf2_record.ID)
            for k in range(0, 50):
                mut_bases.append(bases[random.randrange(0, 4)])
        mod_read = list(fl_read)
        for j in range(0, len(mut_positions)):
            mod_read[mut_positions[j]] = mut_bases[j]
        reads.write(f">{read}\n" + "".join(mod_read)[padding[idx][0]:] + "\n")
    reads.close()

    # wipe index
    if os.path.exists(f"rm {workdir}*.fai"):
        os.system(f"rm {workdir}*.fai")
    if os.path.exists(f"rm {workdir}*.bai"):
        os.system(f"rm {workdir}*.bai")
    if os.path.exists(f"rm {workdir}*.csi"):
        os.system(f"rm {workdir}*.csi")

    # mapping and variant calling
    os.system(f"samtools faidx {workdir}ref.fa")
    subprocess.run(
        f"minimap2 -ax map-ont {workdir}ref.fa {workdir}alt.fa | samtools view -b - | samtools sort - > {workdir}alt_to_ref.bam",
        shell=True, stderr=subprocess.PIPE)
    os.system(f"samtools index {workdir}alt_to_ref.bam")
    delly_return = subprocess.run(f"delly lr -t DEL -g {workdir}ref.fa {workdir}alt_to_ref.bam -o {workdir}out.bcf", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if delly_return.returncode != 0:
        print(f"Error {cyvcf2_record.ID}")
        print(f"Error {delly_return.stderr.decode('ASCII')}")
    if not os.path.exists(f"{workdir}out.bcf"):
        print(f"Nonexisting file {cyvcf2_record.ID}")
    os.system(f"bcftools convert -O v {workdir}out.bcf -o {workdir}out.vcf")

    # extract homlen and write into variant
    homlen = -1
    for record in cyvcf2.VCF(f"{workdir}out.vcf"):
        homlen = record.INFO["HOMLEN"]
    if homlen > int(del_cutoff / 2):
        homlen = int(del_cutoff / 2)
    if homlen > (end - start):
        homlen = (end - start)

    return homlen

def get_ins_microhomology(cyvcf2_record, reference, flanks, ins_cutoff, workdir):
    svlen = len(cyvcf2_record.ALT[0])
    bases = ["A", "C", "G", "T"]
    padding = [[0, 100], [50, 50], [100, 0]]

    # build reference
    ref_fo = open(f"{workdir}ref.fa", "w")
    ref_fo.write(
        ">ref\n" + flanks[0] +
        reference.fetch(cyvcf2_record.CHROM, cyvcf2_record.POS - 2000, cyvcf2_record.POS + 2000) +
        flanks[1] + "\n")
    ref_fo.close()

    # build reads
    if svlen > ins_cutoff:
        ins_seq = cyvcf2_record.ALT[0][1:][:int(ins_cutoff / 2)] + cyvcf2_record.ALT[0][1:][-int(ins_cutoff / 2):]
    else:
        ins_seq = cyvcf2_record.ALT[0][1:]
    fl_read = flanks[0] + reference.fetch(cyvcf2_record.CHROM, cyvcf2_record.POS - 2000,
                                          cyvcf2_record.POS) + ins_seq + reference.fetch(cyvcf2_record.CHROM,
                                                                                         cyvcf2_record.POS,
                                                                                         cyvcf2_record.POS + 2000) + \
              flanks[1]
    reads = open(f"{workdir}alt.fa", "w")
    for idx, read in enumerate(["D", "E", "F"]):
        mut_positions = []
        mut_bases = []
        for k in range(0, 50):
            mut_positions.append(random.randrange(0, 1500 + len(flanks[0])))
            mut_bases.append(bases[random.randrange(0, 4)])
            mut_positions.append(random.randrange(len(flanks[0]) + 2000 + len(ins_seq) + 500, len(fl_read) - 1))
            mut_bases.append(bases[random.randrange(0, 4)])
        mod_read = list(fl_read)
        for j in range(0, len(mut_positions)):
            mod_read[mut_positions[j]] = mut_bases[j]
        reads.write(f">{read}\n" + "".join(mod_read)[padding[idx][0]:] + "\n")
    reads.close()

    # wipe index
    if os.path.exists(f"rm {workdir}*.fai"):
        os.system(f"rm {workdir}*.fai")
    if os.path.exists(f"rm {workdir}*.bai"):
        os.system(f"rm {workdir}*.bai")
    if os.path.exists(f"rm {workdir}*.csi"):
        os.system(f"rm {workdir}*.csi")

    # mapping and variant calling
    os.system(f"samtools faidx {workdir}ref.fa")
    subprocess.run(
        f"minimap2 -ax map-ont {workdir}ref.fa {workdir}alt.fa | samtools view -b - | samtools sort - > {workdir}alt_to_ref.bam",
        shell=True, stderr=subprocess.PIPE)
    os.system(f"samtools index {workdir}alt_to_ref.bam")
    delly_return = subprocess.run(f"delly lr -t INS -g {workdir}ref.fa {workdir}alt_to_ref.bam -o {workdir}out.bcf", shell=True,
                   stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if delly_return.returncode != 0:
        print(f"Error {cyvcf2_record.ID}")
    if not os.path.exists(f"{workdir}out.bcf"):
        print(f"Nonexisting file {cyvcf2_record.ID}")
    os.system(f"bcftools convert -O v {workdir}out.bcf -o {workdir}out.vcf")

    # extract homlen and write into variant
    homlen = -1
    for record in cyvcf2.VCF(f"{workdir}out.vcf"):
        homlen = record.INFO["HOMLEN"]
    if homlen > int(ins_cutoff / 2):
        homlen = int(ins_cutoff / 2)
    if homlen > svlen:
        homlen = svlen

    return homlen

def get_inv_microhomology(cyvcf2_record, reference, flanks, inv_cutoff, workdir):

    start = cyvcf2_record.POS
    end = cyvcf2_record.INFO["END"]
    bases = ["A", "C", "G", "T"]
    padding = [[0, 100], [50, 50], [100, 0]]

    # build reference
    if end - start > inv_cutoff:
        inv_seq = reference.fetch(cyvcf2_record.CHROM, start, start + int(inv_cutoff / 2)) + reference.fetch(cyvcf2_record.CHROM, end - int(del_cutoff / 2), end)
    else:
        inv_seq = reference.fetch(cyvcf2_record.CHROM, start, end)

    ref_fo = open(f"{workdir}ref.fa", "w")
    ref_fo.write(">ref\n" + flanks[0] + reference.fetch(cyvcf2_record.CHROM, start - 4000, start) + inv_seq + reference.fetch(cyvcf2_record.CHROM, end, end + 4000) + flanks[1] + "\n")
    ref_fo.close()

    # build reads
    fl_read = flanks[0] + reference.fetch(cyvcf2_record.CHROM, start - 4000, start) + Seq(inv_seq).reverse_complement() + reference.fetch(cyvcf2_record.CHROM, end, end + 4000) + flanks[1]
    reads = open(f"{workdir}alt.fa", "w")
    for idx, read in enumerate(["A", "B", "C"]):
        mut_positions = []
        mut_bases = []
        if len(fl_read) - 1 < 0:
            print("Read with sub 0 length")
            print(cyvcf2_record.ID)

        for k in range(0, 50):
            mut_positions.append(random.randrange(0, 3500 + len(flanks[0])))
            mut_bases.append(bases[random.randrange(0, 4)])
            try:
                mut_positions.append(random.randrange(len(flanks[0]) + 4500 + len(inv_seq), len(fl_read) - 1))
                mut_bases.append(bases[random.randrange(0, 4)])
            except ValueError:
                print("ValueError")
                print(cyvcf2_record.ID)

        mod_read = list(fl_read)
        for j in range(0, len(mut_positions)):
            mod_read[mut_positions[j]] = mut_bases[j]
        reads.write(f">{read}\n" + "".join(mod_read)[padding[idx][0]:] + "\n")
    reads.close()

    # wipe index
    if os.path.exists(f"rm {workdir}*.fai"):
        os.system(f"rm {workdir}*.fai")
    if os.path.exists(f"rm {workdir}*.bai"):
        os.system(f"rm {workdir}*.bai")
    if os.path.exists(f"rm {workdir}*.csi"):
        os.system(f"rm {workdir}*.csi")

    # mapping and variant calling
    os.system(f"samtools faidx {workdir}ref.fa")
    subprocess.run(
        f"ngmlr -r {workdir}ref.fa -q {workdir}alt.fa -o {workdir}tmp.sam",
        shell=True, stderr=subprocess.PIPE)
    subprocess.run(
        f"samtools sort {workdir}tmp.sam -o {workdir}sort.sam",
        shell=True, stderr=subprocess.PIPE)
    subprocess.run(
        f"samtools view -bh {workdir}sort.sam -o {workdir}alt_to_ref.bam",
        shell=True, stderr=subprocess.PIPE)
    os.system(f"samtools index {workdir}alt_to_ref.bam")
    delly_return = subprocess.run(
        f"delly lr -t INV -g {workdir}ref.fa {workdir}alt_to_ref.bam -o {workdir}out.bcf", shell=True,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if delly_return.returncode != 0:
        print(f"Error {cyvcf2_record.ID}")
    if not os.path.exists(f"{workdir}out.bcf"):
        print(f"Nonexisting file {cyvcf2_record.ID}")
    os.system(f"bcftools convert -O v {workdir}out.bcf -o {workdir}out.vcf")

    # extract homlen and write into variant
    homlen = -1
    for record in cyvcf2.VCF(f"{workdir}out.vcf"):
        homlen = record.INFO["HOMLEN"]
    if homlen > int(inv_cutoff / 2):
        homlen = int(inv_cutoff / 2)
    if homlen > (end - start):
        homlen = (end - start)

    return homlen