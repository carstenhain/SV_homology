import cyvcf2 # type: ignore
import argparse
from tqdm import tqdm # type: ignore
import os

if __name__ == "__main__":

    ### command line parameter
    parser = argparse.ArgumentParser(description="Quick script to a large vcf file into smaller chunks for parallel processing")
    parser.add_argument("--vcf_file", type=str, help="Path to the VCF file to be processed", required=True)
    parser.add_argument("--num_splits", type=int, help="Number of splits to create", required=True)
    parser.add_argument("--workdir", type=str, help="Path to the working directory", required=True)
    args = parser.parse_args()

    ### count variants
    num_variants = 0
    for record in tqdm(cyvcf2.VCF(args.vcf_file), desc="Counting variants in VCF"):
        num_variants += 1

    ### size of each split
    bin_size = int(num_variants / args.num_splits) + 1

    ### report
    print(f"Found {num_variants}, splitting into {args.num_splits} files with {bin_size} variants each")

    ### determine vcf format
    vcf_format = "NA"
    if args.vcf_file.endswith(".vcf"):
        vcf_format = "vcf"
    elif args.vcf_file.endswith(".vcf.gz"):
        vcf_format = "vcf.gz"
    else:
        raise ValueError("Only vcf and vcf.gz files are supported")

    ### initialize vcf streams
    vcf = cyvcf2.VCF(args.vcf_file)
    vcf_writer = []
    for i in range(0, args.num_splits):
        vcf_writer.append(cyvcf2.Writer(os.path.join(args.workdir, f"split_{i + 1}.{vcf_format}"), vcf))

    ### write in the respective stream
    for i, record in enumerate(tqdm(vcf, desc="Splitting VCF", total=num_variants)):
        vcf_writer[int(i / bin_size)].write_record(record)

    ### close all streams
    for writer in vcf_writer:
        writer.close()