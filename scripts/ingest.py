import sys
import argparse
import vcf
import pandas as pd
import parse_sv_vcf
import parse_sv_tsv


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("vcf", help="Path to SV vcf file")
    parser.add_argument("-samplename", help="Sample name if multisampleVCF")

    args = parser.parse_args()

    print("[INFO] ### Start ingesting files")
    print("[INFO] file: " + str(args.vcf))
    variants = []

    # Account for: .vcf, .vcf.gz, .tsv
    if args.vcf[-3:] == ".gz":
        vcf_reader = vcf.Reader(filename=args.vcf)
    elif args.vcf[-4:] == ".vcf":
        vcf_reader = vcf.Reader(open(args.vcf, 'r'))
    elif args.vcf[-4:] == ".tsv":
        vcf_reader = open(args.vcf, 'r')
    else:
        sys.exit("[ERROR] Do not recognize type of input vcf " + args.vcf + ". Exiting.")


    if hasattr(vcf_reader, 'metadata'):
        # Input file type is VCF. Start parsing.
        variants, nr_of_vars, nr_filtered = parse_sv_vcf.parse(vcf_reader, args.samplename)

    else:
        # Input file type is TSV. Start parsing.
        print("[INFO] File format: TSV file")
        variants, nr_of_vars, nr_filtered = parse_sv_tsv.parse(vcf_reader)

    # Create dataframe from all variant data
    if variants:
        columns = ["start_chrom", "start", "end_chrom", "end", "ref", "alt" ,"length", "type"]
        data = pd.DataFrame(variants, columns=columns)
        print(data)
    else:
        sys.exit("No data parsed, please check your input data.")

    '''
    dummy_data = [[1, 100, 'X', 200, 'C', "[[CHRX:144760323[C]", 100, "BND"],
                  [3, 5000, 'X', 10000, 'G', "[G]CHRX:153909144]]", 5000, "BND"],
                  [5, 20, 5, 200, 'G', "[<DEL>]", 90, "DEL"]]
    dummy_df = pd.DataFrame(dummy_data, columns=columns)
    print(dummy_df)
    '''

    # Todo: what about fusions?

    print("Total number of variants processed: " + str(nr_of_vars))
    print("Filtered variants of total: " + str(nr_filtered))

if __name__ == "__main__":
    main()
