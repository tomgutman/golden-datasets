import sys
import argparse
import vcf


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("vcf", help="Path to SV vcf file")

    args = parser.parse_args()

    print("### Start process")
    # Account for: .vcf, .vcf.gz, .tsv
    if args.vcf[-3:] == ".gz":
        vcf_reader = vcf.Reader(filename=args.vcf)
    elif args.vcf[-4:] == ".vcf":
        vcf_reader = vcf.Reader(open(args.vcf, 'r'))
    elif args.vcf[-4:] == ".tsv":
        vcf_reader = open(args.vcf, 'r')
    else:
        sys.exit("Error: do not recognize type of input vcf " + args.vcf + ". Exiting.")

    if hasattr(vcf_reader, 'metadata'):
        # Input file type is VCF. Start parsing.

        # Get reference genome version
        if 'reference' in vcf_reader.metadata:
            ref_version = vcf_reader.metadata['reference']
        else:
            ref_version = "Unknown"
        print("Reference version: " + ref_version)

        # Get VCF format version
        if 'fileformat' in vcf_reader.metadata:
            vcf_format = vcf_reader.metadata['fileformat']
        else:
            vcf_format = "Unknown"
        print("VCF format: " + "vcf_format")

        # Start reading in records and add to dataframe
        for record in vcf_reader:
            print(record)

            # Convert into table with: chr, start-pos, stop-pos, length, type
            sys.exit()
    else:
        # Input file type is TSV. Start parsing.
        for line in vcf_reader:
            line = line.rstrip().split()
            if line[0] == "##Reference":
                print("Found reference!")
                print(line)


        #print("Reference version: " + ref_version)
        print("VCF format: TSV file")



if __name__ == "__main__":
    main()
