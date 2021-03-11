import sys
import argparse
import ingest_sv
import os

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-sv", help="Path to SV vcf file")
    parser.add_argument("-snv", help="Sample name if multisampleVCF")
    parser.add_argument("outputdir", help="Path to directory that will contain the final results")

    args = parser.parse_args()

    outputdir = args.outputdir

    if not os.path.exists(outputdir):
        os.makedirs(outputdir)
    else:
        sys.exit("Cannot create output directory because it already exists.")


    if args.sv:
        ingest_sv.main()



if __name__ == "__main__":
    main()