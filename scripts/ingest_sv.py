import sys
import argparse
import vcf
import pandas as pd


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("vcf", help="Path to SV vcf file")
    parser.add_argument("-samplename", help="Sample name if multisampleVCF")

    args = parser.parse_args()

    print("[INFO] ### Start process")
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

    nr_of_vars = 0
    nr_filtered = 0

    if hasattr(vcf_reader, 'metadata'):
        # Input file type is VCF. Start parsing.
        sample = None

        # Get reference genome version
        if 'reference' in vcf_reader.metadata:
            ref_version = vcf_reader.metadata['reference']
        else:
            ref_version = "Unknown"
        print("[INFO] Reference version: " + ref_version)

        # Get VCF format version
        if 'fileformat' in vcf_reader.metadata:
            vcf_format = vcf_reader.metadata['fileformat']
        else:
            vcf_format = "Unknown"
        print("[INFO] File format: " + vcf_format)

        # Start reading in records and add to dataframe
        for record in vcf_reader:
            nr_of_vars += 1

            # Only use 'PASS' calls
            if record.FILTER:

                # If we want to keep vars with certain filters, list here
                if record.FILTER != ['INFERRED']:
                    print('[DEBUG] Found filter: ' + str(record.FILTER))
                    nr_filtered += 1
                    continue

            # Select only the relevant sample and skip this step if performed once
            if sample == None:
                if len(record.samples) > 1:

                    # Need to filter, so need optional argument
                    if args.samplename == None:
                        sys.exit("[ERROR] Multisample VCF found and no samplename given. Please input the tumor samplename with "
                                 "-samplename to make filtering possible.")
                    else:
                        sample = args.samplename
                        print("[INFO] Filtering for the following sample in the VFC: " + sample)
                else:
                    sample = record.samples[0]
                    print("[INFO] Using only sample in the VCF: " + sample)

                # First get index of sample
                index_sample = None
                for ind, call in enumerate(record.samples):
                    if call.sample == sample:
                        index_sample = ind
                if not index_sample:
                    sys.exit("[ERROR] Cannot find sample in the VCF file. Exiting.")

            start_chrom = record.CHROM.replace("CHR", "").replace("chr", "")
            start = record.POS
            '''
            if record.POS != record.start:
                print("WEIRD")
                print(record)
                print(record.POS)
                print(record.start)
            '''

            # Check if variant is SV
            if record.is_sv:
                if len(record.ALT) > 1:
                    sys.exit("[DEBUG] MORE THAN ONE ALT ALLELE!")

                # Get end coordinate from file if registered in INFO field
                if 'END' in record.INFO:
                    end = record.sv_end + 1
                else:
                    # Parse end from ALT
                    end_pos = str(record.ALT[0]).replace("CHR", "").replace("chr", "").replace("A", "").replace("G", "").replace("T", "").replace("C", "").replace("]", "").replace("[", "")
                    if "." in end_pos:
                        '''
                        # SV is probably BND?, possible alts: ".str", "str.", "N.", ".N"
                        print("[DEBUG] IS SV BND?")
                        print(record)
                        print(record.samples)
                        print(record.INFO)
                        '''
                        end_pos = start
                        end_chrom = start_chrom
                    else:
                        #print(record)
                        end_chrom, end = end_pos.split(":")

                if 'SVLEN' in record.INFO:
                    length = record.INFO['SVLEN']
                else:
                    if end_chrom == start_chrom:
                        length = int(end) - int(start)
                    else:
                        # SV starts and ends on different chromosomes
                        length = None
                ref = record.REF
                alt = record.ALT[0]

                sv_type = None
                if 'SVTYPE' in record.INFO: # Support Hartwig
                    sv_type = record.INFO['SVTYPE']
                else:
                    sv_type = record.var_subtype


            elif record.is_snv:
                sys.exit("[DEBUG] SNV in SV file.")


            # Convert into table with: start_chrom, start, end_chrom, end, ref, alt, length, sv_type, genotype
            #print([start_chrom, start, end_chrom, end, ref, alt, length, sv_type])
            variants.append([start_chrom, start, end_chrom, end, ref, alt, length, sv_type])
    else:
        # Input file type is TSV. Start parsing.
        print("[INFO] File format: TSV file")
        ref_version = "Unknown"
        for line in vcf_reader:
            line = line.rstrip().split()
            if line[0] == "##Reference":
                print("Found reference!")
                ref_version = "\t".join(line)


        print("[INFO] Reference version: " + ref_version)


    # Create dataframe from all variant data
    columns = ["start_chrom", "start", "end_chrom", "end", "ref", "alt" ,"length", "type"]
    data = pd.DataFrame(variants, columns=columns)
    print(data)
    '''
    dummy_data = [[1, 100, 'X', 200, 'C', "[[CHRX:144760323[C]", 100, "BND"],
                  [3, 5000, 'X', 10000, 'G', "[G]CHRX:153909144]]", 5000, "BND"],
                  [5, 20, 5, 200, 'G', "[<DEL>]", 90, "DEL"]]
    dummy_df = pd.DataFrame(dummy_data, columns=columns)
    print(dummy_df)
    '''

    print("Total number of variants processed: " + str(nr_of_vars))
    print("Filtered variants of total: " + str(nr_filtered))

if __name__ == "__main__":
    main()
