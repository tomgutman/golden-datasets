import sys

def parse(vcf_reader, samplename):
    """
    Function to parse SNV vcf file (call or truth file)
    require the vcf path and the samplename if multisample vcf
    """

    variants_snv = []
    variants_indel = []
    sample = None
    nr_of_vars = 0
    nr_filtered = 0

    for record in vcf_reader:
        #print(record.var_type)
        nr_of_vars += 1

        # Only use 'PASS' calls
        if record.FILTER == []:

            # Select only the relevant sample and skip this step if performed once
            if sample == None:
                if len(record.samples) > 1:

                    # Need to filter, so need optional argument
                    if samplename == None:
                        sys.exit("[ERROR] Multisample VCF found and no samplename given. Please input the tumor samplename with "
                                 "-samplename to make filtering possible.")
                    else:
                        sample = samplename
                        print("[INFO] Filtering for the following sample in the VFC: " + sample)
                else:
                    sample = str(record.samples[0])
                    print("[INFO] Using only sample in the VCF: " + sample)

            chrom = record.CHROM
            pos = record.POS
            ref = record.REF
            alt = record.ALT

            if record.var_type == "indel":
                variants_indel.append([chrom, pos, ref, alt])
                #print("indel: {} {} {}".format(pos,ref,alt))

            elif record.var_type == "snp":
                #print ("snp: {} {} {}".format(pos,ref,alt))
                variants_snv.append([chrom, pos, ref, alt])
            else:
                print("unknown: {} {} {} {}".format(record.var_type,pos,ref,alt))



        elif record.FILTER != []:
            nr_filtered += 1
            continue

    #print(variants)

    return(variants_snv, variants_indel, nr_of_vars, nr_filtered)
