#!/usr/bin/env python
# -*- coding: utf-8 -*-

from sys import argv
from cyvcf2 import VCF


def setup_var(vcf_file):
    """

    Args:
        vcf_file:

    Returns:

    """
    var_list = []

    for variant in VCF(vcf_file):
        sv_dict = {'chr': variant.CHROM, 'pos': variant.POS, 'id': variant.ID, 'ref': variant.REF, 'alt': variant.ALT,
                   'qual': variant.QUAL, 'svtype': variant.INFO.get('SVTYPE'), 'svclass': variant.INFO.get('SVCLASS')}
        var_list.append(sv_dict)
    return var_list


def setup_truth(truth_file):
    """

    Args:
        truth_file:

    Returns:

    """
    truth_list = []

    for variant in VCF(truth_file):
        sv_dict = {'chr': variant.CHROM, 'pos': variant.POS, 'id': variant.ID, 'ref': variant.REF, 'alt': variant.ALT,
                   'qual': variant.QUAL, 'svtype': variant.INFO.get('SVTYPE'), 'svlen': variant.INFO.get('SVLEN')}
        # print(sv_dict)
        truth_list.append(sv_dict)
    return truth_list


def main(args):
    """

    Args:
        args ([str]): command line parameter list

    Returns:

    """
    var_test, truth_test = args
    print(var_test[34])
    print(truth_test[34])

    for var1 in var_test:
        for var2 in truth_test:
            # print(f"{var1.get('pos')}, {var2.get('pos')}")
            # print("Var1 is {}, var2 is {}".format(var1.get('pos'), var2.get('pos')))
            if var2.get('pos') - 100 < var1.get('pos') < var2.get('pos') + 100:
                print(var1.get('pos'), var2.get('pos'))
                print(var1.get('svtype'), var2.get('svtype'))
                TP += 1
        break

    # for variant in VCF('/data/tmp/tgutman/SeqOIA/TMB/insilico_3_sv_gridss.vcf'):
    #     ref = variant.REF
    #     alt = variant.ALT
    #     pos = variant.POS
    #     chr = variant.CHROM
    #     asc = variant.INFO.get('ASC')
    #     assr = variant.format('ASSR')
    #     id = variant.ID
    #     filter = variant.FILTER
    #     print(chr, pos, ref, alt, filter, asc, assr, id)
    #     break


if __name__ == "__main__":

    main(argv[1:])
