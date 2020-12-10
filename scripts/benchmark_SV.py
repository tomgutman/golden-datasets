#!/home/utilisateur/miniconda3/envs/happy_cyvcf2/bin/python

from cyvcf2 import VCF

def setup_var(vcf_file):
    var_list = []

    for variant in VCF(vcf_file):
        sv_dict = {}
        sv_dict['chr'] = variant.CHROM
        sv_dict['pos'] = variant.POS
        sv_dict['id'] = variant.ID
        sv_dict['ref'] = variant.REF
        sv_dict['alt'] = variant.ALT
        sv_dict['qual'] = variant.QUAL
        sv_dict['svtype'] = variant.INFO.get('SVTYPE')
        sv_dict['svclass'] = variant.INFO.get('SVCLASS')
        var_list.append(sv_dict)
    return(var_list)

def setup_truth(truth_file):
    truth_list = []

    for variant in VCF(truth_file):
        sv_dict = {}
        sv_dict['chr'] = variant.CHROM
        sv_dict['pos'] = variant.POS
        sv_dict['id'] = variant.ID
        sv_dict['ref'] = variant.REF
        sv_dict['alt'] = variant.ALT
        sv_dict['qual'] = variant.QUAL
        sv_dict['svtype'] = variant.INFO.get('SVTYPE')
        sv_dict['svlen'] = variant.INFO.get('SVLEN')
        # print(sv_dict)
        truth_list.append(sv_dict)
    return(truth_list)


if __name__=="__main__":


    var_test = setup_var('/home/utilisateur/Documents/remoteDir/Documents/Tom/EUCANCan/Benchmark/results_dream/CURIE/head_insilico_1_sv.vcf')
    truth_test = setup_truth('/home/utilisateur/Documents/remoteDir/Documents/Tom/EUCANCan/Benchmark/truth.SV.synthetic.challenge.set1.vcf')

    print(var_test[34])
    print(truth_test[34])

    for var1 in var_test:
        for var2 in truth_test:
            #print(f"{var1.get('pos')}, {var2.get('pos')}")
            #print("Var1 is {}, var2 is {}".format(var1.get('pos'), var2.get('pos')))
            if var2.get('pos') - 100 < var1.get('pos') < var2.get('pos') + 100:
                print(var1.get('pos'), var2.get('pos'))
                print(var1.get('svtype'), var2.get('svtype'))
                TP+=1
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
