import sys
import argparse
import pandas as pd
import numpy as np

#TO DO
#Make sure both indel and SV file are normalized in the same way beforehand
#Stablish proper formats (Columns), indel calling do not yields end chr and end position and sv type etc.
#Deal with real indel file, now it is tested with dummy files
#Make sure to re-evaluate the events with ref length == alt length
'''
# COLUMNS SNV: chrom,position,ref,alt
# COLUMNS SV: index,start_chrom,start,end_chrom,end,ref,alt,length,type

From INDEL/SNV to SV
# 1: end chrom = start chrom
# 2: position end = postion start + length
# 3: length = max(length(ref), length(alt))
# 4: type = del/ins/delins/other?

ref = AGC
alt = CTAATA
type = delins
length = ? > length(alt) - length(ref)

From SV to INDEL/SNV
1: chrom = start_chrom
2: position = start
3: ref = ref
4: alt = alt
'''

#TERMINAR ESTO CON LOS DUMMYS QUE HE EHCHO
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("df_indels", help="Path to indels dataframe file")
    parser.add_argument("df_sv", help="Path to SV dataframe file")
    parser.add_argument("new_df_indels", help="Path to new indels dataframe file variants <50 bp")
    parser.add_argument("new_df_sv", help="Path to new sv dataframe file variants >50 bp")   

    args = parser.parse_args()
   
    print("[INFO] ### Evaluating variants <50bp as indels and variants >=50bp as structural variants")


    df_indels = pd.read_csv(args.df_indels, dtype={'start_chrom': str, 'start': int, 'end_chrom': str, 'end': str, 'ref': str, 'alt': str, 'length': 'Int32', 'type': str})
    df_sv = pd.read_csv(args.df_sv, index_col=0, dtype={'start_chrom': str, 'start': int, 'end_chrom': str, 'end': str, 'ref': str, 'alt': str, 'length': 'Int32', 'type': str})
    print(df_indels)
    print(df_sv)


    if len(df_indels.index)>0 and len(df_sv.index)>0 : 
        threshold = 50
        #df_indels["length"] = df_indels.to_numeric(df_indels["length"])
        #df_sv["length"] = df_sv.to_numeric(df_sv["length"])

        df_indels['ref_len'] = df_indels['ref'].str.len()
        df_indels['alt_len'] = df_indels['alt'].str.len()
        df_indels['length'] = df_indels[['ref_len', 'alt_len']].max(axis=1)

        indels_longer_than_50 =  df_indels.loc[(df_indels['length'] >= threshold) | (pd.isna(df_indels['length']))] # These should be moved to SV
        sv_smaller_than_50 =  df_sv.loc[(df_sv['length'] < threshold) | (pd.isna(df_sv['length']))] # These should be moved to indels/snv
        if isinstance(indels_longer_than_50, pd.DataFrame):
            print("Is dataframe")
            indels_longer_than_50 = indels_longer_than_50.rename(columns={'chrom': 'start_chrom', 'pos': 'start'})
        else:
            print("Is series")
            indels_longer_than_50 = indels_longer_than_50.rename({'chrom': 'start_chrom', 'pos': 'start'})
        # type + end_chrom + end
        indels_longer_than_50['type'] = np.where(indels_longer_than_50['ref_len'] == indels_longer_than_50['alt_len'], 'delins', np.where(indels_longer_than_50['alt_len'] > indels_longer_than_50['ref_len'], 'INS', 'DEL'))
        indels_longer_than_50['end_chrom'] = indels_longer_than_50['start_chrom']
        indels_longer_than_50['end'] = indels_longer_than_50['start'] + indels_longer_than_50['length']
        print(indels_longer_than_50)
        print(list(indels_longer_than_50))
        print(sv_smaller_than_50)

        #Next step: sort the indels_longer_than_50 into the right order
        #Next step 2: get the correct columns of sv_smaller_than_50 df

        df_indels_true = df_indels.loc[df_indels['length'] < threshold]
        df_sv_true = df_sv.loc[df_sv['length'] >= threshold]

        df_indels = df_indels.drop(['ref_len', 'alt_len'], axis=1)

        print(df_indels_true)
        print(df_sv_true)

        # Here we first have to realign the columns
        final_indels = df_indels.append(sv_smaller_than_50, ignore_index=True)
        final_svs = df_sv.append(indels_longer_than_50, ignore_index=True)
    
        if args.new_df_indels and args.new_df_sv:
                final_indels.to_csv(args.new_df_indels)
                final_svs.to_csv(args.new_df_sv)
    
        return final_indels, final_svs

    else:
        sys.exit("No data parsed, please check your input data.")



if __name__ == "__main__":
    main()
