import sys
import argparse
import pandas as pd

#TO DO
#Make sure both indel and SV file are normalized in the same way beforehand
#Stablish proper formats (Columns), indel calling do not yields end chr and end position and sv type etc.
#Deal with real indel file, now it is tested with dummy files

#TERMINAR ESTO CON LOS DUMMYS QUE HE EHCHO
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("df_indels", help="Path to indels dataframe file")
    parser.add_argument("df_sv", help="Path to SV dataframe file")
    parser.add_argument("new_df_indels", help="Path to new indels dataframe file variants <50 bp")
    parser.add_argument("new_df_sv", help="Path to new sv dataframe file variants >50 bp")   

    args = parser.parse_args()
   
    print("[INFO] ### Evaluating variants <50bp as indels and variants >50bp as structural variants")

    
    df_indels = pd.read_csv(args.df_indels, index_col=0, dtype={'start_chrom': str, 'start': int, 'end_chrom': str, 'end': str, 'ref': str, 'alt': str, 'length': 'Int32', 'type': str})
    df_sv = pd.read_csv(args.df_sv, index_col=0, dtype={'start_chrom': str, 'start': int, 'end_chrom': str, 'end': str, 'ref': str, 'alt': str, 'length': 'Int32', 'type': str})

    
    if len(df_indels.index)>0 and len(df_sv.index)>0 : 
    
        #df_indels["length"] = df_indels.to_numeric(df_indels["length"])
        #df_sv["length"] = df_sv.to_numeric(df_sv["length"])
        
        indels_longer_than_50 =  df_indels['length']>50
        sv_smaller_than_50 =  df_sv['length']<50
    
        df_indels = df_indels[df_indels.length < 50]
        df_sv = df_sv[df_sv.length > 50]
    
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
