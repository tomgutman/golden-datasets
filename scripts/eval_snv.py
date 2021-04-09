import sys
import argparse
import vcf
import pandas as pd
import csv

# Load data
print("test")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-a", "--test_df", type=str, help="Csv file to compare to truth csv (.csv)")
    parser.add_argument("-b", "--truth_df", type=str, help="Csv of the truth csv (.csv)")
    parser.add_argument("-o", "--outputfile", type=str, help="Output file name")

    args = parser.parse_args()
    print(args)
    #return (args)

# Combine ref_alt

# compare two dataframes

# Compute metrics

if __name__ == "__main__":
    main()
    #print(args)

    #test_df = pd.read_csv(args.test_df)
    #truth_df = pd.read_csv(args.truth_df)
    #merge_df = pd.merge(test_df,truth_df, on 'chrom')

    #TP = truth_df[~truth_df.isin(test_df)].dropna()
    #print(test_df, truth_df)
