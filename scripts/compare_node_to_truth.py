#EVALUATION SCRIPT
import sys
import argparse
import pandas as pd


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("df_node", help="Path to SV dataframe file")
    parser.add_argument("df_truth", help="Path to SV dataframe truth file")
    parser.add_argument("-metrics", help="Save metrics in a file")
    #parser.add_argument("-metrics", help="If dataframe should be saved to CSV, give output file path here")

    args = parser.parse_args()

    node = pd.read_csv(args.df_node)
    truth = pd.read_csv(args.df_truth)
    fp_df = []
    tp_df = []

    similar_sv_types = {'TRA': ['TRA', 'BND'],
                        'BND': ['BND', 'TRA']}

    for index, row in node.iterrows():
        #print(row)
        # Evaluate if line is in truth df
        # First get all entries from truth df with same starting position +- 1kb
        matches = ""

        #CHECK DIFFERENT FIELDS
        #POS
        matches = truth.loc[(truth['start'] > row['start'] - 1000) & (truth['start'] < row['start'] + 1000)]

        #CHROM
        matches = matches.loc[matches['start_chrom'] == row['start_chrom']]

        #LENGTH
        matches = matches.loc[(matches['length'] > row['length'] * 0.9) & (matches['length'] < row['length'] * 1.1)]

        #SV TYPE
        if row['type'] in ["TRA", "BND"]:
            pass
        else:
            matches = matches.loc[matches['type'] == row['type']]


        if matches.empty:
            # Call as false positive
            fp_df.append(row.tolist())
        if not matches.empty:
            # Call as true positive
            tp_df.append(row.tolist())
            # If complete match, remove entry from truth dataframe
            # check if matches contains more then one line
            truth = truth.drop(index=matches.index.values.tolist())
            truth = truth.reset_index(drop=True)


    print(fp_df)
    print("FP:" + str(len(fp_df)))
    print(tp_df)
    print("TP:" + str(len(tp_df)))


if __name__ == "__main__":
    main()