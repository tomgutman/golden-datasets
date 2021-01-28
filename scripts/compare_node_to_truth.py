#EVALUATION SCRIPT
import sys
import argparse
import pandas as pd

'''
inv_a <-VC_tab_formated[grepl("*]$",VC_tab_formated$ALT),]
      inv_a$SV_TYPE <- "INV"
      inv_b <- VC_tab_formated[grepl("^\\Q[\\E",VC_tab_formated$ALT),]
      inv_b$SV_TYPE <- "INV"
      
      dup_b <- VC_tab_formated[grepl("^\\Q]\\E",VC_tab_formated$ALT),]
      dup_b$SV_TYPE <- "DUP"
      
      del_a <- VC_tab_formated[grepl("\\Q[\\E$",VC_tab_formated$ALT),]
      del_a$SV_TYPE <- "DEL"
'''

def calculate_characteristics(truth, test, window):
    """Return a list that contains information about the comparison between truth and test.

    """
    #If two truth events are being compared to one test variant, decide here which one are we selecting based on start site or overlapping rules...

    if truth['start_chrom'] == test['start_chrom']:
        same_chrom_start = True
    else:
        same_chrom_start = False

    if truth['end_chrom'] == test['end_chrom']:
        same_chrom_end = True
    else:
        same_chrom_end = False

    if abs(truth['start'] - test['start']) <  window:
        var_in_truth_within_window = True
    else:
        var_in_truth_within_window = False

    diff_start_pos = int(test['start']) - int(truth['start'])

    diff_end_pos = int(test['end']) - int(truth['end'])

    diff_length = test['length'] - truth['length']

    norm_start_pos = (int(test['start']) - int(truth['start'])) / truth['length']

    norm_end_pos = (int(test['end']) - int(truth['end'])) / truth['length']

    length_ratio = min(test['length'], truth['length']) / max(test['length'], truth['length'])

    print(test['type'] + "\t" + truth['type'])
    if test['type'] == truth['type']:
        match_type = "yes"
    elif test['type'] == "BND":
        match_type = "BND_test"
    elif truth['type'] == "BND":
        match_type = "BND_truth"
    elif truth['type'] == "DUP" and test['type'] == "INS":
        match_type = "similar"
    elif truth['type'] == "INS" and test['type'] == "DUP":
        match_type = "similar"
    elif test['type'] == "tandem-DUP" and truth['type'] == "DUP":
        match_type = "similar"
    else:
        #Make it more complete in the near future... similar, BND test, BND truth etc....
        print("We should check: ")
        print(truth)
        print(test)
        print(test['type'] + "\t" + truth['type'])
        match_type = "NO"


    dup_truth = truth['times_checked']

    print([same_chrom_start, same_chrom_end, var_in_truth_within_window, diff_start_pos, diff_end_pos, diff_length, norm_start_pos, norm_end_pos, length_ratio, match_type, dup_truth])

    return([same_chrom_start, same_chrom_end, var_in_truth_within_window, diff_start_pos, diff_end_pos, diff_length, norm_start_pos, norm_end_pos, length_ratio, match_type, dup_truth])

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("df_node", help="Path to SV dataframe test file")
    parser.add_argument("df_truth", help="Path to SV dataframe truth file")
    parser.add_argument("-metrics", help="Save metrics in a file")
    #parser.add_argument("-metrics", help="If dataframe should be saved to CSV, give output file path here")

    args = parser.parse_args()

    node = pd.read_csv(args.df_node)
    truth = pd.read_csv(args.df_truth)
    truth['times_checked'] = 0
    sv_df = []
    columns_sv_df = ['same_chrom_start', 'same_chrom_end', 'var_in_truth_within_window', 'diff_start_pos', 'diff_end_pos', 'diff_length', 'norm_start_pos', 'norm_end_pos', 'length_ratio', 'match_type', 'dup_truth']
    truth_matched_df = []
    window = 1000

    similar_sv_types = {'TRA': ['TRA', 'BND'],
                        'BND': ['BND', 'TRA']}

    # For each row in the test file:
    for index, row in node.iterrows():
        #print(row)
        # Check if there is a variant in truth that matches start_chr and end_chr and start position within 1kb
        matches = ""
        matches = truth.loc[(truth['start_chrom'] == row['start_chrom']) & (truth['end_chrom'] == row['end_chrom']) \
                            & (abs(truth['start'] - row['start']) <  window)]

        if not matches.empty:
            # Add matches to truth_matched_df to evaluate if we've seen them before
            if matches.shape[0] > 1:	# Skip over the rows with multiple matches in truth
                continue
            for j, match in enumerate(matches.values.tolist()):
                # Start processing singles
                if match not in truth_matched_df: #Here it is comparing every field in matches, it should only compare the ones that are used to check if the test variant "has a truth one" some lines above
                    #print(match)
                    # Get the corresponding row from the dataframe to process below
                    sv_df.append(calculate_characteristics(matches.iloc[j], row, window))
                    #Compare the test variant (row) to every possible truth variant within the window. Then we will select the most appropiate one based on start site or overlapping rules...
                    truth_matched_df.append(match)
                    ## TODO: increase times_checked number of this truth thingy
                else:
                    print("#### We've seen this one before")
                    #Add a TRUE to the column "dup_variant" of the test dataframe
        #print(truth_matched_df)
        #print(row)
        #print(matches)



    '''
    #CHECK DIFFERENT FIELDS
    #POS
    #matches = truth.loc[(truth['start'] > row['start'] - 1000) & (truth['start'] < row['start'] + 1000)]

    #CHROM
    #matches = matches.loc[matches['start_chrom'] == row['start_chrom']]

    #LENGTH
    #matches = matches.loc[(matches['length'] > row['length'] * 0.9) & (matches['length'] < row['length'] * 1.1)]

    #SV TYPE
    #if row['type'] in ["TRA", "BND"]:
    #    pass
    #else:
    #    matches = matches.loc[matches['type'] == row['type']]


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


'''

if __name__ == "__main__":
    main()