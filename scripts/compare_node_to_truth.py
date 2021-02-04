#EVALUATION SCRIPT
import sys
import argparse
import pandas as pd
import numpy as np

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

    if pd.isnull(truth['length']) or pd.isnull(test['length']):
        diff_length = np.NaN
    else:
        diff_length = test['length'] - truth['length']

    if pd.isnull(truth['length']):
        length_truth = np.NaN
        norm_start_pos = np.NaN
        norm_end_pos = np.NaN
    else:
        length_truth = truth['length']
        norm_start_pos = (int(test['start']) - int(truth['start'])) / truth['length']
        norm_end_pos = (int(test['end']) - int(truth['end'])) / truth['length']

    if pd.isnull(truth['length']) or pd.isnull(test['length']):
        length_ratio = np.NaN
    else:
        length_ratio = min(test['length'], truth['length']) / max(test['length'], truth['length'])

    # Compare the types of the SVs
    #print(test['type'] + "\t" + truth['type'])
    if test['type'] == truth['type']:
        match_type = "YES"
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
        print("We should check, type mismatch?: ")
        print(truth)
        print(test)
        print(test['type'] + "\t" + truth['type'])
        match_type = "NO"

    if truth['times_checked'] < 1:
        dup_truth = False
    else:
        dup_truth = True

    # Debugging
    #print([same_chrom_start, same_chrom_end, var_in_truth_within_window, diff_start_pos, diff_end_pos, diff_length, length_truth, norm_start_pos, norm_end_pos, length_ratio, match_type, dup_truth])

    return([same_chrom_start, same_chrom_end, var_in_truth_within_window, diff_start_pos, diff_end_pos, diff_length, length_truth, norm_start_pos, norm_end_pos, length_ratio, match_type, dup_truth])

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("df_node", help="Path to SV dataframe test file")
    parser.add_argument("df_truth", help="Path to SV dataframe truth file")
    parser.add_argument("-metrics", help="Save metrics in a file")
    #parser.add_argument("-metrics", help="If dataframe should be saved to CSV, give output file path here")

    args = parser.parse_args()

    node = pd.read_csv(args.df_node, index_col=0, dtype={'start_chrom': str, 'start': int, 'end_chrom': str, 'end': str, 'ref': str, 'alt': str, 'length': 'Int32', 'type': str})
    truth = pd.read_csv(args.df_truth, index_col=0, dtype={'start_chrom': str, 'start': int, 'end_chrom': str, 'end': str, 'ref': str, 'alt': str, 'length': 'Int32', 'type': str})
    truth['times_checked'] = 0
    sv_comp_list = []
    sv_fp_list = []
    columns_sv_comp_list = ['same_chrom_start', 'same_chrom_end', 'var_in_truth_within_window', 'diff_start_pos', 'diff_end_pos', 'diff_length', 'length_truth', 'norm_start_pos', 'norm_end_pos', 'length_ratio', 'match_type', 'dup_truth']
    window = 1000

    # For each row in the test file:
    for index, row in node.iterrows():
        #print(row)
        # Check if there is a variant in truth that matches start_chr and end_chr and start position within 1kb
        truth_matches = ""
        truth_matches = truth.loc[(truth['start_chrom'] == row['start_chrom']) & (truth['end_chrom'] == row['end_chrom']) \
                                  & (abs(truth['start'] - row['start']) <  window)]

        if not truth_matches.empty:
            #print("[DEBUG] Match found.")
            #print(truth_matches)
            # Add matches to truth_matched_df to evaluate if we've seen them before
            # Find the best match and drop the other one
            best = None
            for j, truth_match in enumerate(truth_matches.values.tolist()):
                # Start processing singles
                # Get the corresponding row from the dataframe to process below
                results = calculate_characteristics(truth_matches.iloc[j], row, window)
                #Compare the test variant (row) to every possible truth variant within the window. Then we will select the most appropiate one based on start site or overlapping rules...
                if not best:
                    best = [j, truth_match, results]
                else:
                    #evaluate if this entry is better: DEFINITION: start site diff is smallest
                    if abs(results[3]) < abs(best[2][3]):
                        print("[DEBUG] NEW BEST RESULT")
                        best = [j, truth_match, results]
            sv_comp_list.append(best[2])
            truth.loc[truth_matches.iloc[best[0]].name,'times_checked'] += 1 #UPDATE
            #print(truth_matched_df)
        else:
            #print("[DEBUG] No match found")
            sv_fp_list.append(list(row))
        #print(row)
        #print(matches)

    sv_comp_df = pd.DataFrame(sv_comp_list, columns=columns_sv_comp_list)
    sv_fp_df = pd.DataFrame(sv_fp_list, columns=list(node.columns))
    sv_fn_df = truth.loc[truth['times_checked'] == 0]

    print(sv_comp_df)
    print(sv_fp_df)
    print(sv_fn_df)

    #TODO: get the correct TP value
    recall = sv_comp_df.shape[0] / (sv_comp_df.shape[0] + sv_fn_df.shape[0])
    precision = sv_comp_df.shape[0]/(sv_comp_df.shape[0] + sv_fp_df.shape[0])
    f1 = 2*(recall * precision)/(recall + precision)

    print("Recall:\t" + str(round(recall,2)))
    print("Precision:\t" + str(round(precision,2)))
    print("F1-score:\t" + str(round(f1,2)))

if __name__ == "__main__":
    main()