# Import npz file of masked matrix and write a frequency file '.freq'
# for each nation.

import numpy as np
import pandas as pd
import os
import errno
import argparse
import textwrap
import ray
from collections import defaultdict

argp = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
argp.add_argument('-t', '--data_dir', required=True,
                  help="(REQUIRED) maasMDS data file directory in 'npz'")
argp.add_argument('-s', '--saving_dir', required=True,
                  help="(REQUIRED) saving directory for output files")
argp.add_argument('-o', '--output_dir', required=False, default="data_freq",
                  help=textwrap.dedent('''\
                    a new directory under the saving directory created for
                    all frequency files'''))
argp.add_argument('-a', '--ancestry', required=True, type=int,
                  help=textwrap.dedent('''\
                  (REQUIRED) target ancestry number (int)
                  NOTE: Please check the log file in 'txt' under the same
                  directory and find the parameter value for 'ANCESTRY'.
                  Or check any msp file in the related Gnomix run output
                  files and search for the subpopulation order/code for
                  the target ancestry.'''))
args = argp.parse_args()

data_dir = args.data_dir
saving_dir = args.saving_dir
output_dir = args.output_dir
ancestry = str(args.ancestry)

data = np.load(data_dir, allow_pickle=True)

# Create a directory (if not exist) for the data.
save_file_path = os.path.join(saving_dir, output_dir)
try:
    os.makedirs(save_file_path)
except OSError as e:
    if e.errno != errno.EEXIST:
        raise

# Entries in data have column names in the following format: 
# ['masks', 'rs_ID_list', 'ind_ID_list', 'labels', 'weights']
# # Use the following 2 lines to print "keys" in NpzFile (data.files).
# for k in data.files:
#     print(k)

rs_ID_list = np.sort(data['rs_ID_list']) # SNP loci
ind_ID_list = data['ind_ID_list'] # individual
labels = data['labels']
df = pd.DataFrame(data['masks'][0][ancestry], index=list(rs_ID_list),
                  columns=list(ind_ID_list))
df = df.sort_index()
print(df.head(5))

# Create a map of nations with respect to column indices.
nation_dict = defaultdict(list)
for j in range(len(labels)):
    nation_dict[labels[j]].append(j)

@ray.remote
def create_freq(nation, nation_dict, df, rs_ID_list, save_file_path):
    '''
    Create a frequency file (.freq) for one specific nation.
    Args:
        nation                  A nation
        nation_dict             The nation map with key-value pairs
                                    (nation, [nation_sample_column_indices])
        df                      A dataframe for SNP x individual table
        rs_ID_list              A list of SNP IDs
        save_file_path          The saving path for new frequency file
    '''
    df_nation = df.iloc[:, nation_dict[nation]]
    df_sum = np.sum(np.nan_to_num(df_nation.values, nan=0), axis=1)
    df_count = np.sum(~np.isnan(df_nation.values), axis=1)
    
    # df_sum[i] <= df_count[i] (as df_sum[i] \in {0, 1})
    # Trivial Case:
    # When df_count = 0, df_sum = 0. To avoid division-by-zero error,
    # we change the entries of df_count from 0 to 1。
    df_freq = np.divide(df_sum, df_count, np.zeros(df_sum.shape),
                        where=(df_count != 0))
    np_data = np.vstack((rs_ID_list, df_freq, df_count)).T
    df_final = pd.DataFrame(np_data, columns = ["CHROM_IDX", "FREQ", "CT"])
    df_final.to_csv(os.path.join(save_file_path, f"{nation:s}.freq"),
                    index = False)

# Create frequency files for all nations.
remote_elements = [create_freq.remote(nation, nation_dict, df, rs_ID_list,
                                      save_file_path) for nation in nation_dict.keys()]
ray.get(remote_elements)

# Create a text file including all nations.
np.savetxt("pop_list.txt", sorted(list(nation_dict.keys())), fmt='%s')
