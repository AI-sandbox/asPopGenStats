# Given a file directory for population files with a list of population filenames
# create an aggregated population file. If the derived allele frequency (DAF)
# is not -1 (by default), we will compute the derived allele SNP positions for
# psi statistics.

import numpy as np
import pandas as pd
import sys
import os

if len(sys.argv) < 4:
    # We accept multple population files (.freq) with the last argument
    # as the derived allele frequency (DAF in float/int).
    # If we want to make aggregated outgroup population for F3 statistics,
    # without a given DAF, please use -1 as default value for DAF.
    print("Usage: python3 pop_aggr.py <data_directory> <population_file.freq> ... "
          "<derived_allele_frequency (float/int)> <aggreated_population_filename> "
          "<derived_allele_frequency_output_filename>")
    sys.exit(0)

data_dir = sys.argv[1].rstrip('/') # e.g., "data/cluster2"
files = sys.argv[2:-3]
DAF = float(sys.argv[-3])
aggr_filename = sys.argv[-2]
DA_filename = sys.argv[-1]

colname = ['FREQ', 'CT']
# Aggregation the population files.
curr_pop = pd.read_csv(os.path.join(data_dir, files[0]),
                       header=0).to_numpy()
chrom_idx = curr_pop[:, 0]
res = np.stack((curr_pop[:, 1] * curr_pop[:, 2], curr_pop[:, 2]), axis=-1) # ['FREQ_CT', 'CT']

for i in range(1, len(files)):
    curr_pop = pd.read_csv(os.path.join(data_dir, files[i]),
                           header=0)[colname].to_numpy()
    res += np.stack((curr_pop[:, 0] * curr_pop[:, 1], curr_pop[:, 1]), axis=-1)
res[:, 0] = np.divide(res[:, 0], res[:, 1], out=np.zeros(res.shape[0]), where=(res[:, 1] != 0)) # ['FREQ', 'CT']

df_final = pd.DataFrame(np.hstack((chrom_idx.reshape(res.shape[0], 1), res)),
                        columns = ["CHROM_IDX", "FREQ", "CT"])
df_final.to_csv(f"{data_dir:s}/{aggr_filename:s}", index=False)

# Prepare a derived allele SNP position file.
if DAF != -1.0:
    if abs(DAF - 0.5) <= 1e-5:
        drv_frq = np.ones(res.shape[0])
        stream_output = [f"No Polarization. Number of SNPs: {res.shape[0]}"]
    else:
        drv_frq = 1 * (res[:, 0] <= DAF)
        drv_frq += 2 * (res[:, 0] >= (1-DAF))
        stream_output = []
        stream_output.append(f"Number of SNPs: {res.shape[0]}")
        stream_output.append(f"Number of DAF {DAF} (TYPE 1): {np.sum(res[:, 0] <= DAF)}")
        stream_output.append(f"Number of DAF {1-DAF} (TYPE 2): {np.sum(res[:, 0] >= (1-DAF))}")
    for i in range(len(stream_output)):
        print(stream_output[i])
    DAF_str = str(DAF).split(".")[-1][:2] # 0.05 -> "05"
    DAF_str += "0" * (2 - len(DAF_str)) # 0.0 -> "00"
    np.savetxt(f"{data_dir:s}/{DA_filename:s}",
               drv_frq, fmt = "%d", delimiter = "\n",
               header = ". ".join(stream_output))
