# Given a file directory for population files with a list of population filenames
# create an aggregated population file. If the derived allele frequency (DAF)
# is not -1 (by default), we will compute the derived allele SNP positions for
# psi statistics.

import numpy as np
import pandas as pd
import os


def construct_aggr_file(data_dir, stats_name, files, aggr_file_abspath,
                        DA_file_abspath, DAF):
    '''
    Aggregate the reference populations to a new population and save to a
    new frequency file.
    Args:
        data_dir            Absolute path to input data's directory
        stats_name          Genetic statistics name
        files               Population group filenames
        aggr_file_abspath   Absolute path to aggregated population frequency file
        DA_file_abspath     Absolute path to derived allele SNP position file
        DAF                 Derived allele frequency with range [0.0, 0.5]
    '''
    colname = ['FREQ', 'CT']
    # Load and stack the population files.
    curr_pop = pd.read_csv(os.path.join(data_dir, files[0]),
                        header=0).to_numpy()
    chrom_idx = curr_pop[:, 0]
    # data frame with column name ['FREQ_CT', 'CT']
    res = np.stack((curr_pop[:, 1] * curr_pop[:, 2], curr_pop[:, 2]), axis=-1)

    # Aggregate the reference populations to a new population and save to a
    # new frequency file.
    if len(files) > 1:
        for i in range(1, len(files)):
            curr_pop = pd.read_csv(os.path.join(data_dir, files[i]),
                                header=0)[colname].to_numpy()
            res += np.stack((curr_pop[:, 0] * curr_pop[:, 1], curr_pop[:, 1]), axis=-1)
        res[:, 0] = np.divide(res[:, 0], res[:, 1], out=np.zeros(res.shape[0]),
                            where=(res[:, 1] != 0)) # ['FREQ', 'CT']

        df_final = pd.DataFrame(np.hstack((chrom_idx.reshape(res.shape[0], 1), res)),
                                columns = ["CHROM_IDX", "FREQ", "CT"])
        df_final.to_csv(f"{aggr_file_abspath:s}", index=False)

    if stats_name == 'psi':
        if abs(DAF + 1.0) > 1e-5:
            print("Invalid derived allele frequency threshold (in range [0.0, 0.5])!")
        # (Psi only) construct a derived allele SNP position file.
        construct_DA_file(res, data_dir, DA_file_abspath, DAF)

    
def construct_DA_file(df, data_dir, DA_file_abspath, DAF):
    '''
    Construct a derived allele SNP position file based on the aggregated population.
    This file marks the derivied allele SNP loci using integers (0, 1, 2) in
    the following format:
    |__________|____________________________|___________|
    0         DAF                         1-DAF         1
    ---- 1 ---- ------------ 0 ------------- ---- 2 ----
    where DAF represents the derived allele frequency threshold [0.0 - 0.5] for
    the reference aggregated population.
    Args:
        df                  Numpy dataframe for aggregated population file with
                                column name ["FREQ", "CT"]
        data_dir            Absolute path to input data's directory
        DA_file_abspath     Absolute path to derived allele SNP position file
        DAF                 Derived allele frequency with range [0.0, 0.5]
    '''
    if abs(DAF - 0.5) <= 1e-5:
        drv_frq = np.ones(df.shape[0])
        stream_output = [f"No Polarization. Number of SNPs: {df.shape[0]}"]
    else:
        drv_frq = 1 * (df[:, 0] <= DAF)
        drv_frq += 2 * (df[:, 0] >= (1-DAF))
        stream_output = []
        stream_output.append(f"Number of SNPs: {df.shape[0]}")
        stream_output.append(f"Number of DAF {DAF} (TYPE 1): "
                             f"{np.sum(df[:, 0] <= DAF)}")
        stream_output.append(f"Number of DAF {1-DAF} (TYPE 2): "
                             f"{np.sum(df[:, 0] >= (1-DAF))}")
    for i in range(len(stream_output)):
        print(stream_output[i])
    DAF_str = str(DAF).split(".")[-1][:2] # 0.05 -> "05"
    DAF_str += "0" * (2 - len(DAF_str)) # 0.0 -> "00"
    np.savetxt(f"{data_dir:s}/{DA_file_abspath:s}", drv_frq, fmt='%d',
               delimiter='\n', header=". ".join(stream_output))
