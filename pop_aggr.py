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
    # ["CHR", "SNP", "A1", "A2", "MAF", "NCHROBS", "MAF_GNX", "NCHROBS_GNX"]
    colname = ["MAF", "NCHROBS", "MAF_GNX", "NCHROBS_GNX"]
    colname_chrom = ["CHR", "SNP", "A1", "A2"]
    # Load and stack the population files.
    curr_pop = pd.read_csv(os.path.join(data_dir, files[0]),
                           header=0, sep='\t')
    chrom_idx = curr_pop[colname_chrom]
    # data frame with column name ["MAF", "NCHROBS", "MAF_GNX", "NCHROBS_GNX"]
    curr_pop_stat = curr_pop[colname].to_numpy()
    res = np.stack((curr_pop_stat[:, 0] * curr_pop_stat[:, 1], curr_pop_stat[:, 1],
                    curr_pop_stat[:, 2] * curr_pop_stat[:, 3], curr_pop_stat[:, 3]),
                    axis=-1)

    # Aggregate the reference populations to a new population and save to a
    # new frequency file.
    if len(files) > 1:
        for i in range(1, len(files)):
            curr_pop_stat = pd.read_csv(os.path.join(data_dir, files[i]),
                                        header=0, sep='\t')[colname].to_numpy()
            res += np.stack((curr_pop_stat[:, 0] * curr_pop_stat[:, 1],
                             curr_pop_stat[:, 1],
                             curr_pop_stat[:, 2] * curr_pop_stat[:, 3],
                             curr_pop_stat[:, 3]),
                             axis=-1)
        res[:, 0] = np.divide(res[:, 0], res[:, 1], out=np.zeros(res.shape[0]),
                              where=(res[:, 1] != 0))  # ["MAF", "NCHROBS"]
        res[:, 2] = np.divide(res[:, 2], res[:, 3], out=np.zeros(res.shape[0]),
                              where=(res[:, 3] != 0))  # ["MAF_GNX", "NCHROBS_GNX"]
        res = res.astype('float32')

        df_final = pd.DataFrame(np.hstack((chrom_idx, res)),
                                columns=["CHR", "SNP", "A1", "A2", "MAF",
                                         "NCHROBS", "MAF_GNX", "NCHROBS_GNX"])
        df_final.to_csv(f"{aggr_file_abspath:s}", index=False)

    if stats_name == 'psi':
        if abs(DAF - 0.25) > 0.25:
            print("Invalid derived allele frequency threshold (in range [0.0, 0.5])!")
        # (Psi only) construct a derived allele SNP position file.
        construct_DA_file(res[:, :2], data_dir, DA_file_abspath, DAF)

    
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
                                column name ["MAF", "NCHROBS"]
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
    np.savetxt(f"{DA_file_abspath:s}", drv_frq, fmt='%d',
               delimiter='\n', header=". ".join(stream_output))
