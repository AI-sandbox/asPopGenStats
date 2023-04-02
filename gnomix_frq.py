import os
import numpy as np
import pandas as pd
import argparse
from collections import defaultdict
import errno
from re import match
import time
import ray


data_dir = "/home/markpenj/hawaii/gnomix_run_4pbs"


argp = argparse.ArgumentParser()
argp.add_argument("-a", "--alt_group", type=str, required=True,
                  help="alternative subpopulation group for alleles")
argp.add_argument("-r", "--ref_group", type=str, required=True,
                  help="reference subpopulation group for alleles",)
args = argp.parse_args()


# Read vcf file.
def read_skip_comment(filename, start_pattern, **kwargs):
    if os.stat(filename).st_size == 0:
        raise ValueError(f"Empty file at {filename:s}.")
    with open(filename, 'r') as f:
        line_idx = 0
        curr_line = f.readline()
        while curr_line.startswith(start_pattern):
            line_idx = f.tell()
            curr_line = f.readline()
        f.seek(line_idx)
        return pd.read_csv(f, **kwargs)


def get_famid_col_idx(col_names):
    col_idx = 0
    while (col_idx < len(col_names)) and \
          (not match("(\d|_)", col_names[col_idx])):
        col_idx += 1
    if col_idx == len(col_names):
        raise ValueError("Cannot find starting column index of samples!")
    return col_idx


def parse_subpop_code(subpop2code_str, subpop_info_):
    # parse_str = "#Subpopulation order/codes: American=0	Papuan=1	Polynesian=2	European=3	EastAsian=4	African=5"
    subpop2code_list = subpop2code_str.split(" ")[-1].split("\t")
    for subpop2code in subpop2code_list:
        subpop, subpop_code = subpop2code.split("=")
        subpop_info_[subpop] = int(subpop_code)


def get_subpop_map(subpop_info_, famid_id2pop_, chr_idx=1):
    filename = os.path.join(data_dir, f"chr{chr_idx}", "query_results.msp")
    # filename = "query_results.msp"
    if os.stat(filename).st_size == 0:
        raise ValueError(f"Empty file at {filename:s}.")
    with open(filename, 'r') as f:
        subpop_info_line = f.readline()
        parse_subpop_code(subpop_info_line, subpop_info_)

        group_names_ = set()
        sample_info_line = f.readline()
        sample_columns = sample_info_line.split('\t')
        for sample in sample_columns:
            sample_id = sample[:-2]  # remove ".1" or ".0"
            if sample_id in famid_id2pop_:
                group_names_.add(famid_id2pop_[sample_id])
        return group_names_


def check_subpop_input(subpop_info_, alt_group_, ref_group_):
    if alt_group_ not in subpop_info_:
        raise ValueError("Please provide correct name for alternative group!")
    elif ref_group_ not in subpop_info_:
        raise ValueError("Please provide correct name for reference group!")
    elif alt_group_ == ref_group_:
        raise ValueError("The alternative group and the reference group "
                         "should be different!")



def snp2int_calc(lst):
    return sum([int(i) for i in lst if type(i) == int or i.isdigit()])


def snp2int(row):
    return row.str.split("|").apply(snp2int_calc)


def get_vcf_data(chr_idx, famid_id2pop_):
    vcf_df = read_skip_comment(os.path.join(data_dir, f"chr{chr_idx}",
                                            "query_file_phased.vcf"),
                               "##", header=0, sep='\t')
    # vcf_df = read_skip_comment("query_file_phased.vcf", "##",
                            #    header=0, sep='\t')
    ## This line is commented if "POS" column in vcf file grows monotonically.
    # vcf_df.sort_values(by="POS", inplace=True)

    # Starting column index for samples, e.g., 1_4781017
    sample_sidx = get_famid_col_idx(list(vcf_df.columns))
    sample_cols = vcf_df.columns[sample_sidx:]
    # Get a dictionary "group2colidx_" with "population: [col_idx_of_sample]".
    group2colidx_ = defaultdict(list)
    for col_idx, sample in enumerate(sample_cols):
        if sample in famid_id2pop_:
            group2colidx_[famid_id2pop_[sample]].append(col_idx)
    
    # Convert SNP format to integer count.
    tmp_df = vcf_df.iloc[:, sample_sidx:]
    tmp_df = tmp_df.apply(snp2int).to_numpy(dtype=int)

    # Aggregate samples from one population group.
    aggr_alt = pd.DataFrame({pop: np.sum(tmp_df[:, colidx_list], axis=1) \
                             / (2 * len(colidx_list))
                             for pop, colidx_list in group2colidx_.items()},
                             dtype=np.float32)
    aggr_tot =  pd.DataFrame({f"{pop:s}_CT": np.ones(len(tmp_df)) \
                              * (2 * len(colidx_list))
                              for pop, colidx_list in group2colidx_.items()},
                              dtype=np.float32)
    # vcf columns: CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT
    aggr_data = pd.concat([vcf_df[["#CHROM", "POS", "ALT", "REF"]], aggr_alt, aggr_tot], axis=1)

    return aggr_data, vcf_df["POS"].to_numpy()


def get_msp_data(chr_idx, famid_id2pop_):
    msp_df = pd.read_csv(os.path.join(data_dir, f"chr{chr_idx}", "query_results.msp"),
                         header=0, skiprows=[0], sep='\t')
    # msp_df = pd.read_csv("query_results.msp", header=0, skiprows=[0], sep='\t')
    # Starting column index for samples, e.g., 1_4781017.0
    sample_sidx = get_famid_col_idx(list(msp_df.columns))
    sample_cols = msp_df.columns[sample_sidx:]
    # Get a dictionary "group2colidx_" with "population: [col_idx_of_sample]".
    group2colidx_ = defaultdict(list)
    for col_idx, sample in enumerate(sample_cols):
        sample_id = sample[:-2]  # remove ".1" or ".0"
        if sample_id in famid_id2pop_:
            group2colidx_[famid_id2pop_[sample_id]].append(col_idx)

    # Filter out the data with irrelevant populations.
    data = msp_df.iloc[:, sample_sidx:].to_numpy(dtype=int)
    data_alt = np.isclose(data, subpop_info[args.alt_group]).astype(float)
    data = data_alt + np.isclose(data, subpop_info[args.ref_group]).astype(float)

    # Aggregate samples from one population group.
    aggr_alt = pd.DataFrame({f"{pop:s}_GNX": \
                             np.divide(np.sum(data_alt[:, colidx_list], axis=1),
                                       np.sum(data[:, colidx_list], axis=1),
                                       out=np.zeros(data_alt.shape[0]),
                                       where=(np.sum(data[:, colidx_list], axis=1) != 0))
                             for pop, colidx_list in group2colidx_.items()},
                             dtype=np.float32)
    aggr_tot =  pd.DataFrame({f"{pop:s}_GNX_CT": \
                              np.sum(data[:, colidx_list], axis=1)
                              for pop, colidx_list in group2colidx_.items()},
                              dtype=np.float32)
    # msp columns: chm, spos, epos, sgpos, egpos, n snps
    aggr_data = pd.concat([aggr_alt, aggr_tot], axis=1)
    return aggr_data, msp_df["spos"].to_numpy()


def calculate_gnomix_row_repeat(pos_col_, spos_col_gnx_):
    spos2count_ = np.zeros_like(spos_col_gnx_)
    spos_col_gnx_ = np.append(spos_col_gnx_, len(pos_col_))
    prev_idx = 0
    curr_idx = 0
    curr_idx_gnx = 1
    gnx_value = spos_col_gnx_[curr_idx_gnx]
    while (curr_idx < len(pos_col_) and curr_idx_gnx < len(spos_col_gnx_) - 1):
        if (curr_idx != 0 and pos_col_[curr_idx] == gnx_value):
            spos2count_[curr_idx_gnx - 1] = curr_idx - prev_idx
            prev_idx = curr_idx
            curr_idx_gnx += 1
            gnx_value = spos_col_gnx_[curr_idx_gnx]
        curr_idx += 1
    spos2count_[curr_idx_gnx - 1] = len(pos_col_) - prev_idx
    return spos2count_

@ray.remote
def get_chrom_data(chr_idx, famid_id2pop_):
    vcf_data, pos_col       = get_vcf_data(chr_idx, famid_id2pop_)
    msp_data, spos_col_gnx  = get_msp_data(chr_idx, famid_id2pop_)

    spos2count = calculate_gnomix_row_repeat(pos_col, spos_col_gnx)
    msp_data_new = pd.DataFrame(np.repeat(msp_data.to_numpy(),
                                          repeats=spos2count, axis=0),
                                columns=msp_data.columns)
    aggr_data = pd.concat([vcf_data, msp_data_new], axis=1)
    return aggr_data


@ray.remote
def write_data(group_, chrom_data_, save_file_path_):
    include_cols = ["#CHROM", "POS", "ALT", "REF", group_, f"{group_}_CT",
                        f"{group_}_GNX", f"{group_}_GNX_CT"]
    rename_cols = {"#CHROM":"CHR",
                    "POS":"SNP",
                    "ALT":"A1",
                    "REF":"A2",
                    group_:"MAF",
                    f"{group_}_CT":"NCHROBS",
                    f"{group_}_GNX":"MAF_GNX",
                    f"{group_}_GNX_CT":"NCHROBS_GNX"}
    df_group = chrom_data_[include_cols]
    df_group = df_group.rename(columns=rename_cols)
    df_group.to_csv(f"{save_file_path_:s}/{group_:s}.freq",
                    index=False, sep='\t')


if __name__ == "__main__":
    alt_group = args.alt_group
    ref_group = args.ref_group

    # Create a dictionary "pop_info_dict" with "famid_id: pop"
    pop_info = pd.read_csv("hawaiiPopInfo.csv", header=0, usecols=["famid_id", "pop"])
    famid_id2pop = {data.values[0]: data.values[1] for _, data in pop_info.iterrows()}

    subpop_info = {}
    group_names = get_subpop_map(subpop_info, famid_id2pop)
    check_subpop_input(subpop_info, alt_group, ref_group)

    # Concatenate data from all chromosomes 1-22.
    print("Start to read vcf and msp files from each chromosome!")
    chrom_data_remote = [get_chrom_data.remote(i, famid_id2pop)
                         for i in range(1, 23)]
    chrom_data_list = ray.get(chrom_data_remote)
    chrom_data = pd.concat(chrom_data_list)
    print("Data aggregation finished!")

    # Save data to the given path.
    save_file_path = os.path.join(os.path.abspath('.'), 
                                  f"gnomix_data_{args.alt_group}_{args.ref_group}")
    try:
        os.makedirs(save_file_path)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise

    print("Start to write files!")
    write_data_remote = [write_data.remote(group, chrom_data, save_file_path)
                         for group in group_names]
    ray.get(write_data_remote)
