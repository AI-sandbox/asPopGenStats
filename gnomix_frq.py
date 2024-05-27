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
argp.add_argument(
    "-m",
    "--mask_group",
    type=str,
    default="Polynesian",
    help="mask subpopulation group for alleles",
)
args = argp.parse_args()


# Read vcf file.
def read_skip_comment(filename, start_pattern, **kwargs):
    if os.stat(filename).st_size == 0:
        raise ValueError(f"Empty file at {filename:s}.")
    with open(filename, "r") as f:
        line_idx = 0
        curr_line = f.readline()
        while curr_line.startswith(start_pattern):
            line_idx = f.tell()
            curr_line = f.readline()
        f.seek(line_idx)
        return pd.read_csv(f, **kwargs)


def get_famid_col_idx(col_names):
    col_idx = 0
    while (col_idx < len(col_names)) and (
        match("^#?[A-Za-z]+( [A-Za-z]+)?$", col_names[col_idx])
    ):
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
    with open(filename, "r") as f:
        subpop_info_line = f.readline()
        parse_subpop_code(subpop_info_line, subpop_info_)

        group_names_ = set()
        sample_info_line = f.readline()
        sample_columns = sample_info_line.split("\t")
        for sample in sample_columns:
            sample_id = sample[:-2]  # remove ".1" or ".0"
            if sample_id in famid_id2pop_:
                group_names_.add(famid_id2pop_[sample_id])
        return group_names_


def check_subpop_input(subpop_info_, mask_group):
    if mask_group not in subpop_info_:
        raise ValueError("Please provide correct name for mask group!")


def snp2int(col):
    col_df = col.str.split("|", expand=True)
    col_df.columns = [f"{col.name}.0", f"{col.name}.1"]
    return col_df


def get_vcf_data(chr_idx, famid_id2pop_):
    vcf_df = read_skip_comment(
        os.path.join(data_dir, f"chr{chr_idx}", "query_file_phased.vcf"),
        "##",
        header=0,
        sep="\t",
    )
    # vcf_df = read_skip_comment("query_file_phased.vcf", "##",
    #    header=0, sep='\t')
    ## This line is commented if "POS" column in vcf file grows monotonically.
    # vcf_df.sort_values(by="POS", inplace=True)

    # Starting column index for samples, e.g., 1_4781017
    sample_sidx = get_famid_col_idx(list(vcf_df.columns))

    # Convert SNP format to integer count.
    data_df = vcf_df.iloc[:, sample_sidx:]
    data_df = pd.concat([snp2int(data_df[col]) for col in data_df.columns], axis=1)
    data_df = data_df[sorted(data_df.columns)]
    sample_cols = data_df.columns
    data = data_df.to_numpy(dtype=int)

    # Get a dictionary "group2colidx" with "population: [col_idx_of_sample]".
    group2colidx = defaultdict(list)
    for col_idx, sample in enumerate(sample_cols):
        sample_id = sample[:-2]  # remove ".1" or ".0"
        if sample_id in famid_id2pop_:
            group2colidx[famid_id2pop_[sample_id]].append(col_idx)

    # Aggregate samples from one population group.
    aggr_alt = pd.DataFrame(
        {
            pop: np.sum(data[:, colidx_list], axis=1) / len(colidx_list)
            for pop, colidx_list in group2colidx.items()
        },
        dtype=np.float32,
    )
    aggr_tot = pd.DataFrame(
        {
            f"{pop:s}_CT": np.ones(len(data)) * len(colidx_list)
            for pop, colidx_list in group2colidx.items()
        },
        dtype=np.float32,
    )
    # vcf columns: CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT
    aggr_data = pd.concat(
        [vcf_df[["#CHROM", "POS", "ALT", "REF"]], aggr_alt, aggr_tot], axis=1
    )

    return aggr_data, data_df, vcf_df["POS"].to_numpy()


def get_msp_data(chr_idx, famid_id2pop_):
    msp_df = pd.read_csv(
        os.path.join(data_dir, f"chr{chr_idx}", "query_results.msp"),
        header=0,
        skiprows=[0],
        sep="\t",
    )
    # msp_df = pd.read_csv("query_results.msp", header=0, skiprows=[0], sep='\t')
    # Starting column index for samples, e.g., 1_4781017.0
    sample_sidx = get_famid_col_idx(list(msp_df.columns))

    # Filter out the data with irrelevant populations.
    data_df = msp_df.iloc[:, sample_sidx:]
    data_df = data_df[sorted(data_df.columns)]
    sample_cols = data_df.columns
    data = data_df.to_numpy(dtype=int)
    data_mask = data != subpop_info[args.mask_group]

    return [data_mask, sample_cols], msp_df["spos"].to_numpy()


def calculate_gnomix_row_repeat(pos_col, spos_col_msk):
    spos2count = np.zeros_like(spos_col_msk)
    spos_col_msk = np.append(spos_col_msk, len(pos_col))
    prev_idx = 0
    curr_idx = 0
    curr_idx_msk = 1
    msk_value = spos_col_msk[curr_idx_msk]
    while curr_idx < len(pos_col) and curr_idx_msk < len(spos_col_msk) - 1:
        if curr_idx != 0 and pos_col[curr_idx] == msk_value:
            spos2count[curr_idx_msk - 1] = curr_idx - prev_idx
            prev_idx = curr_idx
            curr_idx_msk += 1
            msk_value = spos_col_msk[curr_idx_msk]
        curr_idx += 1
    spos2count[curr_idx_msk - 1] = len(pos_col) - prev_idx
    assert np.sum(spos2count) == len(pos_col)
    return spos2count


def check_individuals_and_order(vcf_cols_, msp_cols_):
    return sorted(list(set(vcf_cols_).intersection(set(msp_cols_))))


def get_mask_data(vcf_ind_data, msp_mask_data, famid_id2pop_):
    vcf_cols = vcf_ind_data.columns
    msp_cols = msp_mask_data.columns
    new_cols = check_individuals_and_order(vcf_cols, msp_cols)
    if len(new_cols) < len(vcf_cols):
        print("WARNING: vcf file has fewer samples after intersecting with msp file!")
        vcf_ind = vcf_ind[new_cols]
    if len(new_cols) < len(msp_cols):
        print("WARNING: msp file has fewer samples after intersecting with vcf file!")
        msp_mask = msp_mask[new_cols]

    # Get a dictionary "group2colidx" with "population: [col_idx_of_sample]".
    group2colidx = defaultdict(list)
    for col_idx, sample in enumerate(new_cols):
        sample_id = sample[:-2]  # remove ".1" or ".0"
        if sample_id in famid_id2pop_:
            group2colidx[famid_id2pop_[sample_id]].append(col_idx)

    mask_snp_mtx = np.ma.MaskedArray(
        vcf_ind_data.to_numpy(dtype=int), mask=msp_mask_data.to_numpy(), fill_value=0
    )
    mask_mtx = 1 - msp_mask_data.to_numpy(dtype=int)  # mask = filtered
    aggr_frq = pd.DataFrame(
        {
            f"{pop:s}_MSK": np.sum(mask_snp_mtx[:, colidx_list], axis=1).filled(0)
            for pop, colidx_list in group2colidx.items()
        },
        dtype=np.float32,
    )
    aggr_tot = pd.DataFrame(
        {
            f"{pop:s}_MSK_CT": np.sum(mask_mtx[:, colidx_list], axis=1)
            for pop, colidx_list in group2colidx.items()
        },
        dtype=np.float32,
    )

    aggr_frq = pd.DataFrame(
        np.divide(
            aggr_frq.to_numpy(),
            aggr_tot.to_numpy(),
            out=np.zeros(aggr_frq.shape),
            where=(aggr_tot.to_numpy() != 0),
        ),
        dtype=np.float32,
        columns=aggr_frq.columns,
    )

    aggr_data = pd.concat([aggr_frq, aggr_tot], axis=1)
    return aggr_data


@ray.remote
def get_chrom_data(chr_idx, famid_id2pop_):
    vcf_data, vcf_ind_df, pos_col = get_vcf_data(chr_idx, famid_id2pop_)
    msp_mask_data, spos_col = get_msp_data(chr_idx, famid_id2pop_)

    spos2count = calculate_gnomix_row_repeat(pos_col, spos_col)
    msp_mask_df = pd.DataFrame(
        np.repeat(msp_mask_data[0], repeats=spos2count, axis=0),
        columns=msp_mask_data[1],
    )
    mask_data = get_mask_data(vcf_ind_df, msp_mask_df, famid_id2pop_)
    aggr_data = pd.concat([vcf_data, mask_data], axis=1)
    return aggr_data


@ray.remote
def write_data(group_, chrom_data_, save_file_path_):
    include_cols = [
        "#CHROM",
        "POS",
        "ALT",
        "REF",
        group_,
        f"{group_}_CT",
        f"{group_}_MSK",
        f"{group_}_MSK_CT",
    ]
    rename_cols = {
        "#CHROM": "CHR",
        "POS": "SNP",
        "ALT": "A1",
        "REF": "A2",
        group_: "MAF",
        f"{group_}_CT": "NCHROBS",
        f"{group_}_MSK": "MAF_MSK",
        f"{group_}_MSK_CT": "NCHROBS_MSK",
    }
    df_group = chrom_data_[include_cols]
    df_group = df_group.rename(columns=rename_cols)
    df_group.to_csv(
        f"{save_file_path_:s}/{group_:s}.freq", na_rep="NaN", index=False, sep="\t"
    )


if __name__ == "__main__":
    mask_group = args.mask_group

    # Create a dictionary "pop_info_dict" with "famid_id: pop"
    pop_info = pd.read_csv("hawaiiPopInfo.csv", header=0, usecols=["famid_id", "pop"])
    famid_id2pop = {data.values[0]: data.values[1] for _, data in pop_info.iterrows()}

    subpop_info = {}
    group_names = get_subpop_map(subpop_info, famid_id2pop)
    check_subpop_input(subpop_info, mask_group)

    # Concatenate data from all chromosomes 1-22.
    print("Start to read vcf and msp files from each chromosome!")
    chrom_data_remote = [get_chrom_data.remote(i, famid_id2pop) for i in range(1, 23)]
    chrom_data_list = ray.get(chrom_data_remote)
    chrom_data = pd.concat(chrom_data_list)
    print("Data aggregation finished!")

    # Save data to the given path.
    save_file_path = os.path.join(os.path.abspath("."), f"mask_data_{args.mask_group}")
    try:
        os.makedirs(save_file_path)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise

    print("Start to write files!")
    write_data_remote = [
        write_data.remote(group, chrom_data, save_file_path) for group in group_names
    ]
    ray.get(write_data_remote)
