# Compute F3/pi/psi/F_ST statistics of any two populations from a given
# population list file. The program accept three or four parameters, in which
# the last two parameters indicate the target population list files.
# 1. If we provide one population list file, the program will compute the given
# statistics for (n choose 2) pairs (x_i, x_j) in that population list.
# 2. If we provide two population list files, the program will compute the given
# statistics for all pairs (x_i, y_j) in the two population lists.

import os
import sys
import errno
import numpy as np
import pandas as pd
import subprocess as sp
import argparse
from tqdm.auto import trange

argp = argparse.ArgumentParser()
argp.add_argument('stat', help="the statistics to compute",
                  choices=["F3", "Pi", "Psi", "F_ST"])
argp.add_argument('-f', '--file', nargs='+', required=True,
                  help="(REQUIRED) files containing a list of population names for computation, "\
                       "one file for self-match pairs, two files for cross-match pairs")
argp.add_argument('-b', '--blocksize', type=int, default=500,
                  help="block size for block bootstrap")
argp.add_argument('-r', '--ref_group', nargs='+', metavar='GROUP',
                  help='(F3/Psi only) reference population group for F3/Psi statistics, '\
                       'if multiple population group names, create an aggregated population')
argp.add_argument('-D', '--DAF', type=float, default=-1.0,
                  help='(Psi only) derived allele frequency threshold, '\
                       '(a float number between [0.0, 0.5])')
argp.add_argument('-d', '--downsample_size', type=int, default=2,
                  help='(Psi only) downsampling size for Psi statistics')
args = argp.parse_args()


# Define a saving file path.
stats_name = args.stat.lower().replace('_', '')
save_file_path = os.path.join(sys.path[0], f"{stats_name:s}_output")
output_file_name = f"{stats_name:s}_stat.txt"
# Create a directory (if not exist) for the output of the statistics.
try:
    os.makedirs(save_file_path)
except OSError as e:
    if e.errno != errno.EEXIST:
        raise
# TODO: Define file directory for data. (We will change this when adding npz2freq.py)
data_dir = "data_migrate"

# Read a list of population names.
def read_population_name(populus_file_):
    # Return a list of population names.
    # Args:
    #   populus_file_    A given file containing population names.
    populus_list_ = []
    with open(populus_file_, 'r') as f:
        for populus in f:
            if populus.strip() != "":
                populus_list_.append(populus.strip())
    return populus_list_

populus_file1 = args.file[0]
populus_list1 = read_population_name(populus_file1)
if len(args.file) > 1:
    # Take the second population name list.
    populus_file2 = args.file[1]
    populus_list2 = read_population_name(populus_file2)

# Interactive cell for user to specify the outgroup populations (O) if "stats_name"
# is "F3" or to specify the reference populations (determining SNP positions
# for derived alleles) if "stats_name" is "psi".
aggr_filename = ""
DA_filename = ""
if stats_name in ['f3', 'psi']:
    try:
        if args.ref_group != None:
            outgroup_name = args.ref_group
        else:
            print("Please specify the reference (outgroup) populations for F3/psi "
                "statistics (first-letter capitalized, separated by space, if many, "
                "e.g., Atayal Bunun Paiwan)")
            print("(If there is more than one population, we will create a new aggregated "
                "population based on the inputs, named by the first two letters of "
                "each population)")
            outgroup_name = input().strip().split()
            if len(outgroup_name) == 0:
                raise("Invalid input for outgroup populations.")
        outgroup_name = ["_".join([part.lower().capitalize() for part in pop.split('_')])
                         for pop in outgroup_name]

        thres = args.DAF
        if stats_name == 'psi' and thres == -1.0:
            print("Please specify the derived allele frequency threshold "
                  "(a float number in [0.0, 0.5]): ")
            thres = float(input().strip().split()[0])
        
        # Create a name format of the aggregated population file (e.g., "AtBuPa.freq").
        if len(outgroup_name) == 1:
            aggr_filename = outgroup_name[0] + ".freq"
        else:
            aggr_filename = "".join([pop.strip().split(".")[0][:2].capitalize()
                                 for pop in outgroup_name]) + ".freq"
        
        # Create a name format of the derived allele SNP position file
        # (e.g., "psi_atbupa_frq_05.freq").
        DAF_str = str(thres).split(".")[-1][:2] # 0.05 -> "05"
        DAF_str += "0" * (2 - len(DAF_str)) # 0.0 -> "00"
        DA_filename = f"psi_{aggr_filename.lower().split('.')[0]}_frq_{DAF_str}.txt"

        pop_data_files = [pop + ".freq" for pop in outgroup_name]
        if (not os.path.exists(os.path.join(data_dir, aggr_filename))) or \
           (stats_name == 'psi' and not os.path.exists(os.path.join(data_dir, DA_filename))):
            print("Multiple reference population groups detected without "
                  "an aggregated population file or a DAF-SNP position file.")
            print("Creating aggregated population file...")
            os.system(f"python3 pop_aggr.py {data_dir:s} {' '.join(pop_data_files):s} {thres:f} "
                      f"{aggr_filename:s} {DA_filename:s}")
            print("Creating task completed!")
            print()
    except Exception as e:
        raise

def rscript_input(pop1, pop2, stats_name, block_size,
                  data_dir, save_file_path, output_file_name,
                  outgroup_filename="", DA_filename="",
                  downsample_size=args.downsample_size):
    # Create a command line passed into R.
    basic_format = f"Rscript {stats_name:s}_ts.R {pop1:s}.freq {pop2:s}.freq {block_size:d} "\
                   f"{data_dir:s} {os.path.join(save_file_path, output_file_name):s} "
    if stats_name == 'f3':
        basic_format += f"{outgroup_filename:s}"
    if stats_name == 'psi':
        basic_format += f"{DA_filename:s} {downsample_size:d}"
    return basic_format.strip()


# Compute the given statistics with all combinations of two population groups.
os.system(f"touch {os.path.join(save_file_path, output_file_name):s}")
if len(args.file) == 1:
    # Create an output matrix for the given statistics.
    out_mtx = np.zeros((len(populus_list1), len(populus_list1)))

    # Compute the given statistics for (n choose 2) pairs (x_i, x_j) in one
    # population list.
    for i in trange(len(populus_list1)):
        populus1 = populus_list1[i]
        for j in range(i + 1, len(populus_list1)):
            populus2 = populus_list1[j]
            command_line = rscript_input(populus1, populus2, stats_name, args.blocksize,
                                         data_dir, save_file_path, output_file_name,
                                         outgroup_filename=aggr_filename, DA_filename=DA_filename,
                                         downsample_size=args.downsample_size)
            x = sp.check_output(command_line.split())
            out_mtx[i][j] = float(x.decode("utf-8").strip().split('\n')[-1])
            print(f"{populus1:s}-{populus2:s} {out_mtx[i][j]}")
    if stats_name == 'psi':
        # skew-symmetric matrix for Psi (as psi_{a, b} = -psi_{b, a})
        out_mtx += -out_mtx.T
    else:
        out_mtx += out_mtx.T
else:
    # Create an output matrix for the given statistics.
    out_mtx = np.zeros((len(populus_list1), len(populus_list2)))

    # Compute the given statistics for all pairs (x_i, y_j) in two
    # population lists.
    for i in trange(len(populus_list1)):
        populus1 = populus_list1[i]
        for j in range(len(populus_list2)):
            populus2 = populus_list2[j]
            command_line = rscript_input(populus1, populus2, stats_name, args.blocksize,
                                         data_dir, save_file_path, output_file_name,
                                         outgroup_filename=aggr_filename, DA_filename=DA_filename,
                                         downsample_size=args.downsample_size)
            x = sp.check_output(command_line.split())
            out_mtx[i][j] = float(x.decode("utf-8").strip().split('\n')[-1])
            print(f"{populus1:s}-{populus2:s} {out_mtx[i][j]}")

# Save the output statistic matrix into the given "save_file_path".
if len(args.file) > 1:
    df = pd.DataFrame(out_mtx, index = populus_list1, columns = populus_list2)
else:
    df = pd.DataFrame(out_mtx, index = populus_list1, columns = populus_list1)
df.to_csv(f"{save_file_path:s}/{stats_name:s}_mtx.txt")
