# Compute F3/Pi/Psi/F_ST/Heterozygosity statistics of any two populations from a given
# population list file. The program accept three or four parameters, in which
# the last two parameters indicate the target population list files.
# 1. If we provide one population list file, the program will compute the given
# statistics for (n choose 2) pairs (x_i, x_j) in that population list.
# 2. If we provide two population list files, the program will compute the given
# statistics for all pairs (x_i, y_j) in the two population lists.

import os
import glob
import errno
import numpy as np
import pandas as pd
import argparse
from tqdm.auto import trange

import gene_stats as gsr

# +-*=+-*=+-*=+-*=+-*=+-*=+-*=+-*=+*-=+-*=+*-=+-*=+-*=+-*=+-*=+-*= #
#                           Reading Inputs                         #
# +-*=+-*=+-*=+-*=+-*=+-*=+-*=+-*=+*-=+-*=+*-=+-*=+-*=+-*=+-*=+-*= #

argp = argparse.ArgumentParser()
argp.add_argument('stat', help="the statistics to compute (H stands for heterozygosity)",
                  choices=["F3", "Pi", "Psi", "F_ST", "H"])
argp.add_argument('-f', '--file', nargs='+', required=True,
                  help="(REQUIRED) files containing a list of population names for computation, "\
                       "each file has one column of names, one file for self-match pairs, "\
                       "two files for cross-match pairs")
argp.add_argument('-t', '--data_dir', required=True,
                  help="(REQUIRED) frequency data file directory")
argp.add_argument('-o', '--output_dir', required=True,
                  help="(REQUIRED) output file directory")
argp.add_argument('-b', '--blocksize', type=int, default=500,
                  help="block size for block bootstrap")
argp.add_argument('-n', '--num_replicates', type=int, default=100,
                  help="number of replicates for block bootstrap")
argp.add_argument('-r', '--ref_group', nargs='+', metavar='GROUP',
                  help='(F3/Psi only) reference population group for F3/Psi statistics, '\
                       'if multiple population group names, create an aggregated population')
argp.add_argument('-D', '--DAF', type=float, default=-1.0,
                  help='(Psi only) derived allele frequency threshold, '\
                       '(a float number between [0.0, 0.5])')
argp.add_argument('-d', '--downsample_size', type=int, default=2,
                  help='(Psi only) downsampling size for Psi statistics')
argp.add_argument('--rm_DA_files', type=bool, default=False,
                  help='(Boolean) whether we will remove existing files that indicates SNP '\
                       'positions for derived alleles at the new run')
args = argp.parse_args()

# Define file directory for data.
data_dir = args.data_dir
output_dir = args.output_dir

# Define a saving file path.
stats_name = args.stat.lower().replace('_', '')
if stats_name == 'h':
    stats_name = 'heterozygosity'
    pop_size = 1
elif stats_name == 'f3':
    pop_size = 3
else:
    pop_size = 2
save_file_path = os.path.join(output_dir, f"{stats_name:s}_output")

# Create a directory (if not exist) for the output of the statistics.
try:
    os.makedirs(save_file_path)
except OSError as e:
    if e.errno != errno.EEXIST:
        raise


def read_population_name(populus_file_):
    '''
    Read a list of population names in a given file.
    Args:
        populus_file_       A given file containing population names
    '''
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


# +-*=+-*=+-*=+-*=+-*=+-*=+-*=+-*=+*-=+-*=+*-=+-*=+-*=+-*=+-*=+-*= #
#          Information from Reference Population (F3/Psi)          #
# +-*=+-*=+-*=+-*=+-*=+-*=+-*=+-*=+*-=+-*=+*-=+-*=+-*=+-*=+-*=+-*= #

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

        thres = args.DAF
        if stats_name == 'psi' and thres == -1.0:
            print("Please specify the derived allele frequency threshold "
                  "(a float number in [0.0, 0.5]): ")
            thres = float(input().strip().split()[0])
        
        # Create a name format of the aggregated population file (e.g., "AtaBunPai.freq").
        if len(outgroup_name) == 1:
            aggr_filename = outgroup_name[0] + ".freq"
        else:
            aggr_filename = "".join([pop.strip().split(".")[0][:3].capitalize()
                                     for pop in outgroup_name]) + ".freq"
        
        # Create a name format of the derived allele SNP position file
        # (e.g., "psi_atbupa_frq_05.freq").
        DAF_str = str(thres).split(".")[-1][:2] # 0.05 -> "05"
        DAF_str += "0" * (2 - len(DAF_str)) # 0.0 -> "00"
        DA_filename = f"psi_{aggr_filename.split('.')[0]:s}_frq_{DAF_str:s}.txt"

        pop_data_files = [pop + ".freq" for pop in outgroup_name]
        if (not os.path.exists(os.path.join(data_dir, aggr_filename))) or \
           (stats_name == 'psi' and not os.path.exists(os.path.join(data_dir, DA_filename))):
            print("Reference population groups detected without "
                  "an aggregated population file or a DAF-SNP position file.")
            print("Creating aggregated population file...")
            os.system(f"python3 pop_aggr.py {data_dir:s} {' '.join(pop_data_files):s} {thres:f} "
                      f"{aggr_filename:s} {DA_filename:s}")
            print("Creating task completed!")
            print()
    except Exception as e:
        raise

# +-*=+-*=+-*=+-*=+-*=+-*=+-*=+-*=+*-=+-*=+*-=+-*=+-*=+-*=+-*=+-*= #
#                  Computation of the Statistics                   #
# +-*=+-*=+-*=+-*=+-*=+-*=+-*=+-*=+*-=+-*=+*-=+-*=+-*=+-*=+-*=+-*= #

def stats_input(pops, stats_name, block_size, num_replicates, data_dir,
                output_filename, outgroup_filename="", DA_filename="",
                downsample_size=args.downsample_size):
    '''
    Create a new object of the given genetic statistics based on the given
    configuration and compute the mean and standard error of the statistics.
    Args:
        pops                A list of population names
        stats_name          Genetic statistics
        block_size          Block size for block bootstrap
        num_replicates      Number of replicates for the bootstrap
        data_dir            Directory for frequency files (.freq)
        output_filename     Name for the output text file (.txt) of the statistics
        outgroup_filename   (F3/Psi) Frequency file name for the outgroup populations
        DA_filename         (Psi) Derived allele position file name
        downsample_size     (Psi) Downsampling size (int)
    Return:
        stats_mean          The mean of the statistics
        stats_std_error     The standard error of the statistics
    ''' 
    pops = [pop + ".freq" for pop in pops]
    if stats_name == 'heterozygosity':
        config = gsr.UnaryInputConfig(data_dir=data_dir, pop1_file=pops[0],
                                      block_size=block_size,
                                      num_replicates=num_replicates,
                                      output_file=output_filename)
        stats_obj = gsr.Heterozygosity(config)
    elif stats_name in ['fst', 'pi']:
        config = gsr.BinaryInputConfig(data_dir=data_dir, pop1_file=pops[0],
                                       pop2_file=pops[1], block_size=block_size,
                                       num_replicates=num_replicates,
                                       output_file=output_filename)
        if stats_name == 'fst':
            stats_obj = gsr.Fst(config)
        else:
            stats_obj = gsr.Pi(config)
    elif stats_name == 'f3':
        config = gsr.TrinaryInputConfig(data_dir=data_dir, pop1_file=pops[0],
                                        pop2_file=pops[1],
                                        block_size=block_size,
                                        num_replicates=num_replicates,
                                        output_file=output_filename,
                                        pop3_file=outgroup_filename)
        stats_obj = gsr.F3(config)
    else:
        config = gsr.PsiConfig(data_dir=data_dir, pop1_file=pops[0],
                               pop2_file=pops[1], block_size=block_size,
                               num_replicates=num_replicates,
                               output_file=output_filename,
                               da_file=DA_filename,
                               downsample_size=downsample_size)
        stats_obj = gsr.Psi(config)
    
    stats_mean = stats_obj.get_mean()
    stats_std_error = stats_obj.get_std_error()
    return stats_mean, stats_std_error


# Standardize output file name.
if stats_name == 'f3':
    output_file_suffix = f"_{aggr_filename.split('.')[0]:s}"
elif stats_name == 'psi':
    output_file_suffix = f"_{aggr_filename.split('.')[0]:s}_frq_{DAF_str:s}_{args.downsample_size:d}"
else:
    output_file_suffix = ""
output_filename = os.path.join(save_file_path,
                               f"{stats_name:s}_stat{output_file_suffix:s}.txt")

# Compute the given statistics with all combinations of two population groups.
os.system(f"touch {output_filename:s}")
if pop_size == 1: # e.g., heterozygosity
    assert len(args.file) == 1
    # Compute the given statistics for (x_i) in one population list.
    for i in trange(len(populus_list1)):
        populus1 = populus_list1[i]
        stats_mean, stats_std_error = stats_input([populus1], stats_name, args.blocksize,
                                                  args.num_replicates, data_dir, output_filename,
                                                  outgroup_filename=aggr_filename,
                                                  DA_filename=DA_filename,
                                                  downsample_size=args.downsample_size)
        print(f"{populus1:s}: {stats_mean:f}, {stats_std_error:f}")
else:
    if len(args.file) == 1:
        # self join: Group A - Group A
        # Compute the given statistics for (n choose 2) pairs (x_i, x_j) in one
        # population list.
        #
        # Create an output matrix for the given statistics.
        out_mtx = np.zeros((len(populus_list1), len(populus_list1)))
        out_mtx_se = np.zeros((len(populus_list1), len(populus_list1)))
        for i in trange(len(populus_list1)):
            populus1 = populus_list1[i]
            for j in range(i + 1, len(populus_list1)):
                populus2 = populus_list1[j]
                stats_mean, stats_std_error = stats_input([populus1, populus2], stats_name,
                                                          args.blocksize, args.num_replicates,
                                                          data_dir, output_filename,
                                                          outgroup_filename=aggr_filename,
                                                          DA_filename=DA_filename,
                                                          downsample_size=args.downsample_size)
                out_mtx[i][j], out_mtx_se[i][j] = stats_mean, stats_std_error
                print(f"{populus1:s}-{populus2:s}: {out_mtx[i][j]:f}, {out_mtx_se[i][j]:f}")
        if stats_name == 'psi':
            # skew-symmetric matrix for Psi (as psi_{a, b} = -psi_{b, a})
            out_mtx += -out_mtx.T
        else:
            out_mtx += out_mtx.T
        out_mtx_se += out_mtx_se.T

        # Save the output statistic matrix into the given "save_file_path".
        df = pd.DataFrame(out_mtx, index = populus_list1, columns = populus_list1)
        df_se = pd.DataFrame(out_mtx_se, index = populus_list1, columns = populus_list1)
    else:
        # inner join: Group A - Group B
        # Compute the given statistics for all pairs (x_i, y_j) in two
        # population lists.
        #
        # Create an output matrix for the given statistics.
        out_mtx = np.zeros((len(populus_list1), len(populus_list2)))
        out_mtx_se = np.zeros((len(populus_list1), len(populus_list2)))
        for i in trange(len(populus_list1)):
            populus1 = populus_list1[i]
            for j in range(len(populus_list2)):
                populus2 = populus_list2[j]
                stats_mean, stats_std_error = stats_input([populus1, populus2], stats_name,
                                                          args.blocksize, args.num_replicates,
                                                          data_dir, output_filename,
                                                          outgroup_filename=aggr_filename,
                                                          DA_filename=DA_filename,
                                                          downsample_size=args.downsample_size)
                out_mtx[i][j], out_mtx_se[i][j] = stats_mean, stats_std_error
                print(f"{populus1:s}-{populus2:s}: {out_mtx[i][j]:f}, {out_mtx_se[i][j]:f}")

        # Save the output statistic matrix into the given "save_file_path".
        df = pd.DataFrame(out_mtx, index = populus_list1, columns = populus_list2)
        df_se = pd.DataFrame(out_mtx_se, index = populus_list1, columns = populus_list2)

    df.to_csv(f"{save_file_path:s}/{stats_name:s}_mtx{output_file_suffix:s}.csv")
    print(f"{save_file_path:s}/{stats_name:s}_mtx{output_file_suffix:s}.csv")
    df_se.to_csv(f"{save_file_path:s}/{stats_name:s}_mtx_se{output_file_suffix:s}.csv")

# (Psi) Remove any intermediate files generated during the computation.
if args.rm_DA_files and stats_name == 'psi':
    try:
        file_list = glob.glob(f'{data_dir}/psi_*.txt')
        for file in file_list:
            os.remove(file)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise
