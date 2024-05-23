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
import subprocess
import sys
from tqdm.auto import trange

import gene_stats as gst
import pop_aggr as paggr

# +-*=+-*=+-*=+-*=+-*=+-*=+-*=+-*=+*-=+-*=+*-=+-*=+-*=+-*=+-*=+-*= #
#                           Reading Inputs                         #
# +-*=+-*=+-*=+-*=+-*=+-*=+-*=+-*=+*-=+-*=+*-=+-*=+-*=+-*=+-*=+-*= #

argp = argparse.ArgumentParser()
argp.add_argument('stat', help="the statistics to compute (H stands for heterozygosity)",
                  choices=["F2", "F3", "F4", "Pi", "Psi", "F_ST", "H"])
argp.add_argument('-f', '--file', nargs='+', required=True,
                  help="(REQUIRED) files containing a list of population names for "\
                       "computation, each file has one column of names. For a "\
                       "binary statistics such as Psi, one input file computes "\
                       "self-match pairs, two input files computes cross-match pairs")
argp.add_argument('-t', '--data_dir', required=True,
                  help="(REQUIRED) frequency data file directory")
argp.add_argument('-o', '--output_dir', required=True,
                  help="(REQUIRED) output file directory")
argp.add_argument('-b', '--blocksize', type=int, default=500,
                  help="block size for block bootstrap")
argp.add_argument('-n', '--num_replicates', type=int, default=100,
                  help="number of replicates for block bootstrap")
argp.add_argument('-r', '--ref_group', nargs='+', metavar='REF_GROUP',
                  help='(F3/Psi only) reference population group for F3/Psi '\
                       'statistics, if multiple population group names, create '\
                       'an aggregated population')
argp.add_argument('-D', '--DAF', type=float, default=-1.0,
                  help='(Psi only) derived allele frequency threshold, '\
                       '(a float number between [0.0, 0.5], 0.5 represents '\
                       'no polarization)')
argp.add_argument('-d', '--downsample_size', type=int, default=2,
                  help='(Psi only) downsampling size for Psi statistics')
argp.add_argument('-m', '--mask', type=bool, default=False,
                  help='Use masked data for genetic statistics computation')
argp.add_argument('--rm_DA_files', type=bool, default=False,
                  help='(Boolean) whether we will remove existing files that '\
                       'indicates SNP positions for derived alleles at the new run')


class StatisticsParameters(object):
    def __init__(self, args_):
        # Absolute path to input data's directory
        self.data_dir = os.path.abspath(args_.data_dir)
        # Standardized name for parsing genetic statistics
        self.stats_name = args_.stat.lower().replace('_', '')
        # The number of input populations for computation
        self.pop_size = 0
        if self.stats_name == 'h':
            self.stats_name = 'heterozygosity'
            self.pop_size = 1
        elif self.stats_name == 'f3':
            self.pop_size = 3
        elif self.stats_name == 'f4':
            self.pop_size = 4
        else:
            self.pop_size = 2

        # Absolute path to output statistics files
        self.save_dir = os.path.join(args_.output_dir,
                                     f"{self.stats_name:s}_output")
        self.save_dir = os.path.abspath(self.save_dir)
        # Create a directory (if not exist) for the output of the statistics.
        try:
            os.makedirs(self.save_dir)
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise
        
        # Read lists of population names for computing statistics.
        num_pop_files = len(args_.file)
        if ((self.pop_size == 1) and num_pop_files != 1) or \
            ((self.pop_size == 2) and (num_pop_files not in [1, 2])) or \
            ((self.stats_name == 'f3') and (num_pop_files not in [1, 2])) or \
            ((self.stats_name == 'f4') and num_pop_files != 3):
            print("Please specify a correct number of population list files for "
                  "computing genetic statistics:\n"
                  "Unary input statistics St(A): 1 file,\n"
                  "Binary input statistics St(A, B): 1 file (self-match), "
                  "2 files (cross-match),\n"
                  "F3(O; A, B) (with reference population O specified): "
                  "1 file (self-match), 2 files (cross-match),\n"
                  "F4(A, B; C, D) (with reference population C specified): "
                  "3 files (cross-match for A, B, and D).\n")
            sys.exit(0)
        self.populi_list = []
        for populus_file in args_.file:
            self.populi_list.append(self.read_population_name(populus_file))

        self.aggr_file_abspath = ""
        self.DA_file_abspath = ""
        self.output_file_suffix = ""
        self.output_file_abspath = ""
        self.blocksize = args_.blocksize
        self.num_replicates = args_.num_replicates
        self.downsample_size = args_.downsample_size
        self.is_masked = args_.mask


    def read_population_name(self, populus_file):
        '''
        Read a list of population names in a given file.
        Args:
            populus_file_       A given file containing population names
        Returns:
            populus_list_       A list of population names
        '''
        populus_list = []
        with open(populus_file, 'r') as f:
            for populus in f:
                if populus.strip() != "":
                    populus_list.append(populus.strip())
        return populus_list


# +-*=+-*=+-*=+-*=+-*=+-*=+-*=+-*=+*-=+-*=+*-=+-*=+-*=+-*=+-*=+-*= #
#         Information from Reference Population (F3/F4/Psi)        #
# +-*=+-*=+-*=+-*=+-*=+-*=+-*=+-*=+*-=+-*=+*-=+-*=+-*=+-*=+-*=+-*= #

def aggregate_ref_pop(data_dir, stats_name, pop_data_files, aggr_filename,
                      DA_filename, thres):
    '''
    Aggregate specified reference population into one frequency file.
    Args:
        data_dir            Absolute path to input data's directory
        stats_name          Genetic statistics name
        pop_data_files      A list of names for reference populations
        aggr_filename       Name for aggregated population frequency file
        DA_filename         Name for derived allele SNP position file
        thres               A float for derived allele frequency threshold
    '''
    aggr_file_abspath = os.path.join(data_dir, aggr_filename)
    DA_file_abspath = os.path.join(data_dir, DA_filename)
    if (not os.path.exists(aggr_file_abspath)) or \
        (stats_name == 'psi' and not os.path.exists(DA_file_abspath)):
        print("Reference population groups detected without "
              "an aggregated population file or for Psi, a derived allele SNP "
              "position file.")
        print("Creating the aggregated population file (and for Psi, its "
              "corresponding derived allele SNP position file)...")
        paggr.construct_aggr_file(data_dir, stats_name, pop_data_files,
                                  aggr_file_abspath, DA_file_abspath, thres)
        print("Creating task completed!\n")
    return aggr_file_abspath, DA_file_abspath


def read_ref_pop_info(params, args_):
    '''
    Read information from reference populations for certain statistics (F3/F4/Psi).
    Interactive cell for user to specify the outgroup populations (O) if
    "stats_name" is "F3"/"F4", or to specify the reference populations
    (determining SNP positions for derived alleles) if "stats_name" is "psi".
    Args:
        params          The parameters for computing statistics
        args_           The input arguments loaded by argument parser
    '''
    stats_name = params.stats_name
    data_dir = params.data_dir
    output_file_suffix = ""
    if stats_name in ['f3', 'f4', 'psi']:
        try:
            # Specify reference populations.
            if args_.ref_group != None:
                outgroup_name = args_.ref_group
            else:
                # Specify customerized reference populations.
                print("Please specify the reference (outgroup) populations "
                        "for F3/F4/psi statistics (first-letter capitalized, "
                        "separated by space, if many, e.g., Atayal Bunun Paiwan)")
                print("(If there is more than one population, we will create "
                        "a new aggregated population based on the inputs, "
                        "named by the first three letters of each population)")
                outgroup_name = input().strip().split()
                if len(outgroup_name) == 0:
                    raise("Invalid input for outgroup populations.")

            # (Psi only) Specify derived allele frequency threshold.
            thres = -1.0
            if stats_name == 'psi':
                thres = args_.DAF
                while thres < -0.0 or thres > 0.5:
                    print("Please specify derived allele frequency threshold "
                          "(a float in [0.0, 0.5], where 0.5 indicates no "
                          "polarization): ")
                    thres = float(input().strip().split()[0])
            
            # Create a name format of the aggregated population file
            # (e.g., "AtaBunPai.freq").
            if len(outgroup_name) == 1:
                aggr_filename = outgroup_name[0] + ".freq"
            else:
                aggr_filename = "".join([pop.strip().split(".")[0][:3]\
                    .capitalize() for pop in outgroup_name]) + ".freq"
            
            # (Psi only) Create a name format of the derived allele SNP position
            # file (e.g., "psi_AtaBunPai_frq_05.freq").
            DA_filename = ""
            if stats_name == 'psi':
                DAF_str = str(thres).split(".")[-1][:2] # 0.05 -> "05"
                DAF_str += "0" * (2 - len(DAF_str))     # 0.0 -> "00"
                DA_filename = f"psi_{aggr_filename.split('.')[0]:s}_frq_"\
                              f"{DAF_str:s}.txt"
            
            # Aggregate specified reference population into one frequency file.
            pop_data_files = [pop + ".freq" for pop in outgroup_name]
            aggr_file_abspath, DA_file_abspath = aggregate_ref_pop(data_dir,
                                                                   stats_name,
                                                                   pop_data_files,
                                                                   aggr_filename,
                                                                   DA_filename,
                                                                   thres)
            params.aggr_file_abspath = aggr_file_abspath
            params.DA_file_abspath = DA_file_abspath

            # (Psi only) Form output file suffix for statistics.
            if stats_name == 'psi':
                if abs(thres - 0.5) <= 1e-5:
                    # No polarization (usually a window-based SNP using Gnomix
                    # output files)
                    output_file_suffix = f"_None_"
                else:
                    output_file_suffix = f"_{aggr_filename.split('.')[0]:s}_"
                output_file_suffix += f"frq_{DAF_str:s}_{args_.downsample_size:d}"
            else:
                output_file_suffix = f"_{aggr_filename.split('.')[0]:s}"
            
        except Exception as e:
            raise
    
    # Standardize output file name.
    params.output_file_suffix = output_file_suffix
    params.output_file_abspath = os.path.join(params.save_dir,
                                              f"{stats_name:s}_stat"
                                              f"{output_file_suffix:s}.txt")
    print(params.output_file_abspath)


# +-*=+-*=+-*=+-*=+-*=+-*=+-*=+-*=+*-=+-*=+*-=+-*=+-*=+-*=+-*=+-*= #
#                  Computation of the Statistics                   #
# +-*=+-*=+-*=+-*=+-*=+-*=+-*=+-*=+*-=+-*=+*-=+-*=+-*=+-*=+-*=+-*= #

def stats_input(pops, params):
    '''
    Create a new object of the given genetic statistics based on the given
    configuration and compute the mean and standard error of the statistics.
    Args:
        pops                A list of population names
        params              A struct that contains input parameters of statistics
    Return:
        stats_mean          The mean of the statistics
        stats_std_error     The standard error of the statistics
    ''' 
    pops = [pop + ".freq" for pop in pops]          # A list of population names
    stats_name = params.stats_name                  # Genetic statistics
    block_size = params.blocksize                   # Block size for block bootstrap
    num_replicates = params.num_replicates          # Number of replicates for the bootstrap
    data_dir = params.data_dir                      # Directory for frequency files (.freq)
    output_filename = params.output_file_abspath    # Name for the output text file (.txt) of the statistics
    outgroup_filename = params.aggr_file_abspath    # (F3/Psi) Frequency file name for the outgroup populations
    DA_filename = params.DA_file_abspath            # (Psi) Derived allele position file name
    downsample_size = params.downsample_size        # (Psi) Downsampling size (int)
    is_masked = params.is_masked                    # Whether to use masked frequency data

    if stats_name == 'heterozygosity':
        config = gst.UnaryInputConfig(data_dir=data_dir,
                                      pop1_file=pops[0],
                                      block_size=block_size,
                                      num_replicates=num_replicates,
                                      output_file=output_filename,
                                      is_masked=is_masked)
        stats_obj = gst.Heterozygosity(config)
    elif stats_name in ['f2', 'fst', 'pi']:
        config = gst.BinaryInputConfig(data_dir=data_dir,
                                       pop1_file=pops[0],
                                       pop2_file=pops[1],
                                       block_size=block_size,
                                       num_replicates=num_replicates,
                                       output_file=output_filename,
                                       is_masked=is_masked)
        if stats_name == 'f2':
            stats_obj = gst.F2(config)
        elif stats_name == 'fst':
            stats_obj = gst.Fst(config)
        else:
            stats_obj = gst.Pi(config)
    elif stats_name == 'f3':
        config = gst.TernaryInputConfig(data_dir=data_dir,
                                        pop1_file=pops[0],
                                        pop2_file=pops[1],
                                        block_size=block_size,
                                        num_replicates=num_replicates,
                                        output_file=output_filename,
                                        pop3_file=outgroup_filename,
                                        is_masked=is_masked)
        stats_obj = gst.F3(config)
    elif stats_name == 'f4':
        config = gst.QuaternaryInputConfig(data_dir=data_dir,
                                           pop1_file=pops[0],
                                           pop2_file=pops[1],
                                           block_size=block_size,
                                           num_replicates=num_replicates,
                                           output_file=output_filename,
                                           pop3_file=outgroup_filename,
                                           pop4_file=pops[2],
                                           is_masked=is_masked)
        stats_obj = gst.F4(config)
    else:
        config = gst.PsiConfig(data_dir=data_dir,
                               pop1_file=pops[0],
                               pop2_file=pops[1],
                               block_size=block_size,
                               num_replicates=num_replicates,
                               output_file=output_filename,
                               da_file=DA_filename,
                               downsample_size=downsample_size,
                               is_masked=is_masked)
        stats_obj = gst.Psi(config)
    
    stats_mean = stats_obj.get_mean()
    stats_std_error = stats_obj.get_std_error()
    return stats_mean, stats_std_error


# Compute the given statistics with all combinations of two population groups.
class UnaryComputation(object):
    def __init__(self, params):
        num_pop_files = len(params.populi_list)
        if num_pop_files != 1:
            print("The unary input genetic statistics requires 1 population "
                  "list file.")
            sys.exit(0)
        self.self_match_computation(params)

    def self_match_computation(self, params):
        populi_list = [np.array(pop_list) for pop_list in params.populi_list]
        populus_list1 = populi_list[0]
        num_pop1 = len(populus_list1)
        out_mtx = np.zeros(num_pop1, )
        out_mtx_se = np.zeros(num_pop1, )
        for i in trange(num_pop1):
            pop1 = populus_list1[i]
            stats_mean, stats_std_error = stats_input([pop1], params)
            out_mtx[i], out_mtx_se[i] = stats_mean, stats_std_error
            print(f"{pop1:s}: {stats_mean:f}, {stats_std_error:f}")
        self.save_output_file(params, out_mtx, out_mtx_se, populus_list1)
    
    def save_output_file(self, params, out_mtx, out_mtx_se, populus_list1):
        stat_mtx_file = f"{params.save_dir:s}/{params.stats_name:s}_mtx"\
                        f"{params.output_file_suffix:s}"
        print(f"Output statistics matrix file: {stat_mtx_file:s}.npz")
        np.savez(stat_mtx_file, stat_mtx=out_mtx, stat_mtx_se=out_mtx_se,
                 pop1_list=populus_list1)


class GeneralComputation(object):
    def __init__(self, params):
        num_pop_files = len(params.populi_list)
        if num_pop_files not in [1, 2]:
            print("The input genetic statistics requires 1 or 2 population "
                  "list files.")
            sys.exit(0)
        self.validate_inputs(params)

        if num_pop_files == 1:
            self.self_match_computation(params)
        else:
            self.cross_match_computation(params)

    def validate_inputs(self, params):
        """
        Confirm that all .freq files contain the exact same SNPs and alternate alleles (A1).
        """
        pop_dir = params.data_dir
        freq_files = set()
        # Get all population frequency files
        for pop_file in params.populi_list:
            for pop in pop_file:
                freq_files.add(os.path.join(pop_dir, pop + ".freq"))
        freq_files = list(freq_files)

        ref = pd.read_csv(freq_files[0], usecols=["CHR", "SNP", "A1"], sep='\t', header=0)
        for i in range(1, len(freq_files)):
            curr_pop = pd.read_csv(freq_files[i], usecols=["CHR", "SNP", "A1"], sep='\t', header=0)
            if not ref.equals(curr_pop):
                print("All population frequency files must contain the exact same SNPs and minor alleles. "
                      f"Differing row(s) found in {freq_files[i].split('/')[-1]}.")
                sys.exit(0)

    def self_match_computation(self, params):
        # self join: Group A - Group A
        # Compute the given statistics for (n choose 2) pairs (x_i, x_j) in one
        # population list.
        populi_list = [np.array(pop_list) for pop_list in params.populi_list]
        populus_list1 = populi_list[0]
        num_pop1 = len(populus_list1)
        out_mtx = np.zeros((num_pop1, num_pop1))
        out_mtx_se = np.zeros((num_pop1, num_pop1))
        for i in trange(num_pop1):
            pop1 = populus_list1[i]
            for j in range(i + 1, num_pop1):
                pop2 = populus_list1[j]
                stats_mean, stats_std_error = stats_input([pop1, pop2], params)
                out_mtx[i][j], out_mtx_se[i][j] = stats_mean, stats_std_error
                print(f"{pop1:s}-{pop2:s}: "\
                      f"{out_mtx[i][j]:f}, {out_mtx_se[i][j]:f}")
        if params.stats_name == 'psi':
            # skew-symmetric matrix for Psi (as psi_{a, b} = -psi_{b, a})
            out_mtx += -out_mtx.T
        else:
            out_mtx += out_mtx.T
        out_mtx_se += out_mtx_se.T
        self.save_output_file(params, out_mtx, out_mtx_se, populus_list1,
                              populus_list1)

    def cross_match_computation(self, params):
        # inner join: Group A - Group B
        # Compute the given statistics for all pairs (x_i, y_j) in two
        # population lists.
        populi_list = [np.array(pop_list) for pop_list in params.populi_list]
        populus_list1, populus_list2 = populi_list
        num_pop1, num_pop2 = len(populus_list1), len(populus_list2)
        out_mtx = np.zeros((num_pop1, num_pop2))
        out_mtx_se = np.zeros((num_pop1, num_pop2))
        for i in trange(num_pop1):
            pop1 = populus_list1[i]
            for j in range(num_pop2):
                pop2 = populus_list2[j]
                stats_mean, stats_std_error = stats_input([pop1, pop2], params)
                out_mtx[i][j], out_mtx_se[i][j] = stats_mean, stats_std_error
                print(f"{pop1:s}-{pop2:s}: "\
                      f"{out_mtx[i][j]:f}, {out_mtx_se[i][j]:f}")
        self.save_output_file(params, out_mtx, out_mtx_se, populus_list1,
                              populus_list2)
    
    def save_output_file(self, params, out_mtx, out_mtx_se, populus_list1,
                         populus_list2):
        stat_mtx_file = f"{params.save_dir:s}/{params.stats_name:s}_mtx"\
                        f"{params.output_file_suffix:s}"
        print(f"Output statistics matrix file: {stat_mtx_file:s}.npz")
        np.savez(stat_mtx_file, stat_mtx=out_mtx, stat_mtx_se=out_mtx_se,
                 pop1_list=populus_list1, pop2_list=populus_list2)


class F4Computation(object):
    def __init__(self, params):
        num_pop_files = len(params.populi_list)
        if num_pop_files != 3:
            print("The genetic statistics F4 requires 3 population "
                  "list files representing A, B, and D in F4(A, B; C, D).")
            sys.exit(0)
        self.validate_inputs(params)
        self.cross_match_computation(params)

    def validate_inputs(self, params):
        """
        Confirm that all .freq files contain the exact same SNPs and alternate alleles (A1).
        """
        pop_dir = params.data_dir
        freq_files = set()
        # Get all population frequency files
        for pop_file in params.populi_list:
            for pop in pop_file:
                freq_files.add(os.path.join(pop_dir, pop + ".freq"))
        freq_files = list(freq_files)

        ref = pd.read_csv(freq_files[0], usecols=["CHR", "SNP", "A1"], sep='\t', header=0)
        for i in range(1, len(freq_files)):
            curr_pop = pd.read_csv(freq_files[i], usecols=["CHR", "SNP", "A1"], sep='\t', header=0)
            if not ref.equals(curr_pop):
                print("All population frequency files must contain the exact same SNPs and minor alleles. "
                      f"Differing row(s) found in {freq_files[i].split('/')[-1]}.")
                sys.exit(0)

    def cross_match_computation(self, params):
        # inner join: Group A - Group B - Group D
        # Compute the given statistics for all pairs (x_i, y_j, z_k) in three
        # population lists.
        populi_list = [np.array(pop_list) for pop_list in params.populi_list]
        populus_list1, populus_list2, populus_list4 = populi_list
        num_pop1, num_pop2, num_pop4 = len(populus_list1), len(populus_list2), \
                                       len(populus_list4)
        out_mtx = np.zeros((num_pop1, num_pop2, num_pop4))
        out_mtx_se = np.zeros((num_pop1, num_pop2, num_pop4))
        for k in trange(num_pop4):
            pop4 = populus_list4[k]
            for i in range(num_pop1):
                pop1 = populus_list1[i]
                for j in range(num_pop2):
                    pop2 = populus_list2[j]
                    stats_mean, stats_std_error = stats_input([pop1, pop2, pop4],
                                                              params)
                    out_mtx[i][j][k], out_mtx_se[i][j][k] = stats_mean, stats_std_error
                    print(f"F4({pop1:s}, {pop2:s}; O, {pop4:s}): "\
                          f"{out_mtx[i][j][k]:f}, {out_mtx_se[i][j][k]:f}")
        self.save_output_file(params, out_mtx, out_mtx_se, populus_list1,
                              populus_list2, populus_list4)
    
    def save_output_file(self, params, out_mtx, out_mtx_se, populus_list1,
                         populus_list2, populus_list4):
        stat_mtx_file = f"{params.save_dir:s}/{params.stats_name:s}_mtx"\
                        f"{params.output_file_suffix:s}"
        print(f"Output statistics matrix file: {stat_mtx_file:s}.npz")
        np.savez(stat_mtx_file, stat_mtx=out_mtx, stat_mtx_se=out_mtx_se,
                 pop1_list=populus_list1, pop2_list=populus_list2,
                 pop4_list=populus_list4)


if __name__ == "__main__":
    args = argp.parse_args()
    stat_params = StatisticsParameters(args)

    read_ref_pop_info(stat_params, args)

    subprocess.check_output([f"touch {stat_params.output_file_abspath:s}"],
                            shell=True)
    if stat_params.stats_name == "heterozygosity":
        compute = UnaryComputation(stat_params)
    elif stat_params.stats_name == "f4":
        compute = F4Computation(stat_params)
    else:
        compute = GeneralComputation(stat_params)

    print("Computation Succeeded!")

    # (Psi only) Remove any intermediate files generated during the computation.
    if args.rm_DA_files and stat_params.stats_name == 'psi':
        try:
            file_list = glob.glob(f'{stat_params.data_dir}/psi_*.txt')
            for file in file_list:
                os.remove(file)
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise
