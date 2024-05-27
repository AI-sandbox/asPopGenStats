from abc import ABC, abstractmethod
from arch.bootstrap import CircularBlockBootstrap
from dataclasses import dataclass
import numpy as np
import os
import pandas as pd
import ray
from scipy.stats import hypergeom
import time

NUM_SNP_RATIO = 5  # Set the least number of SNP threshold by a ratio to block_size

# +-*=+-*=+-*=+-*=+-*=+-*=+-*=+-*=+*-=+-*=+*-=+-*=+-*=+-*=+-*=+-*= #
#                       Statistics Functions                       #
# +-*=+-*=+-*=+-*=+-*=+-*=+-*=+-*=+*-=+-*=+*-=+-*=+-*=+-*=+-*=+-*= #


@ray.remote
def compute_statistics_helper(df, block_size, stat, **kwargs):
    """
    Helper function for constructing one block bootstrap and executing the
    corresponding statistics function.
    Args:
        df              Dataframe for SNP loci x population (frequency & count) table
        block_size      Block size for block bootstrap
        stat            Genetic statistics
    Return:
        _               Computation output of the genetic statistics
    """
    boot = CircularBlockBootstrap(block_size, df)
    if stat == "heterozygosity":
        return boot.apply(heterozygosity, 1)
    elif stat == "f2":
        return boot.apply(f2, 1)
    elif stat == "fst":
        return boot.apply(fst, 1)
    elif stat == "f3":
        return boot.apply(f3, 1)
    elif stat == "f4":
        return boot.apply(f4, 1)
    elif stat == "pi":
        return boot.apply(pi, 1)
    elif stat == "psi":
        if kwargs.get("downsample_size") == None:
            raise ValueError("Please provide downsampling size for Psi computation.")
        else:
            return boot.apply(psi, 1, {"downsample": kwargs["downsample_size"]})
    else:
        raise ValueError("No such statistics is found.")


def heterozygosity(df):
    """
    The statistics function for Heterozygosity. Notice that
    \hat{h} = 1 - (N * \sum{ (p_l)^2 } - 1) / (N - 1)
            = N * (1 - \sum{ (p_l)^2 }) / (N - 1)
            = 2 * p_l * (1 - p_l) * N / (N - 1)                 (*)
    where p_l is the sample frequency of the l-th allele at the SNP locus and
    N is the total number of alleles at that site.
    (*) The last equation holds when there are only two possibilities, with
    probability p_l and (1 - p_l), respectively.

    Reference:
    (1) Ioannidis, A. G., Blanco-Portillo, J., Sandoval, K., Hagelberg, E.,
        Barberena-Jonas, C., Hill, A. V., ... & Moreno-Estrada, A. (2021).
        Paths and timings of the peopling of Polynesia inferred from genomic
        networks. Nature, 597(7877), 522-526.
    """
    # df.colname = ['populus1_frq', 'populus1_ct']
    temp = 2 * df[:, 0] * (1 - df[:, 0])
    heterozyg = np.mean(temp * (df[:, 1] / (df[:, 1] - 1)))
    return heterozyg


def f2(df):
    """
    The statistics function for F2. Notice that
    \hat{F2} = (p_A - p_B)^2 - (p_A * (1 - p_A) / (n_A - 1))
                             - (p_B * (1 - p_B) / (n_B - 1))
    where p_A is the sample allele frequency in the ancestry of interest at the
    SNP locus and n_A is the total observations at that site in population A.

    Reference:
    (1) Ioannidis, A. G., Blanco-Portillo, J., Sandoval, K., Hagelberg, E.,
        Barberena-Jonas, C., Hill, A. V., ... & Moreno-Estrada, A. (2021).
        Paths and timings of the peopling of Polynesia inferred from genomic
        networks. Nature, 597(7877), 522-526.
    """
    # df.colname = ['populus1_frq', 'populus1_ct', 'populus2_frq', 'populus2_ct']
    f2_biased = (df[:, 2] - df[:, 0]) ** 2
    adjusted_pop1 = df[:, 0] * (1 - df[:, 0]) / (df[:, 1] - 1)
    adjusted_pop2 = df[:, 2] * (1 - df[:, 2]) / (df[:, 3] - 1)
    f2 = np.mean(f2_biased - adjusted_pop1 - adjusted_pop2)
    return f2


def fst(df):
    """
    The statistics function for F_ST. Notice that
    \hat{F}_{ST} = ( (p_A - p_B)^2 - (p_A * (1 - p_A) / (n_A - 1))
                                   - (p_B * (1 - p_B) / (n_B - 1)) )
                 / ( p_A * (1 - p_B) + p_B * (1 - p_A) )
                 = \hat{F2} / \hat{\pi}.                 (*)
    where p_A is the sample allele frequency in the ancestry of interest at the
    SNP locus and n_A is the total observations at that site in population A.
    (*) For multiple sites, the numerator (F2) and the denominator (Pi) are
    averaged across all SNPs separately before taking the ratio to create a
    consistent estimator.

    Reference:
    (1) Ioannidis, A. G., Blanco-Portillo, J., Sandoval, K., Hagelberg, E.,
        Barberena-Jonas, C., Hill, A. V., ... & Moreno-Estrada, A. (2021).
        Paths and timings of the peopling of Polynesia inferred from genomic
        networks. Nature, 597(7877), 522-526.
    (2) Bhatia, G., Patterson, N., Sankararaman, S. & Price, A. L. Estimating
        and interpreting FST: the impact of rare variants. Genome Res. 23,
        1514–1521 (2013).
    """
    # df.colname = ['populus1_frq', 'populus1_ct', 'populus2_frq', 'populus2_ct']
    f2_biased = (df[:, 2] - df[:, 0]) ** 2
    adjusted_pop1 = df[:, 0] * (1 - df[:, 0]) / (df[:, 1] - 1)
    adjusted_pop2 = df[:, 2] * (1 - df[:, 2]) / (df[:, 3] - 1)
    f2 = np.sum(f2_biased - adjusted_pop1 - adjusted_pop2)
    pi = np.sum(df[:, 0] * (1 - df[:, 2]) + df[:, 2] * (1 - df[:, 0]))
    return f2 / pi


def f3(df):
    """
    The statistics function for F3. Notice that
    \hat{F3}(C; A, B) = (p_C - p_A) * (p_C - p_B) - \hat{h}_C / n_C
    where p_C is the sample allele frequency in the ancestry of interest at the
    SNP locus and n_C is the total observations at that site in population C.
    \hat{h}_C refers to the heterozygosity estimator at the specific SNP locus
    in population C.

    Reference:
    (1) Ioannidis, A. G., Blanco-Portillo, J., Sandoval, K., Hagelberg, E.,
        Barberena-Jonas, C., Hill, A. V., ... & Moreno-Estrada, A. (2021).
        Paths and timings of the peopling of Polynesia inferred from genomic
        networks. Nature, 597(7877), 522-526.
    """
    # df.colname = ['populus1_frq', 'populus2_frq', 'outgroup_frq', 'outgroup_ct']
    f3_biased = np.mean((df[:, 2] - df[:, 0]) * (df[:, 2] - df[:, 1]))
    adjusted = np.mean(2 * df[:, 2] * (1 - df[:, 2]) / (df[:, 3] - 1))
    return f3_biased - adjusted


def f4(df):
    """
    The statistics function for F4. Notice that
    \hat{F4}(A, B; C, D) = (p_A - p_B) * (p_C - p_D)
    where p_A is the sample allele frequency in the ancestry of interest at the
    SNP locus and n_A is the total observations at that site in population A.

    Reference:
    (1) Ioannidis, A. G., Blanco-Portillo, J., Sandoval, K., Hagelberg, E.,
        Barberena-Jonas, C., Hill, A. V., ... & Moreno-Estrada, A. (2021).
        Paths and timings of the peopling of Polynesia inferred from genomic
        networks. Nature, 597(7877), 522-526.
    """
    # df.colname = ['populus1_frq', 'populus2_frq', 'outgroup_frq', 'populus4_frq']
    f4 = np.mean((df[:, 0] - df[:, 1]) * (df[:, 2] - df[:, 3]))
    return f4


def pi(df):
    """
    The statistics function for Pi. Notice that
    \hat{\pi}(A,B) = p_A * (1 - p_B) + p_B * (1 - p_A)
    where p_A is the sample allele frequency in the ancestry of interest at the
    SNP locus and n_A is the total observations at that site in population A.

    Reference:
    (1) Ioannidis, A. G., Blanco-Portillo, J., Sandoval, K., Hagelberg, E.,
        Barberena-Jonas, C., Hill, A. V., ... & Moreno-Estrada, A. (2021).
        Paths and timings of the peopling of Polynesia inferred from genomic
        networks. Nature, 597(7877), 522-526.
    """
    # df.colname = ['populus1_frq', 'populus2_frq']
    pi = np.mean(df[:, 0] * (1 - df[:, 1]) + df[:, 1] * (1 - df[:, 0]))
    return pi


def psi(df, downsample=2):
    """
    The statistics function for Psi. Notice that
    \hat{\psi}(A,B) = p_A - p_B                         (*)
    where p_A is the sample allele frequency in the ancestry of interest at the
    SNP locus and n_A is the total observations at that site in population A.
    (*) For multiple sites, the Psi values are only averaged across all derived
    allele SNP loci. Each derived allele SNP locus is determined by the reference
    population C (polarization outgroup) where the sample minor allele frequency,
    i.e., min(p_C, 1 - p_C), at this site is not above some frequency threshold.
    See 'pop_aggr.py' for more implementation details.

    The code is deeply indebted to the methods described in
    (1) Peter, B.M. and Slatkin, M., 2013. Detecting range expansions from genetic
        data. Evolution, 67(11), pp.3274-3289.
    (2) Peter, B.M. and Slatkin, M., 2015. The effective founder effect in a
        spatially expanding population. Evolution, 69(3), pp.721-734.
    and the implementation of the methods in
    (1) https://github.com/BenjaminPeter/rangeexpansion

    Reference:
    (1) Ioannidis, A. G., Blanco-Portillo, J., Sandoval, K., Hagelberg, E.,
        Barberena-Jonas, C., Hill, A. V., ... & Moreno-Estrada, A. (2021).
        Paths and timings of the peopling of Polynesia inferred from genomic
        networks. Nature, 597(7877), 522-526.
    """
    # Downsampling size
    n = downsample

    # df.colname = ['populus1_frq', 'populus1_ct', 'populus2_frq', 'populus2_ct']
    cond = (df[:, 0] > 0) & (df[:, 2] > 0) & (df[:, 1] >= n) & (df[:, 3] >= n)
    df = df[cond, :]
    nrow = df.shape[0]
    if nrow == 0:
        return np.nan

    # poly_mat is the normalization constant for psi.
    # (Do not consider the cases when
    #  1. either population has downsampled SNPs with 0 count;
    #  2. both populations have downsampled SNPs with n counts.)
    #
    # B\A    0   1   2  ...  n-1  n
    #        ______________________
    # 0     |0   0   0  ...  0    0
    # 1     |0   1   1  ...  1    1
    # 2     |0   1   1  ...  1    1
    # ...   |
    # n     |0   1   1  ...  1    0
    poly_mat = np.zeros((n + 1, n + 1))
    poly_mat[1:, 1:] = 1
    poly_mat[-1, -1] = 0

    # psi_mat is the contribution to psi for each entry.
    #
    # B\A    0    1    2   ...  n-1    n
    #       ____________________________
    # 0    | 0    0    0   ...   0     0
    # 1    | 0    0    1   ...  n-2   n-1
    # 2    | 0   -1    0   ...  n-3   n-2
    # ...  |
    # n    | 0  -n+1 -n+2  ...   -1    0
    psi_mat = np.arange(n + 1).reshape((1, -1)) - np.arange(n + 1).reshape((-1, 1))
    psi_mat[0, :] = 0
    psi_mat[:, 0] = 0

    # Get the sum of 2D weights (in discrete proportions) from block of allele
    # positions.
    # B\A    0     1    ...    n
    #       ______________________
    # 0    |p_00  p_01  ...   p_0n
    # 1    |p_10  p_11  ...   p_1n
    # ...  |
    # n    |p_n0  p_n1  ...   p_nn
    rvs = hypergeom(df[:, [1, 3]], df[:, [0, 2]], n)  # (nrow,2)
    xs = np.arange(n + 1).repeat(2 * nrow).reshape((n + 1, nrow, 2))  # (n+1,nrow,2)
    xs_prob = rvs.pmf(xs)  # (n+1,nrow,2)
    # (nrow,n+1,1) x (nrow,1,n+1) -> (nrow,n+1,n+1) -> (n+1,n+1)
    resampled_mat = np.einsum("ji,ki->jk", xs_prob[:, :, 1], xs_prob[:, :, 0])
    return np.sum(resampled_mat * psi_mat) / np.sum(resampled_mat * poly_mat) / n


# +-*=+-*=+-*=+-*=+-*=+-*=+-*=+-*=+*-=+-*=+*-=+-*=+-*=+-*=+-*=+-*= #
#                        Statistics Inputs                         #
# +-*=+-*=+-*=+-*=+-*=+-*=+-*=+-*=+*-=+-*=+*-=+-*=+-*=+-*=+-*=+-*= #


@dataclass
class UnaryInputConfig:
    """
    The object input for a unary genetic statistics GS(A).
    """

    data_dir: str
    pop1_file: str
    block_size: int
    num_replicates: int
    output_file: str
    is_masked: bool


@dataclass
class BinaryInputConfig(UnaryInputConfig):
    """
    The object input for a binary genetic statistics GS(A, B).
    """

    pop2_file: str


@dataclass
class TernaryInputConfig(BinaryInputConfig):
    """
    The object input for a ternary genetic statistics GS(A, B, C).
    """

    pop3_file: str


@dataclass
class QuaternaryInputConfig(TernaryInputConfig):
    """
    The object input for a ternary genetic statistics GS(A, B, C, D).
    """

    pop4_file: str


@dataclass
class PsiConfig(BinaryInputConfig):
    """
    The object input for genetic statistics Psi, Psi(A, B).
    """

    da_file: str
    downsample_size: int


# +-*=+-*=+-*=+-*=+-*=+-*=+-*=+-*=+*-=+-*=+*-=+-*=+-*=+-*=+-*=+-*= #
#                         Statistics Objects                       #
# +-*=+-*=+-*=+-*=+-*=+-*=+-*=+-*=+*-=+-*=+*-=+-*=+-*=+-*=+-*=+-*= #


def read_freq(filename, is_masked):
    if is_masked:
        return pd.read_csv(
            filename, sep="\t", header=0, usecols=["MAF_MSK"], dtype={"MAF_MSK": float}
        )
    else:
        return pd.read_csv(
            filename, sep="\t", header=0, usecols=["MAF"], dtype={"MAF": float}
        )


def read_freq_and_ct(filename, is_masked):
    if is_masked:
        return pd.read_csv(
            filename,
            sep="\t",
            header=0,
            usecols=["MAF_MSK", "NCHROBS_MSK"],
            dtype={"MAF_MSK": float, "NCHROBS_MSK": float},
        )
    else:
        return pd.read_csv(
            filename,
            sep="\t",
            header=0,
            usecols=["MAF", "NCHROBS"],
            dtype={"MAF": float, "NCHROBS": float},
        )


class GeneticStatistics(ABC):
    """
    The object for a generic genetic statistics.
    """

    @abstractmethod
    def __init__(self, config):
        self.data_dir = config.data_dir
        self.block_size = int(config.block_size)
        self.num_replicates = int(config.num_replicates)
        self.output_file = self.absolute_path(config.output_file)
        self.is_masked = config.is_masked
        self.num_SNPs_thres = self.block_size * NUM_SNP_RATIO
        self.mean, self.std_error = 0.0, 0.0

    def absolute_path(self, f):
        """
        The function returns the absolute path of the given file.
        Args:
            f           The given file in relative path
        """
        return os.path.abspath(f)

    def absolute_path_pop_file(self, data_dir, f):
        """
        The function returns the absolute path of the given frequency file.
        Args:
            data_dir    The directory of frequency files
            f           The given frequency file in relative path
        """
        return os.path.join(self.absolute_path(data_dir), f)

    def get_pop_name(self, pop_file):
        """
        The function returns the population name by the frequency file name.
        Args:
            pop_file    The frequency file for one population
        """
        return (pop_file.strip().split("/")[-1]).split(".")[0]

    def get_mean(self):
        """
        The function returns the mean of the genetic statistics.
        """
        return self.mean

    def get_std_error(self):
        """
        The function returns the standard error of the genetic statistics.
        """
        return self.std_error

    @abstractmethod
    def filter_data(self, df):
        """
        The function filters the dataframe to be successfully passed into the
        statistics function.
        Args:
            df      The dataframe of SNP loci x population (frequency &
                        count) table
        """
        pass

    @abstractmethod
    def compute_statistics(self, df):
        """
        The function computes the genetic statistics.
        Args:
            df      The dataframe of SNP loci x population (frequency &
                        count) table
        """
        nrow = df.shape[0]
        self.no_compute_se = False
        if nrow < self.num_SNPs_thres:
            self.no_compute_se = True

    @abstractmethod
    def write_output(self, pop1, pop2, nrow):
        """
        The function writes the statistics output to the output text file.
        Args:
            pop1        The name of population 1
            pop2        The name of population 2
            nrow        The number of SNP loci
        """
        output = f"{pop1:s}-{pop2:s}\t{self.mean:.8f}\t{self.std_error:.8f}\t{nrow:d}\n"
        with open(self.output_file, "a") as f:
            f.write(output)


class Heterozygosity(GeneticStatistics):
    def __init__(self, config):
        super().__init__(config)
        self.pop1_file = self.absolute_path_pop_file(self.data_dir, config.pop1_file)

        df_pop1 = read_freq_and_ct(self.pop1_file, self.is_masked)
        df_pop1.columns = ["populus1_frq", "populus1_ct"]
        # ['populus1_frq', 'populus1_ct']
        df = self.filter_data(df_pop1)

        self.mean, self.std_error = self.compute_statistics(df)

        pop1_name = self.get_pop_name(self.pop1_file)
        nrow = df.shape[0]
        self.write_output(pop1_name, nrow)

    def filter_data(self, df):
        return df.dropna().query("populus1_ct > 1").to_numpy()

    def compute_statistics(self, df):
        super().compute_statistics(df)
        if self.no_compute_se:
            print("Not enough number of SNPs for block bootstrap!")
            return heterozygosity(df), np.nan
        else:
            remote_elements = [
                compute_statistics_helper.remote(df, self.block_size, "heterozygosity")
                for _ in range(self.num_replicates)
            ]
            sample_statistics = ray.get(remote_elements)
            sample_statistics = np.concatenate(sample_statistics, axis=None)
            return heterozygosity(df), np.std(sample_statistics, ddof=1)

    def write_output(self, pop1, nrow):
        output = f"{pop1:s}\t{self.mean:.8f}\t{self.std_error:.8f}\t{nrow:d}\n"
        with open(self.output_file, "a") as f:
            f.write(output)


class F2(GeneticStatistics):
    def __init__(self, config):
        super().__init__(config)
        self.pop1_file = self.absolute_path_pop_file(self.data_dir, config.pop1_file)
        self.pop2_file = self.absolute_path_pop_file(self.data_dir, config.pop2_file)

        df_pop1 = read_freq_and_ct(self.pop1_file, self.is_masked)
        df_pop1.columns = ["populus1_frq", "populus1_ct"]
        df_pop2 = read_freq_and_ct(self.pop2_file, self.is_masked)
        df_pop2.columns = ["populus2_frq", "populus2_ct"]
        # ['populus1_frq', 'populus1_ct', 'populus2_frq', 'populus2_ct']
        df = self.filter_data(pd.concat([df_pop1, df_pop2], axis=1))

        self.mean, self.std_error = self.compute_statistics(df)

        pop1_name = self.get_pop_name(self.pop1_file)
        pop2_name = self.get_pop_name(self.pop2_file)
        nrow = df.shape[0]
        self.write_output(pop1_name, pop2_name, nrow)

    def filter_data(self, df):
        return df.dropna().query("populus1_ct > 1 & populus2_ct > 1").to_numpy()

    def compute_statistics(self, df):
        super().compute_statistics(df)
        if self.no_compute_se:
            print("Not enough number of SNPs for block bootstrap!")
            return f2(df), np.nan
        else:
            remote_elements = [
                compute_statistics_helper.remote(df, self.block_size, "f2")
                for _ in range(self.num_replicates)
            ]
            sample_statistics = ray.get(remote_elements)
            sample_statistics = np.concatenate(sample_statistics, axis=None)
            return f2(df), np.std(sample_statistics, ddof=1)

    def write_output(self, pop1, pop2, nrow):
        super().write_output(pop1, pop2, nrow)


class Fst(F2):
    def __init__(self, config):
        super().__init__(config)

    def filter_data(self, df):
        return super().filter_data(df)

    def compute_statistics(self, df):
        super().compute_statistics(df)
        if self.no_compute_se:
            print("Not enough number of SNPs for block bootstrap!")
            return fst(df), np.nan
        else:
            remote_elements = [
                compute_statistics_helper.remote(df, self.block_size, "fst")
                for _ in range(self.num_replicates)
            ]
            sample_statistics = ray.get(remote_elements)
            sample_statistics = np.concatenate(sample_statistics, axis=None)
            return fst(df), np.std(sample_statistics, ddof=1)

    def write_output(self, pop1, pop2, nrow):
        super().write_output(pop1, pop2, nrow)


class F3(GeneticStatistics):
    def __init__(self, config):
        super().__init__(config)
        self.pop1_file = self.absolute_path_pop_file(self.data_dir, config.pop1_file)
        self.pop2_file = self.absolute_path_pop_file(self.data_dir, config.pop2_file)
        self.pop3_file = self.absolute_path_pop_file(self.data_dir, config.pop3_file)

        df_pop1 = read_freq(self.pop1_file, self.is_masked)
        df_pop1.columns = ["populus1_frq"]
        df_pop2 = read_freq(self.pop2_file, self.is_masked)
        df_pop2.columns = ["populus2_frq"]
        df_outgroup = read_freq_and_ct(self.pop3_file, self.is_masked)
        df_outgroup.columns = ["outgroup_frq", "outgroup_ct"]
        # ['populus1_frq', 'populus2_frq', 'outgroup_frq', 'outgroup_ct']
        df = self.filter_data(pd.concat([df_pop1, df_pop2, df_outgroup], axis=1))

        self.mean, self.std_error = self.compute_statistics(df)

        pop1_name = self.get_pop_name(self.pop1_file)
        pop2_name = self.get_pop_name(self.pop2_file)
        pop3_name = self.get_pop_name(self.pop3_file)
        nrow = df.shape[0]
        self.write_output(pop1_name, pop2_name, pop3_name, nrow)

    def filter_data(self, df):
        return df.dropna().query("outgroup_ct > 1").to_numpy()

    def compute_statistics(self, df):
        super().compute_statistics(df)
        if self.no_compute_se:
            print("Not enough number of SNPs for block bootstrap!")
            return f3(df), np.nan
        else:
            remote_elements = [
                compute_statistics_helper.remote(df, self.block_size, "f3")
                for _ in range(self.num_replicates)
            ]
            sample_statistics = ray.get(remote_elements)
            sample_statistics = np.concatenate(sample_statistics, axis=None)
            return f3(df), np.std(sample_statistics, ddof=1)

    def write_output(self, pop1, pop2, pop3, nrow):
        output = (
            f"F3({pop1:s}, {pop2:s}; {pop3:s})\t{self.mean:.8f}\t"
            f"{self.std_error:.8f}\t{nrow:d}\n"
        )
        with open(self.output_file, "a") as f:
            f.write(output)


class F4(GeneticStatistics):
    def __init__(self, config):
        super().__init__(config)
        self.pop1_file = self.absolute_path_pop_file(self.data_dir, config.pop1_file)
        self.pop2_file = self.absolute_path_pop_file(self.data_dir, config.pop2_file)
        self.pop3_file = self.absolute_path_pop_file(self.data_dir, config.pop3_file)
        self.pop4_file = self.absolute_path_pop_file(self.data_dir, config.pop4_file)

        df_pop1 = read_freq(self.pop1_file, self.is_masked)
        df_pop1.columns = ["populus1_frq"]
        df_pop2 = read_freq(self.pop2_file, self.is_masked)
        df_pop2.columns = ["populus2_frq"]
        df_outgroup = read_freq(self.pop3_file, self.is_masked)
        df_outgroup.columns = ["outgroup_frq"]
        df_pop4 = read_freq(self.pop4_file, self.is_masked)
        df_pop4.columns = ["populus4_frq"]
        # ['populus1_frq', 'populus2_frq', 'outgroup_frq', 'populus4_frq']
        df = self.filter_data(
            pd.concat([df_pop1, df_pop2, df_outgroup, df_pop4], axis=1)
        )

        self.mean, self.std_error = self.compute_statistics(df)

        pop1_name = self.get_pop_name(self.pop1_file)
        pop2_name = self.get_pop_name(self.pop2_file)
        pop3_name = self.get_pop_name(self.pop3_file)
        pop4_name = self.get_pop_name(self.pop4_file)
        nrow = df.shape[0]
        self.write_output(pop1_name, pop2_name, pop3_name, pop4_name, nrow)

    def filter_data(self, df):
        return df.dropna().to_numpy()

    def compute_statistics(self, df):
        super().compute_statistics(df)
        if self.no_compute_se:
            print("Not enough number of SNPs for block bootstrap!")
            return f4(df), np.nan
        else:
            remote_elements = [
                compute_statistics_helper.remote(df, self.block_size, "f4")
                for _ in range(self.num_replicates)
            ]
            sample_statistics = ray.get(remote_elements)
            sample_statistics = np.concatenate(sample_statistics, axis=None)
            return f4(df), np.std(sample_statistics, ddof=1)

    def write_output(self, pop1, pop2, pop3, pop4, nrow):
        output = (
            f"F4({pop1:s}, {pop2:s}; {pop3:s}, {pop4:s})\t{self.mean:.8f}\t"
            f"{self.std_error:.8f}\t{nrow:d}\n"
        )
        with open(self.output_file, "a") as f:
            f.write(output)


class Pi(GeneticStatistics):
    def __init__(self, config):
        super().__init__(config)
        self.pop1_file = self.absolute_path_pop_file(self.data_dir, config.pop1_file)
        self.pop2_file = self.absolute_path_pop_file(self.data_dir, config.pop2_file)

        df_pop1 = read_freq(self.pop1_file, self.is_masked)
        df_pop1.columns = ["populus1_frq"]
        df_pop2 = read_freq(self.pop2_file, self.is_masked)
        df_pop2.columns = ["populus2_frq"]
        # ['populus1_frq', 'populus2_frq']
        df = self.filter_data(pd.concat([df_pop1, df_pop2], axis=1))

        self.mean, self.std_error = self.compute_statistics(df)

        pop1_name = self.get_pop_name(self.pop1_file)
        pop2_name = self.get_pop_name(self.pop2_file)
        nrow = df.shape[0]
        self.write_output(pop1_name, pop2_name, nrow)

    def filter_data(self, df):
        return df.dropna().to_numpy()

    def compute_statistics(self, df):
        super().compute_statistics(df)
        if self.no_compute_se:
            print("Not enough number of SNPs for block bootstrap!")
            return pi(df), np.nan
        else:
            remote_elements = [
                compute_statistics_helper.remote(df, self.block_size, "pi")
                for _ in range(self.num_replicates)
            ]
            sample_statistics = ray.get(remote_elements)
            sample_statistics = np.concatenate(sample_statistics, axis=None)
            return pi(df), np.std(sample_statistics, ddof=1)

    def write_output(self, pop1, pop2, nrow):
        super().write_output(pop1, pop2, nrow)


class Psi(GeneticStatistics):
    def __init__(self, config):
        super().__init__(config)
        self.pop1_file = self.absolute_path_pop_file(self.data_dir, config.pop1_file)
        self.pop2_file = self.absolute_path_pop_file(self.data_dir, config.pop2_file)
        self.da_file = self.absolute_path_pop_file(self.data_dir, config.da_file)
        self.downsample_size = int(config.downsample_size)

        df_pop1 = pd.read_csv(self.pop1_file, sep="\t", header=0)
        df_pop1 = read_freq_and_ct(self.pop1_file, self.is_masked)
        df_pop1.columns = ["populus1_frq", "populus1_ct"]
        df_pop2 = read_freq_and_ct(self.pop2_file, self.is_masked)
        df_pop2.columns = ["populus2_frq", "populus2_ct"]
        ref_thres = pd.read_csv(
            self.da_file,
            header=None,
            skiprows=1,
            names=["ref_pop_DA_marker"],
            dtype={"ref_pop_DA_marker": int},
        )
        # ['populus1_frq', 'populus1_ct', 'populus2_frq', 'populus2_ct']
        df = self.filter_data(pd.concat([df_pop1, df_pop2, ref_thres], axis=1))

        self.mean, self.std_error = self.compute_statistics(df)

        pop1_name = self.get_pop_name(self.pop1_file)
        pop2_name = self.get_pop_name(self.pop2_file)
        cond = (
            (df[:, 0] > 0)
            & (df[:, 2] > 0)
            & (df[:, 1] >= self.downsample_size)
            & (df[:, 3] >= self.downsample_size)
        )
        df = df[cond, :]
        nrow = df.shape[0]
        self.write_output(pop1_name, pop2_name, nrow)

    def filter_data(self, df):
        df = df.dropna()[df["ref_pop_DA_marker"] != 0].to_numpy()
        df[:, [0, 2]] = np.where(
            df[:, [-1, -1]] == 2, 1.0 - df[:, [0, 2]], df[:, [0, 2]]
        )
        df[:, [0, 2]] = np.rint(df[:, [0, 2]] * df[:, [1, 3]])
        return df[:, :-1]

    def compute_statistics(self, df):
        super().compute_statistics(df)
        if self.no_compute_se:
            print("Not enough number of SNPs for block bootstrap!")
            return psi(df), np.nan
        else:
            remote_elements = [
                compute_statistics_helper.remote(
                    df, self.block_size, "psi", downsample_size=self.downsample_size
                )
                for _ in range(self.num_replicates)
            ]
            sample_statistics = ray.get(remote_elements)
            sample_statistics = np.concatenate(sample_statistics, axis=None)
            return psi(df), np.std(sample_statistics, ddof=1)

    def write_output(self, pop1, pop2, nrow):
        super().write_output(pop1, pop2, nrow)


if __name__ == "__main__":
    # config = BinaryInputConfig(data_dir="data_American",
    #                            block_size="200",
    #                            num_replicates="200",
    #                            output_file="output",
    #                            pop1_file="CookAtiu.freq",
    #                            pop2_file="CookMauke.freq")
    # config2 = TrinaryInputConfig(data_dir="data_American",
    #                            block_size="200",
    #                            num_replicates="200",
    #                            output_file="output",
    #                            pop1_file="CookAtiu.freq",
    #                            pop2_file="CookMauke.freq",
    #                            pop3_file='Atayal.freq')
    config3 = PsiConfig(
        data_dir="data/cluster2",
        block_size="200",
        num_replicates="100",
        output_file="output",
        pop1_file="CookAtiu.freq",
        pop2_file="CookMauke.freq",
        da_file="psi_atbupa_frq_50.txt",
        downsample_size=2,
    )
    it = time.time()
    # h = Pi(config)
    # h = F3(config2)
    h = Psi(config3)

    print(f"time:{time.time() - it:f}")
    print(f"mean:{h.get_mean()}, std_error:{h.get_std_error()}")
    print(h.pop1_file, h.block_size, h.num_replicates, h.output_file)
