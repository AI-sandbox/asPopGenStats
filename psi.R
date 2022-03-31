#!/usr/bin/env Rscript

# Compute psi statistics for the given two (Oceanian) population groups A and B
# (represented by 1 and 2 in code).
#
# The code is deeply indebted to the methods described in
# (1) Peter, B.M. and Slatkin, M., 2013. Detecting range expansions from genetic
#     data. Evolution, 67(11), pp.3274-3289.
# (2) Peter, B.M. and Slatkin, M., 2015. The effective founder effect in a
#     spatially expanding population. Evolution, 69(3), pp.721-734.
# and the implementation of the methods in
# (1) https://github.com/BenjaminPeter/rangeexpansion

# Import packages.
suppressPackageStartupMessages({
   library(boot)
   library(parallel)
   library(tictoc)
   library(stringr)
})

args = commandArgs(trailingOnly = TRUE)
if (length(args) < 7) {
   stop(paste("Usage: Rscript psi.R <population1.freq> <population2.freq>",
              "<block_size (int)> <num_replicates (int)> <data_directory> <output_file>",
              "<derived_allele_SNP_position_file> <downsampling_size (int)>"),
        call. = FALSE)
}
POP1_FILE = args[1]
POP2_FILE = args[2]
BLOCK_SIZE = strtoi(args[3])
NUM_REPLICATES = strtoi(args[4])
DATA_DIR = args[5]
OUTPUT_FILE = args[6]
DA_FILE = args[7]
DOWNSAMPLE_SIZE = strtoi(args[8])

# Set working directory.
setwd(paste0(getwd(), '/', DATA_DIR))

# Set up a derived-allele indicator file taken from a reference population group
# (e.g., Taiwan indigenous groups) with frequency threshold (e.g., 0, 0.05, or 0.1).
#
# A derived-allele indicator file should have the same length as those of the
# columns in other ".freq" files. Each row can only be "0", "1", or "2", where
# "1" indicates the alternative allele frequency is <= threshold and
# "2" indicates the alternative allele frequency is >= 1 - threshold.
#
ref_thresh = read.table(DA_FILE, header = FALSE, skip = 1)

# Take "FREQ" and "CT" columns from the ".freq" files of two target populations.
df_populus1 = read.table(POP1_FILE, header = TRUE, sep = ",",
                         na.strings="NaN",
                         colClasses = c(rep("NULL", 1), rep("numeric", 2)) )
df_populus2 = read.table(POP2_FILE, header = TRUE, sep = ",",
                         na.strings="NaN",
                         colClasses = c(rep("NULL", 1), rep("numeric", 2)) )

freq_series = ts(cbind(df_populus1, df_populus2, ref_thresh))
# Find all non-NaN SNPs with the reference population's derived alleles.
freq_series = freq_series[ref_thresh != 0 &
                          rowSums(is.nan(freq_series[, c(1, 3)])) == 0, ]
# Switch the target populations allele frequency from alpha to 1 - alpha
# based on the reference group's allele frequency:
# if at certain position, the allele frequency for the reference group is
# alpha_ref >= 1 - thresh, we switch the allele frequency of the two target
# populations alpha_1, alpha_2 to 1 - alpha_1, 1 - alpha_2, respectively.
freq_series[, 1] = ifelse(freq_series[, 5] == 2, 1 - freq_series[, 1],
                          freq_series[, 1])
freq_series[, 3] = ifelse(freq_series[, 5] == 2, 1 - freq_series[, 3],
                          freq_series[, 3])
freq_series = freq_series[, 1:4]
colnames(freq_series) <- c("pop1_drv_frq", "pop1_obs_ct", "pop2_drv_frq",
                           "pop2_obs_ct")

## Change derived allele frequency (DAF) to derived allele count (DAC).
freq_series[, 1] = round(freq_series[, 1] * freq_series[, 2])
freq_series[, 3] = round(freq_series[, 3] * freq_series[, 4])


# Create a function that computes psi statistics.
psi <- function(frq_series, downsample = 2) {
   # Downsampling size
   n = downsample
    
   tbl <- as.data.frame( frq_series )
   colnames(tbl) <- c("fi", "ni", "fj", "nj")
    
   to.exclude <- tbl$fi == 0 | tbl$fj == 0 |
    	           tbl$ni < n | tbl$nj < n
   tbl <- tbl[!to.exclude, ]
   if (nrow(tbl) == 0) { return(NaN) }

   # poly.mat is the normalization constant for psi-stat.
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
   poly.mat <- matrix(0, nrow = n+1, ncol = n+1)
   poly.mat[2:(n+1), 2:(n+1)] <- 1
   poly.mat[n+1, n+1] <- 0

   # psi.mat is the contribution to psi for each entry.
   #
   # B\A    0    1    2   ...  n-1    n
   #       ____________________________
   # 0    | 0    0    0   ...   0     0
   # 1    | 0    0    1   ...  n-2   n-1
   # 2    | 0   -1    0   ...  n-3   n-2
   # ...  | 
   # n    | 0  -n+1 -n+2  ...   -1    0
   psi.mat <- outer(0:n, 0:n, FUN = function(x,y)(y-x))
   psi.mat[1, ] <- 0
   psi.mat[, 1] <- 0

   f.contribution <- function(row, b) {
      a <- 0:b
      f1 <- row[1]
      n1 <- row[2]
      f2 <- row[3]
      n2 <- row[4]

      q1 <- dhyper(a, f1, n1 - f1, b)
      q2 <- dhyper(a, f2, n2 - f2, b)

      return(outer(q2, q1))
   }
   
   # Get the sum of 2D weights (in discrete proportions) from block of allele
   # positions.
   # B\A    0     1    ...    n
   #       ______________________
   # 0    |p_00  p_01  ...   p_0n
   # 1    |p_10  p_11  ...   p_1n
   # ...  |
   # n    |p_n0  p_n1  ...   p_nn
   resampled.mat <- matrix(rowSums(apply(tbl, 1, f.contribution, b = n)),
                           nrow = n+1)
   return( sum(resampled.mat * psi.mat) / sum(resampled.mat * poly.mat) / n )
}


# Start block bootstrap.
num_replicates = NUM_REPLICATES

# Get a sample output from psi function.
print(psi(freq_series, DOWNSAMPLE_SIZE))

tic("psi")
if (nrow(freq_series) >= BLOCK_SIZE) {
   cl <- makeCluster(4)
   clusterExport(cl = cl, c('freq_series'))
   booted_psi = tsboot(freq_series, psi, R = num_replicates, l = BLOCK_SIZE,
                     sim = "fixed", endcorr = TRUE, downsample = DOWNSAMPLE_SIZE,
                     parallel = "snow", ncpus = min(c(detectCores(), 10 )),
                     cl = cl )
   psi_stat <- booted_psi$t0 # PSI matrix
} else{
   psi_stat <- NaN
}
toc()
cat(psi_stat)

# Check the total number of SNPs used for psi statistics and block bootstrap.
std_error_psi = 0
num_of_snp <- sum(rowSums(is.nan(freq_series[, c(1, 3)])) == 0) # &
                  #  freq_series[, 1] >= 1 & freq_series[, 3] >= 1 &
                  #  freq_series[, 2] < downsample_size &
                  #  freq_series[, 4] < downsample_size)
if (num_of_snp >= 1000) {
   std_error_psi = apply(booted_psi$t, 2, sd, na.rm = TRUE)
}
output_data <- paste0(str_replace(args[1], ".freq", ""), "-",
                      str_replace(args[2], ".freq", ""), "\t",
                      ifelse(is.na(psi_stat), NA, psi_stat), "\t",
                      std_error_psi, "\t", num_of_snp)
write(output_data, file = OUTPUT_FILE, append = TRUE)
