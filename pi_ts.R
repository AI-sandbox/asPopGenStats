#!/usr/bin/env Rscript

# Compute pi statistics for the given two (Oceanian) population groups A and B
# (represented by 1 and 2 in code).

# Import packages.
suppressPackageStartupMessages({
   library(boot)
   library(parallel)
   library(tictoc)
   library(stringr)
})

args = commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
   stop(paste("Usage: Rscript pi_ts.R <population1.freq> <population2.freq>",
              "<block_size (int)> <data_directory> <output_file>"),
        call. = FALSE)
}
POP1_FILE = args[1]
POP2_FILE = args[2]
BLOCK_SIZE = strtoi(args[3])
DATA_DIR = args[4]
OUTPUT_FILE = args[5]

# Set working directory.
setwd(paste0(getwd(), '/', DATA_DIR))

# Take "FREQ" from the ".freq" files of two target populations (A and B).
df_populus1 = read.table(POP1_FILE, header = TRUE, sep = ",",
                         na.strings = "NaN",
                         colClasses = c(rep("NULL", 1), rep("numeric", 1),
		                                  rep("NULL", 1)) )
df_populus2 = read.table(POP2_FILE, header = TRUE, sep = ",",
                         na.strings = "NaN",
                         colClasses = c(rep("NULL", 1), rep("numeric", 1),
                                        rep("NULL", 1)) )
freq_series = ts(cbind(df_populus1, df_populus2))
colnames(freq_series) <- c("populus1", "populus2")
# Exclude all NaN value SNPs.
freq_series <- na.omit(freq_series)

# Create a function that computes pi statistics.
Pi <- function(frq_series) {
   num_of_SNPs <- nrow(frq_series)

   PI <- sum( (frq_series[, 1] * (1 - frq_series[, 2])) +
              ((1 - frq_series[, 1]) * frq_series[, 2]) )
   return(PI / num_of_SNPs)
}

# Start block bootstrap.
num_replicates = 100

tic("pi")
cl <- makeCluster(4)
clusterExport(cl = cl, c('freq_series'))
booted_pi = tsboot(freq_series, Pi, R = num_replicates, l = BLOCK_SIZE,
                   sim = "fixed", endcorr = TRUE, parallel = "snow",
                   ncpus = min(c(detectCores(), 10 )), cl = cl )
toc()
pi_stat <- booted_pi$t0 # PI matrix
cat(pi_stat)

# Check the total number of valid (i.e., not-NA) SNPS used for pi statistics
# and block bootstrap.
std_error_pi = 0
# If the number of valid SNPs is less than 1000, we will not calculate the
# standard deviation of f3 statistics from block bootstrap.
if (nrow(freq_series) >= 1000) {
   std_error_pi = apply(booted_pi$t, 2, sd, na.rm = TRUE)
}

# Output format:
#   PopulationA-PopulationB    pi_stat    pi_stat_std    num_of_valid_SNPs
output_data <- paste0(str_replace(args[1], ".freq", ""), "-",
                      str_replace(args[2], ".freq", ""), "\t", pi_stat, "\t",
                      std_error_pi, "\t", nrow(freq_series))
write(output_data, file = OUTPUT_FILE, append = TRUE)
