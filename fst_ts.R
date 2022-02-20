#!/usr/bin/env Rscript

# Compute F_ST statistics for the given two (populusn) population groups A and B
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
   stop(paste("Usage: Rscript fst_ts.R <population1.freq> <population2.freq>",
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

# Take "FREQ" from the ".freq" files of three target populations (A, B and
# outgroup O) and "CT" column only from the outgroup O's ".freq" file.
df_populus1 = read.table(args[1], header = TRUE, sep = ",",
                         na.strings = "NaN",
                         colClasses = c(rep("NULL", 1), rep("numeric", 2)) )
df_populus2 = read.table(args[2], header = TRUE, sep = ",",
                         na.strings = "NaN",
                         colClasses = c(rep("NULL", 1), rep("numeric", 2)) )
freq_series = ts(cbind(df_populus1, df_populus2))
colnames(freq_series) <- c("populus1", "populus1_ct", "populus2", "populus2_ct")
# Exclude all NaN value SNPs.
freq_series <- na.omit(freq_series)
# Exclude all SNPs that has p_A = p_B = 0 (i.e., pi_stat = 0)
freq_series <- freq_series[(freq_series[, 1] != 0) | (freq_series[, 3] != 0), ]

# Create a function that computes F_ST statistics.
fst <- function(frq_series) {
   
   biased_f2 <- ((frq_series[, 3] - frq_series[, 1]) ** 2) # array
   adjusted_a <- (frq_series[, 1] * (1 - frq_series[, 1])) /
                 ((frq_series[, 2] - 1))
   adjusted_b <- (frq_series[, 3] * (1 - frq_series[, 3])) /
                 ((frq_series[, 4] - 1))
   Pi <- (frq_series[, 1] * (1 - frq_series[, 3])) +
         ((1 - frq_series[, 1]) * frq_series[, 3])
   adjusted <- (biased_f2 - adjusted_a - adjusted_b)
   
   return( mean(Pi) )
}

# Start block bootstrap.
num_replicates = 100

tic("F_ST")
cl <- makeCluster(4)
clusterExport(cl = cl, c('freq_series'))
booted_fst = tsboot(freq_series, fst, R = num_replicates, l = BLOCK_SIZE,
                    sim = "fixed", endcorr = TRUE, parallel = "snow",
		              ncpus = min(c(detectCores(), 10 )), cl = cl )

toc()
fst_stat <- booted_fst$t0 # F_ST matrix
cat(fst_stat)

# Check the total number of valid (i.e., not-NA) SNPS used for F_ST statistics and
# block bootstrap.
std_error_fst = 0
# If the number of valid SNPs is less than 1000, we will not calculate the standard
# deviation of F_ST statistics from block bootstrap.
if (nrow(freq_series) >= 1000) {
   std_error_fst = apply(booted_fst$t, 2, sd, na.rm = TRUE)
}

# Output format
#   PopulationA-PopulationB    F_ST_stat    F_ST_stat_std    num_of_valid_SNPs
output_data <- paste0(str_replace(args[1], ".freq", ""), "-",
                      str_replace(args[2], ".freq", ""), "\t", fst_stat, "\t",
                      std_error_fst, "\t", nrow(freq_series))
write(output_data, file = OUTPUT_FILE, append = TRUE)
