#!/usr/bin/env Rscript

# Compute F3 statistics for the given two (Oceanian) population groups A and B
# (represented by 1 and 2 in code) with a given population outgroup O.

# Import packages.
suppressPackageStartupMessages({
   library(boot)
   library(parallel)
   library(tictoc)
   library(stringr)
})

args = commandArgs(trailingOnly = TRUE)
if (length(args) < 6) {
   stop(paste("Usage: Rscript f3_ts.R <population1.freq> <population2.freq>",
              "<block_size (int)> <data_directory> <output_file>",
              "<outgroup.freq>"),
        call. = FALSE)
}
POP1_FILE = args[1]
POP2_FILE = args[2]
BLOCK_SIZE = strtoi(args[3])
DATA_DIR = args[4]
OUTPUT_FILE = args[5]
OUTGROUP_FILE = args[6]

# Set working directory.
setwd(paste0(getwd(), '/', DATA_DIR))

# Take "FREQ" from the ".freq" files of three target populations (A, B and
# outgroup O) and "CT" column only from the outgroup O's ".freq" file.
df_populus1 = read.table(POP1_FILE, header = TRUE, sep = ",",
                         na.strings = "NaN",
                         colClasses = c(rep("NULL", 1), rep("numeric", 1),
		                                  rep("NULL", 1)) )
df_populus2 = read.table(POP2_FILE, header = TRUE, sep = ",",
                         na.strings = "NaN",
                         colClasses = c(rep("NULL", 1), rep("numeric", 1),
                                        rep("NULL", 1)) )
df_outgroup = read.table(OUTGROUP_FILE, header = TRUE, sep = ",",
                         na.strings = "NaN",
                         colClasses = c(rep("NULL", 1), rep("numeric", 2)) )
freq_series = ts(cbind(df_populus1, df_populus2, df_outgroup))
colnames(freq_series) <- c("populus1", "populus2", "outgroup", "outgroup_ct")
# Exclude all NaN value SNPs.
freq_series <- na.omit(freq_series)

# Create a function that computes F3 statistics.
f3 <- function(frq_series) {
   
   biased_f3 <- sum( ((frq_series[, 3] - frq_series[, 2]) *
                      (frq_series[, 3] - frq_series[, 1])) )
   temp <- (frq_series[, 3] * (1 - frq_series[, 3])) /
            (frq_series[, 4] - 1)
   heterozygosity <- 2 * sum(temp * frq_series[, 4]) 
   adjusted <- sum(temp)
   
   return( (biased_f3 - adjusted) / nrow(frq_series) )
}

# Start block bootstrap.
num_replicates = 100

tic("F3")
if nrow(freq_series) >= BLOCK_SIZE {
   cl <- makeCluster(4)
   clusterExport(cl = cl, c('freq_series'))
   booted_f3 = tsboot(freq_series, f3, R = num_replicates, l = BLOCK_SIZE,
                     sim = "fixed", endcorr = TRUE, parallel = "snow",
                     ncpus = min(c(detectCores(), 10 )), cl = cl )
   f3_stat <- booted_f3$t0 # F3 matrix
} else {
   f3_stat <- NaN
}
toc()
cat(f3_stat)

# Check the total number of valid (i.e., not-NA) SNPS used for F3 statistics and
# block bootstrap.
std_error_f3 = 0
# If the number of valid SNPs is less than 1000, we will not calculate the standard
# deviation of F3 statistics from block bootstrap.
if (nrow(freq_series) >= 1000) {
   std_error_f3 = apply(booted_f3$t, 2, sd, na.rm = TRUE)
}

# Output format
#   PopulationA-PopulationB    F3_stat    F3_stat_std    num_of_valid_SNPs
output_data <- paste0(str_replace(args[1], ".freq", ""), "-",
                      str_replace(args[2], ".freq", ""), "\t", f3_stat, "\t",
                      std_error_f3, "\t", nrow(freq_series))
write(output_data, file = OUTPUT_FILE, append = TRUE)
