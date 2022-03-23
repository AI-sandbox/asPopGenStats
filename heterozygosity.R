#!/usr/bin/env Rscript

# Compute the heterozygosity for one (Oceanian) population groups A
# (represented by 1 in code).

# Import packages.
suppressPackageStartupMessages({
   library(boot)
   library(parallel)
   library(tictoc)
   library(stringr)
})

args = commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
   stop(paste("Usage: Rscript heterozygosity.R <population1.freq>",
              "<block_size (int)> <num_replicates (int)> <data_directory> <output_file>"),
        call. = FALSE)
}
POP1_FILE = args[1]
BLOCK_SIZE = strtoi(args[2])
NUM_REPLICATES = strtoi(args[3])
DATA_DIR = args[4]
OUTPUT_FILE = args[5]

# Set working directory.
setwd(paste0(getwd(), '/', DATA_DIR))

# Take "FREQ" and "CT" from the ".freq" files of the target population (A).
df_populus1 = read.table(POP1_FILE, header = TRUE, sep = ",",
                         na.strings = "NaN",
                         colClasses = c(rep("NULL", 1), rep("numeric", 2)) )
freq_series = ts(df_populus1)
colnames(freq_series) <- c("populus1", "populus1_ct")
# Exclude all NaN value SNPs.
freq_series <- na.omit(freq_series)
# Exclude all SNPs that has n_A <= 1 
freq_series <- freq_series[(freq_series[, 2] > 1), ]

# Create a function that computes pi statistics.
Heterozygosity <- function(frq_series) {
   num_of_SNPs <- nrow(frq_series)

   temp <- 2 * (frq_series[, 1] * (1 - frq_series[, 1]))
   Heterozyg <- sum( temp * (frq_series[, 2] / (frq_series[, 2] - 1)) )
   return(Heterozyg / num_of_SNPs)
}

# Start block bootstrap.
num_replicates = NUM_REPLICATES

tic("heterozygosity")
if (nrow(freq_series) >= BLOCK_SIZE) {
   cl <- makeCluster(4)
   clusterExport(cl = cl, c('freq_series'))
   booted_heterozyg = tsboot(freq_series, Heterozygosity, R = num_replicates, l = BLOCK_SIZE,
                             sim = "fixed", endcorr = TRUE, parallel = "snow",
                             ncpus = min(c(detectCores(), 10 )), cl = cl )
   heterozyg <- booted_heterozyg$t0 # Heterozygosity matrix
} else {
   heterozyg <- NaN
}
toc()
cat(heterozyg)

# Check the total number of valid (i.e., not-NA) SNPS used for pi statistics
# and block bootstrap.
std_error_heterozyg = 0
# If the number of valid SNPs is less than 1000, we will not calculate the
# standard deviation of heterozygosity from block bootstrap.
if (nrow(freq_series) >= 1000) {
   std_error_heterozyg = apply(booted_heterozyg$t, 2, sd, na.rm = TRUE)
}

# Output format:
#   PopulationA    heterozygosity    heterozygosity_stat_std    num_of_valid_SNPs
output_data <- paste0(str_replace(args[1], ".freq", ""), "\t", heterozyg, "\t",
                      std_error_heterozyg, "\t", nrow(freq_series))
write(output_data, file = OUTPUT_FILE, append = TRUE)
