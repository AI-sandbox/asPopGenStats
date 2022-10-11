# asPopGenStats

Currently, ancestry-specific versions of “F3”, “F_st”, “Pi”, “Psi”, and Heterozygosity are supported.

usage: main.py [-h] -f FILE [FILE ...] -t DATA_DIR [-b BLOCKSIZE] [-n NUM_REPLICATES]  
               [-r GROUP [GROUP ...]] [-D DAF] [-d DOWNSAMPLE_SIZE] [--rm_DA_files RM_DA_FILES]  
               {F3,Pi,Psi,F_ST,H}  

positional arguments:  
  {F3,Pi,Psi,F_ST,H}    the statistic to compute (H stands for heterozygosity)

optional arguments:  
  -h, --help            show this help message and exit  
  -f FILE [FILE ...], --file FILE [FILE ...]  
                        (REQUIRED) files containing a list of population names for computation,  
                        each file has one column of names, only one file is needed for statistics that are computed  
                        on all pairs of populations (complete graph), two files are needed for statistics that  
                        are computed between only particular pairs of populations (bipartite graph)*  
  -t DATA_DIR, --data_dir DATA_DIR  
                        (REQUIRED) directory containing the snp frequency file**  
  -b BLOCKSIZE, --blocksize BLOCKSIZE  
                        block size for block bootstrap (number of consecutive snps), default=500  
  -n NUM_REPLICATES, --num_replicates NUM_REPLICATES  
                        number of replicates for block bootstrap, default=100  
  -r GROUP [GROUP ...], --ref_group GROUP [GROUP ...]  
                        (F3/Psi only) reference population group for F3/Psi statistics, if  
                        multiple population group names, create a single aggregated reference out of them  
  -D DAF, --DAF DAF     (Psi only) derived allele frequency threshold, (a float number between  
                        [0.0, 0.5])*** , default=-1.0 (i.e. indicates non-Psi stats are being computed)  
  -d DOWNSAMPLE_SIZE, --downsample_size DOWNSAMPLE_SIZE  
                        (Psi only) downsampling size for computing Psi statistics †, default=2  
  --rm_DA_files RM_DA_FILES  
                        (Boolean) whether to will remove the intermediary files that indicate SNP positions  
                        for derived alleles at the end of each run, default=True  

\* all pair example (complete graph): [A1, A2, A3] => A1A2, A1A3, A2A3, cross-match pair example (bipartite graph): [A1, A2], [B1, B2] => A1B1, A1B2, A2B1, A2B2

\*\* in the file directory, there should be frequency files named as \<pop>.freq, for example, “Samoa.freq”, which has the first line denoting the column names of the csv data file (comma-separated format) - CHROM_IDX, FREQ, CT - where “CHROM_IDX” is the label for the column containing chromosome number, “FREQ” is the label for the column containing the alternative allele frequency, and “CT” is the label for the column containing the total number of observations for that SNP across all samples in the population (accounting for samples that are missing or masked at that position).
  
\*\*\* when DAF=0.5, the Psi computation will have no polarization

† be aware that a valid downsampling size should be less than the minimum number of observations for the population groups passed into the statistics

An example of the execution code:
python3 main.py Psi -f pop_list.txt -t data_American -b 50 -n 200 -r Samoa Tonga -D 0.05 -d 2 --rm_DA_files True

(Psi only) The program will first generate an aggregated population file if there are multiple reference population groups and the file name will contain the first two letters of each reference group capitalized. So for the example above, you can find a file “SaTo.freq” in the frequency data file directory. Then it will generate a file specifying the derived allele positions, which is named in the following format: 
psi_<aggr_name>\_frq\_\<DAF>.txt

The outputs of the given statistics are in the directory: \<stat>\_output. The statistics matrix file, \<stat>\_mtx.csv, contains a matrix of the sample mean of the given statistics, where rows and columns represent “pop_list1.txt” and “pop_list2.txt”, respectively. The statistics output file, \<stat>\_stat.txt, contains a list of rows of outputs:  
        popA-popB        \<stat>\_mean        \<stat>\_SE        num_of_SNPs_used  
NOTE: “\<stat>\_mtx.csv” is overwritten at each run, while “\<stat>\_stat.txt” uses appending format for each run.

Cite: Ioannidis, A. G., Blanco-Portillo, J., Sandoval, K., Hagelberg, E., Barberena-Jonas, C., Hill, A. V., ... & Moreno-Estrada, A. (2021). Paths and timings of the peopling of Polynesia inferred from genomic networks. Nature, 597(7877), 522-526.
