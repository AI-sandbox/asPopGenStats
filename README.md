# asPopGenStats


## First Stage
`gnomix_frq.py`: a conversion that aggregates information from Gnomix output files (`.msp`) and vcf files (`.vcf`) and creates frequency files (`.freq`) for population groups masked by specific ancestry (mask_group). 
<mark>Beware to modify variable "**data_dir**" to your own data directory, which includes a set of chromosome data directories, each with vcf and msp files, in `gnomix_frq.py`.
Please also provide a sample information file, e.g.,`hawaiiPopInfo.csv`, that gives sample's unique identifier created by the family ID and sample ID, "**famid_id**", and the sample's population group, "**population**", at the same path of `gnomix_frq.py`.</mark>

<table frame=void><tr><td bgcolor=F2F2F2><font face="monospace" size=4>
usage: gnomix_frq.py [-h] [-m MASK_GROUP]

optional arguments:
  -h, --help            show this help message and exit
  -m MASK_GROUP, --mask_group MASK_GROUP
                        mask subpopulation group for alleles
</font></td></tr></table>


An example of the execution code:
```
python3 gnomix_frq.py -m Polynesian
```

## Second Stage
`main.py`: a main execution file that computes genetic statistics.
Currently, ancestry-specific versions of "**F2**", "**F3**", "**F4**", "**F_st**", "**Pi**", "**Psi**", and "**Heterozygosity**" are supported.


<table frame=void><tr><td bgcolor=F2F2F2><font face="monospace" size=4>
usage: main.py [-h] -f FILE [FILE ...] -t DATA_DIR [-b BLOCKSIZE] [-n NUM_REPLICATES]  
               [-r GROUP [GROUP ...]] [-D DAF] [-d DOWNSAMPLE_SIZE] [--rm_DA_files RM_DA_FILES]  
               {F2,F3,F4,Pi,Psi,F_ST,H} 

positional arguments:  
  {F2,F3,F4,Pi,Psi,F_ST,H}    the statistic to compute (H stands for heterozygosity)

optional arguments:  
  -h, --help            show this help message and exit  
  -f FILE [FILE ...], --file FILE [FILE ...]  
                        (REQUIRED) files containing a list of population names for computation,  
                        each file has one column of names, only one file is needed for statistics that are computed  
                        on all pairs of populations (complete graph), two files are needed for statistics that  
                        are computed between only particular pairs of populations (bipartite graph) [^1]  
  -t DATA_DIR, --data_dir DATA_DIR  
                        (REQUIRED) directory containing the snp frequency file [^2]  
  -b BLOCKSIZE, --blocksize BLOCKSIZE  
                        block size for block bootstrap (number of consecutive snps), default=500  
  -n NUM_REPLICATES, --num_replicates NUM_REPLICATES  
                        number of replicates for block bootstrap, default=100  
  -r REF_GROUP [REF_GROUP ...], --ref_group REF_GROUP [REF_GROUP ...]  
                        (F3/Psi only) reference population group for F3/Psi statistics, if  
                        multiple population group names, create a single aggregated reference out of them  
  -D DAF, --DAF DAF     (Psi only) derived allele frequency threshold, (a float number between  
                        [0.0, 0.5]) [^3] , default=-1.0 (i.e. indicates non-Psi stats are being computed)  
  -d DOWNSAMPLE_SIZE, --downsample_size DOWNSAMPLE_SIZE  
                        (Psi only) downsampling size for computing Psi statistics [^4], default=2
  -m MASK, --mask MASK  Use masked data for genetic statistics computation  
  --rm_DA_files RM_DA_FILES  
                        (Boolean) whether to will remove the intermediary files that indicate SNP positions  
                        for derived alleles at the end of each run, default=True
</font></td></tr></table>

[^1]: All pair example (complete graph):
`[A1, A2, A3] => A1A2, A1A3, A2A3`,
cross-match pair example (bipartite graph):
`[A1, A2], [B1, B2] => A1B1, A1B2, A2B1, A2B2`.
**Heterozygosity** only supports all pair example (1 population list file), whereas **F4** only supports cross-match pair example (3 population list files, i.e., one reference population is specified and fixed). The rest statistics support both format.

[^2]: In the file directory, there should be frequency files named as `<pop>.freq`, for example, `Samoa.freq`, which has the first line denoting the column names of the csv data file (comma-separated format) - **CHR, SNP, A1, A2, MAF, NCHROBS, MAF_MSK, NCHROBS_MSK**:
"**CHROM_IDX**" - chromosome ID
"**SNP**" - SNP physical position
"**A1**" - alternative allele
"**A2**" - reference allele
"**MAF**" - alternative allele frequency
"**NCHROBS**" - total allele observations
"**MAF_MSK**" - alternative allele frequency with ancestry-masked SNP data
"**NCHROBS_MSK**" - total allele observations with ancestry-masked SNP data
  
[^3]: When `DAF=0.5`, the Psi computation will have no polarization.

[^4]: Beware that a valid downsampling size should be less than the minimum number of observations for the population groups passed into the statistics.

An example of the execution code:
```
python3 main.py Psi -f pop_list.txt -t data_American -b 50 -n 200 -r Samoa Tonga -D 0.05 -d 2 -m True --rm_DA_files True
```

<mark>(Psi only)</mark> The program will first generate an aggregated population file if there are multiple reference population groups and the file name will contain the first two letters of each reference group capitalized. So for the example above, you can find a file `SamTon.freq` in the frequency data file directory. Then it will generate a file specifying the derived allele positions, which is named in the following format: `psi_<aggr_name>_frq_<DAF>.txt`.

The outputs of the given statistics are in the directory: `path/to/dir/<pop_gen_stat>_output`. The statistics matrix file, `<pop_gen_stat>_mtx.csv`, contains a matrix of the sample mean of the given statistics, where rows and columns represent `pop_list1.txt` and `pop_list2.txt`, respectively. The statistics output file, `<pop_gen_stat>_stat.txt`, contains a list of rows of outputs:
```
popA-popB, <pop_gen_stat>_mean, <pop_gen_stat>_SE, num_of_SNPs_used
```

NOTE: `<pop_gen_stat>_mtx.csv` is overwritten at each run, while `<pop_gen_stat>_stat.txt` uses appending format for each run.

Cite: Ioannidis, A. G., Blanco-Portillo, J., Sandoval, K., Hagelberg, E., Barberena-Jonas, C., Hill, A. V., ... & Moreno-Estrada, A. (2021). Paths and timings of the peopling of Polynesia inferred from genomic networks. Nature, 597(7877), 522-526.
