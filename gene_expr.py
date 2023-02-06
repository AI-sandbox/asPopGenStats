import os
import numpy as np
import sys
import pandas as pd
import argparse
from collections import defaultdict
import errno
from tqdm import tqdm

# TODO: change the Gnomix path accordingly and dictionary "race_info" based 
# TODO: on Gnomix configuration.
# dest_data_path = '/home/markpenj/hawaii/gnomix_run/output_files/'
race_info = {'American':0, 'Papuan':1, 'Polynesian':2, 'European':3,
             'EastAsian':4, 'African':5}
dest_data_path = '/home/projects/Sandwich_Isles/gnomix_updatedTargets/newRun/'
# race_info = {'American':0, 'Papuan':2, 'AsiaPacific':1, 'European':3,
#              'African':4}

argp = argparse.ArgumentParser()
argp.add_argument('alt_group', help="alternative population group for alleles",
                  choices=race_info.keys())
argp.add_argument('ref_group', help="reference population group for alleles",
                  choices=race_info.keys())
args = argp.parse_args()

# Create a dictionary "pop_info_dict" with "famid_id: pop"
pop_info = pd.read_csv('updated_popinfo.csv', header=0, usecols=['id', 'famid', 'pop'])
pop_info['famid_id'] = pop_info['famid'] + '_' + pop_info['id']
pop_info = (pop_info.drop(['id', 'famid'], axis=1))

# TODO: Remove ids in your target population, e.g., Hawaii. 
# ! We do this step because we receive a new list of ids in the target
# ! population, see line 49: a subset of Hawaii ids with 5% East Asian component.
# ! We will add the filtered ids to "pop_info_dict"
# pop_info = pop_info[pop_info['pop'] != 'RapaNui']
pop_info = pop_info[pop_info['pop'] != 'Hawaii']
print(pop_info.head(5))

pop_info_dict = {data.values[1]: data.values[0] for _, data in pop_info.iterrows()}
# print(list(pop_info_dict.keys())[:5])
# print(pop_info_dict['0_HG02687'])

# ! Add newly filtered Hawaii ids to "pop_info".
# rapanui_info = pd.read_csv('subset_RapaNui_-0.10EUR.ids', header=0,
#                             usecols=['id', 'famid', 'pop'], delimiter=r'\s+',
#                             dtype=str)
# rapanui_info['famid_id'] = rapanui_info['famid'] + '_' + rapanui_info['id']
# rapanui_info = list(rapanui_info['famid_id'])
# for rapanui_id in rapanui_info:
#     pop_info_dict[rapanui_id] = 'RapaNui'
hawaii_info = pd.read_csv('subset_Hawaii_-0.05EAS.ids', header=0,
                            usecols=['id', 'pop'], delimiter=r'\s+',
                            dtype=str)
hawaii_info = list(hawaii_info['id'])
for hawaii_id in hawaii_info:
    pop_info_dict[hawaii_id] = 'Hawaii'
# print([k for k, v in pop_info_dict.items() if v == 'RapaNui'])

def get_chrom_data(idx, pop_info_dict):
    sample_file = pd.read_csv(dest_data_path + os.path.join(f'chr{idx}', 'query_results.msp'),
                              header=0, skiprows=[0], sep='\t')
    # The first 6 column indices are
    # [chm, spos, epos, sgpos, egpos, n snps]
    data_sidx = 6 # starting column index for samples, e.g., 1_4781017.0
    sample_idx = sample_file.columns[data_sidx:]
    # Get a dictionary "group2colidx" with "population_group: col_index_of_sample".
    group2colidx = defaultdict(list)
    for count, sample in enumerate(sample_idx):
        sample_id = sample[:-2] # remove '.1' or '.0'
        if sample_id in pop_info_dict:
            # sample_id from 'chr' not show up in 'pop_info_dict' b/c filtering
            group2colidx[pop_info_dict[sample_id]].append(count)

    # Filter out the data with irrelevant populations.
    data = sample_file.iloc[:, data_sidx:].to_numpy(dtype=int)
    # unique, counts = np.unique(data[:, 0], return_counts=True)
    # print(np.asarray((unique, counts)).T)
    data_alt = np.isclose(data, race_info[args.alt_group]).astype(float)
    data = data_alt + np.isclose(data, race_info[args.ref_group]).astype(float)
    # unique, counts = np.unique(data[:, 0], return_counts=True)
    # print(np.asarray((unique, counts)).T)

    # Aggregate samples from one population group.
    aggr_data_alt = pd.DataFrame({k: np.divide(np.sum(data_alt[:, v], axis=1),
                                               np.sum(data[:, v], axis=1),
                                               out=np.zeros(data_alt.shape[0]),
                                               where=(np.sum(data[:, v], axis=1) != 0))
                                  for k, v in group2colidx.items()})
    assert (aggr_data_alt <= 1.0).all(axis=None) == True
    aggr_data_tot =  pd.DataFrame({f'{k}_CT': np.sum(data[:, v], axis=1) for k, v in group2colidx.items()})
    aggr_data = pd.concat([aggr_data_alt, aggr_data_tot], axis=1)
    for group in group2colidx:
        assert len(group2colidx[group]) % 2 == 0
    aggr_data['CHROM_IDX'] = [f'{idx}.{i}' for i in range(aggr_data.shape[0])]
    return list(group2colidx.keys()), aggr_data

# Concatenate data from all chromosomes 1-22.
chrom_data_list = []
for i in tqdm(range(1, 23)):
    groups, chrom_data_i = get_chrom_data(i, pop_info_dict)
    chrom_data_list.append(chrom_data_i)
chrom_data = pd.concat(chrom_data_list)

# Save data to the given path.
save_file_path = os.path.join(sys.path[0], f'data_{args.alt_group}_{args.ref_group}')
try:
    os.makedirs(save_file_path)
except OSError as e:
    if e.errno != errno.EEXIST:
        raise

for group in groups:
    df_group = chrom_data[['CHROM_IDX', group, f'{group}_CT']]
    df_group = df_group.rename(columns={group:'FREQ', f'{group}_CT':'CT'})
    df_group.to_csv(f'data_{args.alt_group}_{args.ref_group}/{group:s}.freq', index=False)

np.savetxt(f'pop_groups_{args.alt_group}_{args.ref_group}.txt', groups, fmt='%s')