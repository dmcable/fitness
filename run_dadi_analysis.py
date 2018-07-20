from load_data import load_global_allele_df
import pandas as pd
from VariantQuartiles import load_variant_quartiles
from numpy import np


def gen_sfs(strat_name, variant_count_list):
    save_file = 'saved_data/sfs_strat_' + strat_name + '.npy'
    max_freq = 10000
    freq_vect = np.zeros(max_freq + 1)
    not_counted = 0
    for allele_count in variant_count_list:
        if(allele_count <= max_freq):
            freq_vect[allele_count] += 1
        else:
            not_counted += 1
    np.save(save_file, freq_vect)


def create_allele_df_stratified(allele_df, strat_names, variant_quartiles):
    variant_lists = {'0': [], '1': [], '2': [], '3': []}
    mutation_consequences = set(['missense_variant'])
    for key, quartile in variant_quartiles.iteritems():
        if(key != '' and key in allele_df.index):
            allele_row = allele_df.loc[key]
            cond1 = len(mutation_consequences.intersection(allele_row['consequence'])) > 0
            cond2 = allele_row['FILTER'] == 'PASS'
            if(cond1 and cond2):
                variant_lists[quartile].append(key)
    for strat_name in strat_names:
        quartile_df = allele_df.loc[variant_lists[strat_name]]
        gen_sfs(strat_name, quartile_df['AC'].tolist())


if __name__ == '__main__':
    load_global_allele_df()
    allele_df = pd.read_pickle('saved_data/allele_df.pkl')
    strat_names = ['0', '1', '2', '3']
    variant_quartiles = load_variant_quartiles()
    create_allele_df_stratified(allele_df, strat_names, variant_quartiles)
