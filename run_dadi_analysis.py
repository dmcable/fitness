from load_data import load_global_allele_df
import pandas as pd
from VariantQuartiles import load_variant_quartiles
import numpy as np
import pdb
import time

def gen_sfs(strat_name, variant_count_list):
    save_file = 'saved_data/sfs_full_strat_' + strat_name + '.npy'
    max_freq = 10000
    freq_vect = np.zeros(max_freq + 1)
    not_counted = 0
    n = 0
    for allele_count in variant_count_list:
        n += 1
        if(n % 50000 == 0):
            print(n)
        if(allele_count <= max_freq):
            freq_vect[allele_count] += 1
        else:
            not_counted += 1
    np.save(save_file, freq_vect)


def create_allele_df_stratified(allele_df, strat_names, variant_quartiles):
    variant_lists = {'0': [], '1': [], '2': [], '3': []}
    mutation_consequences = set(['missense_variant'])
    n = 0
    for key, quartile_row in variant_quartiles.iterrows():
        n += 1
        if(n % 50000 == 0):
            print(n)
        quartile = quartile_row['quartile']
        if(key != '' and key in allele_df.index):
            allele_row = allele_df.loc[key]
            consequences = allele_row['consequence']
            list_consequences = consequences[consequences.keys()[0]]
            cond1 = len(mutation_consequences.intersection(list_consequences)) > 0
            cond2 = allele_row['FILTER'] == 'PASS'
            if(cond1 and cond2):
                variant_lists[quartile].append(key)
    for strat_name in strat_names:
        print(strat_name)
        quartile_df = allele_df.loc[variant_lists[strat_name]]
        gen_sfs(strat_name, quartile_df['AC'].tolist())


def create_allele_df_synon(allele_df):
    variant_lists = []
    mutation_consequences = set(['synonymous_variant'])
    n = 0
    for key, allele_row in allele_df.iterrows():
        n += 1
        if(n % 50000 == 0):
            print(n)
        if(key != '' and key in allele_df.index):
            consequences = allele_row['consequence']
            list_consequences = consequences[consequences.keys()[0]]
            cond1 = len(mutation_consequences.intersection(list_consequences)) > 0
            cond2 = allele_row['FILTER'] == 'PASS'
            if(cond1 and cond2):
                variant_lists.append(key)
    quartile_df = allele_df.loc[variant_lists]
    gen_sfs('syn', quartile_df['AC'].tolist())


def main_fun(allele_df=None):
    loaded_sfs = False
    if(not loaded_sfs):
        load_global_allele_df()
        if(allele_df is None):
            allele_df = pd.read_pickle('saved_data/allele_df.pkl')
        strat_names = ['0', '1', '2', '3']
        variant_quartiles = load_variant_quartiles()
        create_allele_df_stratified(allele_df, strat_names, variant_quartiles)
        create_allele_df_synon(allele_df)
        time.sleep(10000000)


if __name__ == '__main__':
    main_fun()
