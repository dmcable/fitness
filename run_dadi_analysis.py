from load_data import load_global_allele_df
import pandas as pd
from VariantQuartiles import load_variant_quartiles
import numpy as np
import pdb
import pickle


def aggregate_ac_counts(allele_df, variant_quartiles):
    aggregated = True
    if(aggregated):
        ac_dict = pickle.load(open('ac_dict.pkl', "rb"))
        return ac_dict
    ac_dict = {}
    mutation_consequences = set(['missense_variant'])
    n = 0
    for key, quartile_row in variant_quartiles.iterrows():
        n += 1
        if(n % 50000 == 0):
            print(n)
        quartile = int(quartile_row['quartile'])
        if(key != '' and key in allele_df.index):
            allele_row = allele_df.loc[key]
            consequences = allele_row['consequence']
            if(allele_row['FILTER'] == 'PASS'):
                for gene, cons in consequences.iteritems():
                    if(mutation_consequences.intersection(cons) > 0):
                        if(gene not in ac_dict.keys()):
                            ac_dict[gene] = [[] for i in range(4)]
                        ac_dict[gene][quartile].append(allele_row['AC'])
    pickle.dump(ac_dict, open('ac_dict.pkl', 'wb'))
    return ac_dict


def count_variants(ac_dict):
    variant_counts = {}
    for gene in ac_dict.keys():
        variant_counts[gene] = sum([len(ac_dict[gene][i]) for i in range(4)])
    return variant_counts


def gen_gene_specific_sfs(allele_df, variant_quartiles):
    ac_dict = aggregate_ac_counts(allele_df, variant_quartiles)
    variant_counts = count_variants(ac_dict)
    pickle.dump(variant_counts, open('variant_counts.pkl', 'wb'))
    #threshold = 500  # determined this somehow by looking at variant_counts
    for gene in variant_counts.keys():
        if(True):  # (variant_counts[gene] > threshold):
            for strat in range(4):
                gen_sfs(str(strat), ac_dict[gene][strat], gene)
    pdb.set_trace()


def gen_sfs(strat_name, variant_count_list, prefix='fullest'):
    save_file = 'gene_sfs/sfs_' + prefix + '_strat_' + strat_name + '.npy'
    max_freq = 500000
    freq_vect = np.zeros(max_freq + 1)
    for allele_count in variant_count_list:
        if(allele_count <= max_freq):
            freq_vect[allele_count] += 1
        else:
            freq_vect[-1] += 1
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


def main_fun(allele_df=None, variant_quartiles=None):
    loaded_sfs = False
    if(not loaded_sfs):
        load_global_allele_df()
        if(allele_df is None):
            allele_df = pd.read_pickle('saved_data/allele_df_comb.pkl')
        strat_names = ['0', '1', '2', '3']
        if(variant_quartiles is None):
            variant_quartiles = load_variant_quartiles()
        create_allele_df_stratified(allele_df, strat_names, variant_quartiles)
        create_allele_df_synon(allele_df)


if __name__ == '__main__':
    main_fun()


def combine_sfs():
    for strat_name in ['0', '1', '2', '3', 'syn']:
        main_sfs = np.load('saved_data/sfs_fuller_strat_' + strat_name + '.npy')
        add_sfs = np.load('saved_data/sfs_full_strat_' + strat_name + '.npy')
        for i in range(len(add_sfs)):
            main_sfs[i] = add_sfs[i]
        save_file = 'saved_data/sfs_fullest_strat_' + strat_name + '.npy'
        np.save(save_file, main_sfs)
