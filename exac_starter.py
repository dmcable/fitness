import gzip
import pdb
import pandas as pd
from itertools import compress
import numpy as np
from FitnessEstimatorNewer import FitnessEstimator


def load_consequences(transcripts, allele_properties, gene_list):
    for transcript in transcripts:
        transcript_properties = transcript.split('|')
        allele, consequences_str, _, gene = transcript_properties[0:4]
        # exclude indels and if this variant is not found
        cond1 = gene in gene_list
        cond2 = len(allele) == len(allele_properties['REF'][0])
        cond3 = allele in allele_properties['ALT']
        if(cond1 and cond2 and cond3):
            consequences_set = set(consequences_str.split('&'))
            allele_index = allele_properties['ALT'].index(allele)
            consequence_dict = allele_properties['consequence'][allele_index]
            if(gene not in consequence_dict.keys()):
                consequence_dict[gene] = set()
            consequence_dict[gene] = consequence_dict[gene].union(consequences_set)


def get_allele_properties(line, gene_list):
    allele_properties = {}
    allele_data = line.split('\t')
    allele_properties['ALT'] = allele_data[4].split(',')
    num_alleles = len(allele_properties['ALT'])
    allele_properties['POS'] = [int(allele_data[1]) for allele_index in range(num_alleles)]
    allele_properties['REF'] = [allele_data[3] for allele_index in range(num_alleles)]
    allele_info = allele_data[7].split(';')
    allele_properties['consequence'] = [{} for allele_index in range(num_alleles)]
    for i in range(len(allele_info)):
        entry = allele_info[i].split('=')
        if(entry[0] == 'AC'):
            entry[1] = entry[1].split(',')
            allele_properties['AC'] = [int(ac) for ac in entry[1][0:num_alleles]]
        elif(entry[0] == 'CSQ'):
            transcripts = entry[1].split(',')
            load_consequences(transcripts, allele_properties, gene_list)
    consequences = allele_properties['consequence']
    filter = [len(consequences[allele_index]) > 0 for allele_index in range(num_alleles)]
    return allele_properties, filter


def combine_log_rates(rate1, rate2):
    total_rate = 0
    if(not np.isnan(rate1)):
        total_rate += 10**rate1
    if(not np.isnan(rate2)):
        total_rate += 10**rate2
    return total_rate


def create_gene_df(allele_df, mutation_rates, mutation_type):
    gene_df = pd.DataFrame(index=mutation_rates.index)
    assert(mutation_type in ['mis', 'ptv'])
    if(mutation_type == 'mis'):
        mutation_consequences = set(['missense_variant'])
        gene_df['mut_rate'] = mutation_rates.apply(lambda gene: 10**gene.mis, axis=1)
    else:
        mutation_consequences = set(
            ['stop_gained', 'splice_acceptor_variant', 'splice_donor_variant']
        )
        gene_df['mut_rate'] = mutation_rates.apply(
            lambda gene: combine_log_rates(gene.non, gene.splice_site), axis=1
        )
    gene_df['allele_count'] = 0
    for index, allele_row in allele_df.iterrows():
        gene_consequences = allele_row['consequence']
        for gene in gene_consequences.keys():
            cond1 = gene != ''
            cond2 = len(mutation_consequences.intersection(gene_consequences[gene])) > 0
            cond3 = gene in gene_df.index
            if(cond1 and cond2 and cond3):
                gene_df.loc[gene, 'allele_count'] += allele_row['AC']
    return gene_df


def load_allele_df(vcf_dataset_path, gene_list):
    MAX_ALLELE_COUNT = 6000000000
    n_lines = 0
    passed_header = False
    col_names = ['POS', 'REF', 'ALT', 'AC', 'consequence']
    allele_df_lists = {}
    for name in col_names:
        allele_df_lists[name] = []
    with gzip.open(vcf_dataset_path, 'r') as vcf_reader:
        for curr_line in vcf_reader:
            n_lines += 1
            if(passed_header):
                allele_properties, filter = get_allele_properties(curr_line, gene_list)
                for name, curr_list in allele_df_lists.iteritems():
                    curr_list.extend(list(compress(allele_properties[name], filter)))
            else:
                if(curr_line[0:6] == '#CHROM'):
                    passed_header = True
            if(len(allele_df_lists['POS']) >= MAX_ALLELE_COUNT):
                break
            if(n_lines % 10000 == 0):
                print(str(n_lines) + ' ' + str(len(allele_df_lists['POS'])))
    allele_df = pd.DataFrame(allele_df_lists)
    return allele_df


def load_gene_distribution(save_data, data_file_names):
    mutation_rates = pd.read_excel('mutation_probabilities.xls', index_col=1, sheet_name=1)
    vcf_dataset_path = 'ExAC.vcf.gz'
    allele_df = load_allele_df(vcf_dataset_path, mutation_rates.index)
    allele_df.to_pickle(data_file_names['allele'])
    pdb.set_trace()
    gene_df_missense = create_gene_df(allele_df, mutation_rates, 'mis')
    gene_df_ptv = create_gene_df(allele_df, mutation_rates, 'ptv')
    if(save_data):
        allele_df.to_pickle(data_file_names['allele'])
        gene_df_missense.to_pickle(data_file_names['missense'])
        gene_df_ptv.to_pickle(data_file_names['ptv'])
    return allele_df, gene_df_missense, gene_df_ptv


def load_gene_distribution_saved(data_file_names):
    gene_df_missense = pd.read_pickle(data_file_names['missense'])
    gene_df_ptv = pd.read_pickle(data_file_names['ptv'])
    return gene_df_missense, gene_df_ptv


if __name__ == '__main__':
    cached_data = True
    data_file_names = {
        'allele': 'saved_data/allele_df.pkl', 'missense': 'saved_data/gene_df_missnese.pkl',
        'ptv': 'saved_data/gene_df_ptv.pkl'
    }
    if(not cached_data):
        save_data = True
        allele_df, gene_df_missense, gene_df_ptv = load_gene_distribution(
            save_data, data_file_names
        )
    else:
        gene_df_missense, gene_df_ptv = load_gene_distribution_saved(data_file_names)
    N_genomes = 60706
    estimator_ptv = FitnessEstimator(N_genomes, gene_df_ptv)

    pdb.set_trace()
