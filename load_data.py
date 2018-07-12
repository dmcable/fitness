import gzip
import pdb
import pandas as pd
from itertools import compress
import numpy as np
from genePredExt import parse_genePredExt
import coverage
import VariantQuartiles

N_genomes = 60706


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
    allele_properties['FILTER'] = [allele_data[6] for allele_index in range(num_alleles)]
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
    assert(mutation_type in ['mis', 'ptv', 'syn'])
    if(mutation_type == 'mis'):
        mutation_consequences = set(['missense_variant'])

        def rate_function(gene): return 10**gene.mis
    elif(mutation_type == 'syn'):
        mutation_consequences = set(['synonymous_variant'])

        def rate_function(gene): return 10**gene.syn
    else:
        mutation_consequences = set(
            ['stop_gained', 'splice_acceptor_variant', 'splice_donor_variant']
        )

        def rate_function(gene): return combine_log_rates(gene.non, gene.splice_site)
    gene_df['total_mutations'] = N_genomes * mutation_rates.apply(rate_function, axis=1)
    gene_df['allele_count'] = 0
    for index, allele_row in allele_df.iterrows():
        gene_consequences = allele_row['consequence']
        for gene in gene_consequences.keys():
            cond1 = gene != ''
            cond2 = len(mutation_consequences.intersection(gene_consequences[gene])) > 0
            cond3 = gene in gene_df.index
            cond4 = allele_row['FILTER'] == 'PASS'
            if(cond1 and cond2 and cond3 and cond4):
                gene_df.loc[gene, 'allele_count'] += allele_row['AC']
    return gene_df


def create_gene_df_stratified(allele_df, mutation_rates, strat_names, variant_quartiles):
    gene_df = pd.DataFrame(index=mutation_rates.index)
    mutation_consequences = set(['missense_variant'])
    valid_keys = set(variant_quartiles.keys())
    def rate_function(gene): return 10**gene.mis

    gene_df['total_mutations'] = N_genomes * mutation_rates.apply(rate_function, axis=1) / len(strat_names)
    for name in strat_names:
        gene_df['ac_' + name] = 0
        gene_df['vc_' + name] = 0
        gene_df['vu_' + name] = 0
    for index, allele_row in allele_df.iterrows():
        lookup_key = VariantQuartiles.get_lookup_key(
            str(allele_row['POS']), str(allele_row['REF']), str(allele_row['ALT'])
        )
        gene_consequences = allele_row['consequence']
        for gene in gene_consequences.keys():
            cond1 = gene != ''
            cond2 = len(mutation_consequences.intersection(gene_consequences[gene])) > 0
            cond3 = gene in gene_df.index
            cond4 = allele_row['FILTER'] == 'PASS'
            cond5 = allele_row['AC'] < 0.001 * N_genomes
            cond6 = lookup_key in valid_keys
            if(cond1 and cond2 and cond3 and cond4 and cond5 and cond6):
                gene_df.loc[gene, 'ac_' + variant_quartiles[lookup_key]] += allele_row['AC']
                gene_df.loc[gene, 'vc_' + variant_quartiles[lookup_key]] += 1
                if(allele_row['AC'] >= 1):
                    gene_df.loc[gene, 'vu_' + variant_quartiles[lookup_key]] += 1
    return gene_df


def load_allele_df(vcf_dataset_path, gene_list):
    MAX_ALLELE_COUNT = 6000000000
    n_lines = 0
    passed_header = False
    col_names = ['POS', 'REF', 'ALT', 'AC', 'FILTER', 'consequence']
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
    vcf_dataset_path = 'input_data/ExAC.vcf.gz'
    mut_rates_path = 'input_data/mutation_probabilities.xls'
    mutation_rates = pd.read_excel(mut_rates_path, index_col=1, sheet_name=1)
    allele_df = load_allele_df(vcf_dataset_path, mutation_rates.index)
    gene_df_missense = create_gene_df(allele_df, mutation_rates, 'mis')
    gene_df_ptv = create_gene_df(allele_df, mutation_rates, 'ptv')
    if(save_data):
        allele_df.to_pickle(data_file_names['allele'])
        gene_df_missense.to_pickle(data_file_names['missense'])
        gene_df_ptv.to_pickle(data_file_names['ptv'])
    return allele_df, gene_df_missense, gene_df_ptv


def load_syn():
    mut_rates_path = 'input_data/mutation_probabilities.xls'
    mutation_rates = pd.read_excel(mut_rates_path, index_col=1, sheet_name=1)
    allele_df = pd.read_pickle('saved_data/allele_df.pkl')
    gene_df_syn = create_gene_df(allele_df, mutation_rates, 'syn')
    gene_df_syn.to_pickle('saved_data/gene_df_syn.pkl')
    gene_df_syn = filter_gene_df(gene_df_syn, 'saved_data/gene_df_filtered_syn.pkl')
    pdb.set_trace()


def load_gene_distribution_saved(data_file_names):
    gene_df_missense = pd.read_pickle(data_file_names['missense'])
    gene_df_ptv = pd.read_pickle(data_file_names['ptv'])
    return gene_df_missense, gene_df_ptv


def get_coverage_filter(gene_index):
    gene_filter = pd.Series(index=gene_index)
    gene_filter.loc[:] = False  # default is to drop the gene
    line_count = 0
    gene_coordinate_file = 'input_data/canonical_gencode_gene_structure.txt'
    coverage_file = 'input_data/coverage'
    with open(gene_coordinate_file, 'r') as coordinate_reader:
        for curr_line in coordinate_reader:
            if(line_count % 100 == 0):
                print(line_count)
            line_count += 1
            gene_coordinates = parse_genePredExt(curr_line)
            cov_dict = coverage.open_coverage(coverage_file)
            chrom = getattr(gene_coordinates, 'chrom')
            exons = getattr(gene_coordinates, 'exons')
            gene = getattr(gene_coordinates, 'name2')
            successes = 0
            failures = 0
            if(gene in gene_filter):
                for i in range(len(exons)):
                        successes_to_add, failures_to_add = coverage.good_coverage(
                            cov_dict[chrom], chrom, exons[i][0], exons[i][1]
                        )
                        successes += successes_to_add
                        failures += failures_to_add
                gene_filter[gene] = successes >= failures
    return gene_filter


def filter_gene_df(gene_df, file_name):
    gene_filter = gene_df['allele_count'] <= 0.01 * N_genomes
    gene_df = gene_df[gene_filter]
    gene_filter = get_coverage_filter(gene_df.index)
    gene_df = gene_df[gene_filter]
    gene_df.to_pickle(file_name)
    return gene_df


def load_data():
    cached_data = True
    filtered = True
    data_file_names = {
        'allele': 'saved_data/allele_df.pkl', 'missense': 'saved_data/gene_df_missense.pkl',
        'ptv': 'saved_data/gene_df_ptv.pkl'
    }
    filtered_file_names = {
        'missense': 'saved_data/gene_df_filtered_missense.pkl',
        'ptv': 'saved_data/gene_df_filtered_ptv.pkl'
    }
    if(filtered):
        for key, value in filtered_file_names.iteritems():
            data_file_names[key] = value
    if(not cached_data):
        save_data = True
        allele_df, gene_df_missense, gene_df_ptv = load_gene_distribution(
            save_data, data_file_names
        )
    else:
        gene_df_missense, gene_df_ptv = load_gene_distribution_saved(data_file_names)
    if(not filtered):
        gene_df_ptv = filter_gene_df(gene_df_ptv, filtered_file_names['ptv'])
        gene_df_missense = filter_gene_df(gene_df_missense, filtered_file_names['missense'])
    return gene_df_ptv, gene_df_missense
