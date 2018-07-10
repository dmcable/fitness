import pandas as pd
import pdb
import numpy as np


def load_gene_scores():
    path = 'input_data/primate_scores.txt'
    gene_scores = pd.DataFrame(columns=['gene', 'avg_score'])
    curr_gene = ''
    gene_count = 0.0
    total_score = 0.0
    first_line = True
    with open(path, 'r') as f:
        for line in f:
            if(not first_line):
                elems = line.split('\t')
                if(curr_gene != elems[4]):
                    if(curr_gene != ''):
                        series = pd.Series({'gene':curr_gene, 'avg_score':total_score / gene_count})
                        gene_scores = gene_scores.append(series,ignore_index=True)
                    curr_gene = elems[4]
                    gene_count = 0.0
                    total_score = 0.0

                total_score += float(elems[5][:-1])
                gene_count += 1
            else:
                first_line = False
    series = pd.Series({'gene':curr_gene, 'avg_score':total_score / gene_count})
    gene_scores = gene_scores.append(series,ignore_index=True)
    gene_scores = gene_scores.set_index(gene_scores['gene'], drop=True)
    gene_scores = gene_scores.sort_values('avg_score')
    return gene_scores


if __name__ == '__main__':
    gene_scores_loaded = True
    if(gene_scores_loaded):
        gene_scores = load_gene_scores()
        gene_scores.to_pickle('saved_data/gene_scores.pkl')
    else:
        gene_scores = pd.read_pickle('saved_data/gene_scores.pkl')

    gene_df = pd.DataFrame(columns=['primate_score', 's_mis', 's_ptv', 's_syn'])
    index_names = set(gene_scores.index.tolist())
    mis_df = pd.read_pickle('saved_data/gene_df_missense_means.pkl')
    ptv_df = pd.read_pickle('saved_data/gene_df_ptv_means.pkl')
    syn_df = pd.read_pickle('saved_data/gene_df_syn_means.pkl')
    for gene in index_names:
        primate_score = np.mean(gene_scores.loc[gene]['avg_score'].tolist())
        info = {'primate_score': primate_score}
        if(gene in mis_df.index):
            info['s_mis'] = mis_df.loc[gene]['s_mean']
        if(gene in ptv_df.index):
            info['s_ptv'] = ptv_df.loc[gene]['s_mean']
        if(gene in syn_df.index):
            info['s_syn'] = syn_df.loc[gene]['s_mean']
        new_pt = pd.DataFrame(data=[info], index = [gene])
        gene_df = gene_df.append(new_pt)
    gene_df.to_pickle('saved_data/primate_s_het_table.pkl')
    pdb.set_trace()

    gene_scores = gene_scores.set_index(gene_scores['gene'], drop=True)
    #gene_scores = gene_scores.loc[missense_df.index]
    gene_scores = gene_scores.sort_values('avg_score')
    N = len(gene_scores)
    ranges = [0, 240, 480, 720, 955]
    means = [0, 0, 0, 0]
    stds = [0, 0, 0, 0]
    for index, gene in gene_scores.iterrows():
        if(index in missense_df.index):
            gene_scores.loc[index]['s_missense'] = missense_df.loc[index]['s_mean']
            pdb.set_trace()
    pdb.set_trace()
    for i in range(4):
        num_genes = ranges[i+1] - ranges[i]
        results = np.zeros(num_genes)
        for j in range(ranges[i], ranges[i+1]):
            gene = gene_scores.iloc[j]['gene']
            results[j - ranges[i]] = missense_df.loc[gene]['s_mean']
        means[i] = np.mean(results)
        stds[i] = np.std(results)
