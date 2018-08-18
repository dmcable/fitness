# author: Dylan Cable
# date: June 20, 2018

Forward simulation software packages:
1. Dadi https://bitbucket.org/gutenkunstlab/dadi, http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1000695
2. Fitdadi: https://github.com/LohmuellerLab/fitdadi, https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5419480/pdf/345.pdf

Data File Format:
0. Data is from gnomad (https://console.cloud.google.com/storage/browser/gnomad-public/release/2.0.2/coverage?pli=1), and topmed bravo (got from Sofia on Illumina server)
1. Genomes go in input_data/Genomes folder in VCF format (see doc/Genome_Data_Format image)
2. Coverage goes in input_data/coverage folder (see doc/Coverage_Data_Format image)
3. Mutation Probabilities in 'mutation_probabilities.xls' from Jeremy
4. Primate AI Scores in input_data/primate_scores.txt (get a copy from Kyle)
5. trinucleotide mutaiton rates in 'input_data/rates.txt'

Note: Many times in my code, I will have something like "computed_allele_df = True", which will load the cached object from a file
If you truly want to rerun all of my code (which I recommend), then you should find all these cached points, and set them to False.

Running the steady-state model
1. Just run calculate_s_table.py

Running the forward-simulation model
1. run fitness/run_dadi_analysis.py to generate the sfs in fitdadi-master/sfs. These are the 4 quartile sfs and synon.
NOTE: Be careful (optionally) running the function gen_gene_specific_sfs. Right now, this creates separate sfs files for each 20000 genes; you might want to just run a subset of genes
2. run fitdadi-master/run_sfs_method.py to fit the selection distribution by the forward simulation method. Note that currently as next_N is set to 10000, a lot of downsampling is occurring.

Subfunctions of the steady-state model
1. load_global_allele_df loads the allele dataframe, which contains variants in rows. For each variant, we have the allele_count (AC), the consequences, the location, etc.
2. combine_variants, which is called by load_global_allele_df combines the variants in the allele_df to avoid redundancy because, presumably, there will be redundancy from having aggregated data across multiple sources.
3. load_data calls load_gene_distribution, which calls create_gene_df. create_gene_df creates a gene_df, which for each gene calculates the total number of variants and the total mutation rate. This is done for ptv, missense, synon, etc.
4. VariantQuartiles sorts all the variants for a gene into quartiles based on PrimateAI Score
5. FitnessEstimator fits the prior for the selection coefficients and, given a gene dataframe, estimates a posterior selection coefficient distribution for each gene.
6. Simulate_null_mutations is an attempt to create a forward simulation, but is not really used.
7. simulated_data creates ground truth data to validate the methods of FitnessEstimator
8. evaluate_prior_method evaluates the method of fitting a prior (iterative vs maximum likelihood)
9. analyze_strat aggregates everything into one dataframe and outputs the file 'saved_data/gene_df_strat_means.pkl'. After obtaining this file, I usually convert it into a csv format.

Subfunctions of the forward simulation method (all in the fitdadi-master folder)
(my two main scripts are fit_demog.py and run_sfs_method.py)
1. fit_demog fits a demographic model to the synonomys SFS. The demographic model is in:
2. two_epoch (a function in fit_demog), which is a simple demographic model with two parameters nu and T. For more information, locate the fitdadi manual in section 4 Demographic Inference. fitdadi-master/manual_examples/manual.pdf
3. two_epoch_sel (a function in run_sfs_method) which is a demographic model with an additional parameter for selection.
4. run_sfs_method.py calls fit_demog, then fits a selection coefficient distribution for each quartile SFS.
5. load_specra_and_demog (in run_sfs_method) simulates for a number of selection coefficients, the SFS. The set of these simulated SFS is called the spectra and is computed by Selection.spectra.
6. fit_sel_params fits a selection coefficient distribution to the empirical sfs (data). The result of the model fit is located in 'sel_params_opt', which contains for each quartile the mean and shape parameters of the optimal selection coefficient distribution.
7. fit_gene_sel fits a selection distribution for each gene. This method has not been completely debugged/tested.  
