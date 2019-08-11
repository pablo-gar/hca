import scanpy as sc
import os
import model_preprocessing
import sys
import pandas as pd
import json

DPI = 300

def main():


    h5ad_file = sys.argv[1]
    cellType_label = sys.argv[2]
    out_prefix = sys.argv[-1]

    #h5ad_file = '/Users/pgarcia-nietochanzuckerberg.com/projects/cell_type_transfer/pancreas/data/hca_model_alternative.h5ad'
    #out_prefix = '/Users/pgarcia-nietochanzuckerberg.com/projects/cell_type_transfer/pancreas/results/modeling_alternative_'
    #cellType_label = 'cellType'

    # Setting plotting settings
    sc.settings.autoshow = False
    sc.settings.figdir = os.path.dirname(out_prefix)

    # Read data
    annData = sc.read(h5ad_file)

    # Process data
    annData = model_preprocessing.preprocessing(annData, do_log1p=True, do_select_variable_genes=True, do_min_cells_filter=True,
                                                do_min_genes_filter=True, percent_cells=0.05, min_genes=10, n_var_genes=500)

    sc.pp.normalize_total(annData, exclude_highly_expressed=True)

    # Get available models
    models = [i for i in annData.obs.keys().to_list() if i.startswith(cellType_label + '_')]
    for m in models:
        correct_label = annData.obs[cellType_label] == annData.obs[m].to_list()
        annData.obs[m + '_accuracy'] = ['correct' if i else 'mislabelled' for i in correct_label]


    # Get data ready for plotting
    sc.pp.pca(annData)
    sc.pp.neighbors(annData)
    sc.tl.umap(annData)

    sc.pl.umap(annData, color=cellType_label)
    sc.pl.umap(annData, color=cellType_label, save=os.path.basename(out_prefix) + '_' + cellType_label + '_original_umap.pdf')

    with open(out_prefix + 'performance_metrics_per_model_test.tsv', 'w') as f:
        print('model', 'correct', 'mislablled', 'accuracy', sep="\t", file=f)
        for m in models:
            sc.pl.umap(annData, color=m, save=os.path.basename(out_prefix) + '_' + m + '_predicted_umap.pdf')
            sc.pl.umap(annData, color=m + '_accuracy', save=os.path.basename(out_prefix) + '_' + m + '_missed_umap.pdf')
            correct = sum(annData.obs[m + '_accuracy'] == 'correct')
            miss = sum(annData.obs[m + '_accuracy'] != 'correct')
            print(m, correct, miss, correct / (correct + miss), sep ="\t", file=f)

    # Compile cv results
    cv = {}
    cv_keys = [i for i in annData.uns.keys() if i.startswith('metrics_')]
    for a in cv_keys:
        cv[a] = pd.DataFrame(annData.uns[a])
        cv[a]['model'] = a.split("metrics_")[1]

    # Writing
    cv = pd.concat(cv)
    cv.to_csv(out_prefix + 'performance_metrics_per_model_train.tsv', sep='\t', index=False)

    with open(out_prefix + 'best_pars_per_model.md', 'w') as f:
        best_pars_keys = [i for i in annData.uns.keys() if i.startswith('best_pars_')]
        for b in best_pars_keys:
            print('* ', b.split('best_pars_')[1], sep = '', file=f)
            for key in annData.uns[b]:
                print('  * ', key, ': ', annData.uns[b][key], file=f)
