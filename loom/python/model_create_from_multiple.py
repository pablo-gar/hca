# Merges annDatas using combat

import model_utils
import model_preprocessing
import model_create
import scanpy as sc
import sys
import os

def main():

    model_name = sys.argv[1]
    train_y_name = sys.argv[2]
    conc_type = sys.argv[3] #[combat|concatenate]
    do_log1p = model_preprocessing.str_to_bool(sys.argv[4])
    do_select_variable_genes = model_preprocessing.str_to_bool(sys.argv[5])
    out_train_file = sys.argv[6]
    out_test_file = sys.argv[7]
    out_plots_dir = sys.argv[8]
    test_file = sys.argv[9]
    train_files = sys.argv[10:]

    #model_name = "RandomForest"
    #train_y_name = "cellType"
    #conc_type = "concatenate"
    #do_log1p= True
    #do_select_variable_genes = True
    #out_plots_dir = "/Users/pgarcia-nietochanzuckerberg.com/"
    #test_file = "/Users/pgarcia-nietochanzuckerberg.com/projects/cell_type_transfer/pancreas/data/hca_cellTypes.loom"
    #train_files = ["/Users/pgarcia-nietochanzuckerberg.com/projects/cell_type_transfer/pancreas/data/paper_alternative.loom",
    #               "/Users/pgarcia-nietochanzuckerberg.com/projects/cell_type_transfer/pancreas/data/paper_alternative_GSE86469.loom",
    #               "/Users/pgarcia-nietochanzuckerberg.com/projects/cell_type_transfer/pancreas/data/paper_alternative_GSE84133.loom"
    #               ]


    # Read data
    annDatas = [sc.read(i) for i in train_files + [test_file]]

    for i in annDatas:
        i.var_names_make_unique()

    # Concatenating
    if conc_type == 'concatenate':
        conc = annDatas[0].concatenate(*annDatas[1:])
    elif conc_type == 'combat':
        conc = model_utils.merge_combat(*annDatas)
    else:
        raise Exception('Concatenation type can only be "combat" or "concatenate"')

    # Plotting UMAPs of concatenated data
    # if not os.path.exists(out_plots_dir):
    #    os.makedirs(out_plots_dir)
    # plot_concatenated(conc_basic, out_prefix=os.path.join(out_plots_dir, "concatenated_"))
    # plot_concatenated(conc, out_prefix=os.path.join(out_plots_dir, "combat_"))

    train_ind = conc.obs['batch'].isin(map(str, range(len(train_files))))
    test_ind = conc.obs['batch'] == str(len(train_files))

    train_h5ad = conc[train_ind.to_list(),:]
    test_h5ad = conc[test_ind.to_list(),:]

    train_h5ad_original = train_h5ad.copy()
    test_h5ad_original = test_h5ad.copy()

    train_h5ad_original.uns = annDatas[0].uns
    test_h5ad_original.uns = annDatas[-1].uns

    # Pre-preprocess
    train_h5ad = model_preprocessing.preprocessing(train_h5ad, do_log1p=do_log1p, do_min_cells_filter=True,
                                                   do_select_variable_genes=do_select_variable_genes, flavor='cell_ranger',
                                                   do_min_genes_filter=False)

    test_h5ad = model_preprocessing.preprocessing(test_h5ad, do_log1p=do_log1p, do_min_cells_filter=False, do_select_variable_genes=False, do_min_genes_filter=False)

    #Select same genes
    train_h5ad, test_h5ad = model_utils.subset_annData(train_h5ad, test_h5ad)

    # Get matrix and labels
    if train_y_name in test_h5ad.obs.keys():
        x_train, y_train, x_test, y_test = model_utils.get_mat_and_lab(x=train_h5ad, x_lab=train_y_name, y=test_h5ad,
                                                                       y_lab=train_y_name)
    else:
        x_train, y_train, x_test = model_utils.get_mat_and_lab(x=train_h5ad, x_lab=train_y_name, y=test_h5ad)

    # Train model
    model = model_create.get_model(model_name, x_train, y_train, n_iter=50)


    # Predict on test data and append predictions to original matrix
    y_test_predicted = model.predict(x_test)
    y_train_predicted = model.predict(x_train)

    test_h5ad_original.obs[train_y_name + "_predicted_" + model_name] = y_test_predicted
    train_h5ad_original.obs[train_y_name + "_predicted_" + model_name] = y_train_predicted

    # Add assessment metric for the model as unstructured metadata
    test_h5ad_original.uns['metrics_' + model_name] = model.cv_results_
    train_h5ad_original.uns['metrics_' + model_name] = model.cv_results_

    test_h5ad_original.uns['best_pars_' + model_name] = model.best_params_
    train_h5ad_original.uns['best_pars_' + model_name] = model.best_params_

    # Save data
    test_h5ad_original.write(out_test_file)
    train_h5ad_original.write(out_train_file)


def plot_concatenated(y, out_prefix, batch_key='batch'):

    x = y.copy()


    sc.pp.filter_genes(x, min_cells=int(x.n_obs * 0.05))
    sc.pp.filter_cells(x, min_genes=5)
    sc.pp.highly_variable_genes(x, flavor="cell_ranger")

    sc.pp.pca(x)
    sc.pp.neighbors(x)
    sc.tl.umap(x)

    sc.pl.umap(x, color=batch_key, save='study.pdf')
    sc.pl.umap(x, color='cellType', save='cellTypes.pdf')



if __name__ == '__main__':
    main()
