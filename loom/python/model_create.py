# Creates a model for cell type predictions using a h5ad matrix and add predictions to a another h5ad matrix

import sys
import scanpy as sc
from sklearn import svm
from sklearn import ensemble
from sklearn import linear_model
from sklearn import model_selection
from sklearn import neighbors
from sklearn import neural_network
from sklearn import gaussian_process
import sklearn.metrics
import model_utils
from model_preprocessing import preprocessing


def main():

    model_name = sys.argv[1]
    train_file = sys.argv[2]
    test_file = sys.argv[3]
    train_y_name = sys.argv[4]

    out_train_file = sys.argv[-2]
    out_test_file = sys.argv[-1]

    #DEBUG
    #train_file = "/Users/pgarcia-nietochanzuckerberg.com/projects/cell_type_transfer/pancreas/data/paper_alternative.loom"
    #test_file = "/Users/pgarcia-nietochanzuckerberg.com/projects/cell_type_transfer/pancreas/data/hca_cellTypes.loom"
    #train_y_name = "cellType"
    #model_name = "svm_SVC"
    #model_name = "RandomForest"
    #model_name = "LogisticRegression"
    #model_name = "MLPClassifier"
    #model_name = "GaussianProcess"
    #model_name = "KNeighbors"

    # Read files
    train_h5ad = sc.read(train_file)
    test_h5ad = sc.read(test_file)

    train_h5ad_original = train_h5ad.copy()
    test_h5ad_original = test_h5ad.copy()

    # Process
    train_h5ad = preprocessing(train_h5ad, do_log1p=True, do_min_cells_filter=True, do_select_variable_genes=True, do_min_genes_filter=False)
    test_h5ad = preprocessing(test_h5ad, do_log1p=True, do_min_cells_filter=False, do_select_variable_genes=False, do_min_genes_filter=False)

    # Select same genes
    train_h5ad, test_h5ad = model_utils.subset_annData(train_h5ad, test_h5ad)

    # Get matrix and labels
    if train_y_name in test_h5ad.obs.keys():
        x_train, y_train, x_test,  y_test = model_utils.get_mat_and_lab(x=train_h5ad, x_lab=train_y_name, y=test_h5ad, y_lab=train_y_name)
    else:
        x_train, y_train, x_test = model_utils.get_mat_and_lab(x=train_h5ad, x_lab=train_y_name, y=test_h5ad)

    # Train model
    model = get_model(model_name, x_train, y_train, n_iter=100)

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

def precision_multiclass(y_true, y_pred):
    return sklearn.metrics.precision_score(y_true, y_pred, average='micro')

def recall_multiclass(y_true, y_pred):
    return sklearn.metrics.recall_score(y_true, y_pred, average='micro')


def get_model(model_name: str,
              x_train,
              y_train: list,
              n_iter: int = 20,
              cv: int = 5,
              scoring: str = "accuracy",
              ):

    accepted_models = ["svm_SVC", "RandomForest", "LogisticRegression", "KNeighbors", "MLPClassifier", "GaussianProcess"]
    #metrics = ["accuracy", "precision", "recall"]
    metrics = {
        'accuracy': sklearn.metrics.make_scorer(sklearn.metrics.accuracy_score),
        'precision': sklearn.metrics.make_scorer(precision_multiclass),
        'recall': sklearn.metrics.make_scorer(recall_multiclass)
    }
    do_randomized_grid_search = True

    if model_name not in accepted_models:
        raise Exception("Model " + model_name + " is not recognized. Accepted models are: " + ", ".join(accepted_models))

    # Select grid
    if model_name == "svm_SVC":
        model_obj = svm.SVC()
        par_grid = model_utils.grid_dict_svm_SVC(gamma_factor=1 / (x_train.shape[1] * x_train.var()),
                                                 for_randomized=do_randomized_grid_search)

    if model_name == "RandomForest":
        model_obj = ensemble.RandomForestClassifier()
        par_grid = model_utils.grid_dicts_RandomForest(for_randomized=do_randomized_grid_search)

    if model_name == "LogisticRegression":
        model_obj = linear_model.LogisticRegression()
        par_grid = model_utils.grid_dicts_LogisticRegression(for_randomized=do_randomized_grid_search)

    if model_name == "KNeighbors":
        model_obj = neighbors.KNeighborsClassifier()
        par_grid = model_utils.grid_dicts_KNeighbors(for_randomized=do_randomized_grid_search)

    if model_name == "MLPClassifier":
        model_obj = neural_network.MLPClassifier()
        par_grid = model_utils.grid_dicts_MLPClassifier(for_randomized=do_randomized_grid_search)

    if model_name == "GaussianProcess":
        model_obj = gaussian_process.GaussianProcessClassifier()
        par_grid = model_utils.grid_dicts_GaussianProcess(for_randomized=do_randomized_grid_search)

    # Train model
    if do_randomized_grid_search:
        model = model_selection.RandomizedSearchCV(model_obj, par_grid, cv=cv, scoring=metrics, refit=scoring, n_jobs=-1,
                                                   random_state=42, n_iter=n_iter)
    else:
        model = model_selection.GridSearchCV(model_obj, par_grid, cv=cv, scoring=metrics, refit=scoring, n_jobs=-1)

    model.fit(x_train, y_train)

    # Reporting results back
    print("Best parameters for model", model_name)
    for key in model.best_params_:
        print("\t", key, ":", model.best_params_[key], sep="")

    return model


if __name__ == "__main__":
    main()


