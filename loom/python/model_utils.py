'''
General functions to apply sklearn models on h5ad geneXcell matrices
'''


import scanpy as sc
import sklearn
import bbknn
from anndata import AnnData
from typing import Union, Optional, Tuple, Collection

def merge_knn(*args: AnnData) -> sc.AnnData:
    '''
    Wrapper for merging different datasets using batch balanced knn
    :param args: AnnDatas to merge
    :return: Concatenated AnnData
    '''

    annData = args[0].concatenate(args[1:], batch_key='knn_batch')

    sc.pp.log1p(annData)
    sc.pp.filter_genes(annData, min_cells=int(annData.n_obs * 0.05))
    sc.pp.filter_cells(annData, min_genes=5)
    sc.pp.highly_variable_genes(annData, flavor="seurat")
    sc.pp.normalize_total(annData, exclude_highly_expressed=True)
    sc.pp.pca(annData)

    sc.external.pp.bbknn(annData, batch_key='knn_batch', copy=False)
    sc.tl.umap()

    return annData

def merge_combat(*args: AnnData) -> sc.AnnData:
    '''
    Wrapper for merging different datasets using linear regressions
    uses ComBat
    :param args: AnnDatas to merge
    :return: Concatenated AnnData
    '''

    annData = args[0].concatenate(args[1:], batch_key='batch')
    sc.pp.combat(annData, key='batch')

    return annData

def subset_annData(x: AnnData, y: AnnData):
    '''
    Returns views of x and y in a list with the genes that are shared between the two datasets
    :param x: annData
    :param y: annData
    :return: (x_subset, y_subset

    '''

    x_in_y = [i in y.var_names.to_list() for i in x.var_names.to_list()]
    y_in_x = [i in x.var_names.to_list() for i in y.var_names.to_list()]

    if sum(x_in_y) == 0:
        raise Exception("No genes from x were found in y")

    if sum(y_in_x) == 0:
        raise Exception("No genes from y were found in x")

    x = x[:,x_in_y]
    y = y[:,y_in_x]

    return (x,y)

def validate_annData(x: AnnData, y: AnnData):

    '''
    Checks that two matrices can be used for modeling, usually one is the train data and the other is the test data
    Raises an exception if:
        - x and y don't have the same number of genes
        - x and y don't have the same gene names
    :param x: annData
    :param y: annData
    :return: None

    '''

    if x.n_vars != y.n_vars:
        raise Exception("Number of features in x and y are not the same")

    x_vars = x.var_names.to_list()
    y_vars = y.var_names.to_list()

    x_vars.sort()
    y_vars.sort()

    if x_vars != y_vars:
        raise Exception("Features are not the same")


def get_mat_and_lab(x: AnnData,
                    x_lab: str,
                    y: Optional[AnnData]=None,
                    y_lab: Optional[str]=None
                    ):
    '''
    Returns the matrix(ces) and labels from x and y with matching rows
    :param x: annData
    :param y: annData (Optional)
    :param x_lab: str - Name of observation label to return
    :param y_lab: str - Name of observation label to return
    :return: (x numpy.array, x list) wnen y is undecleated
    :return (x numpy.array, x list, y numpy.array) wnen y_lab is undeclared
    :return (x numpy.array, x list, y numpy.array, y list) wnen y_lab all are declared
    '''

    if y is None:
        return (x.X.toarray(), x.obs[x_lab].to_list())
    elif y_lab is None:
        # Check that matrices are compatible
        validate_annData(x,y)
        x_vars = x.var_names.to_list()
        y = y[:,x_vars]
        return (x.X.toarray(), x.obs[x_lab].to_list(), y.X.toarray())
    else:
        # Check that matrices are compatible
        validate_annData(x,y)
        x_vars = x.var_names.to_list()
        y = y[:,x_vars]
        return (x.X.toarray(), x.obs[x_lab].to_list(), y.X.toarray(), y.obs[y_lab].to_list())


'''
Dictonaries to use for grid searches on sklearn
'''

def grid_dict_svm_SVC(gamma_factor: float, for_randomized: bool = False):
    '''
    Creates dictonary ready to use in paramg_grid of model_selection.GridSearchCV for a SVC model
    :param gamma_factor: float - usually this is 1 / (number of features * variance of matrix)
    :return: list if for_randomized True, else dict
    '''
    x = [{'kernel': ['rbf'],
         'gamma': [gamma_factor * 1, gamma_factor * 1e1, gamma_factor * 1e2, gamma_factor * 1e3,
                   gamma_factor * 1e-1, gamma_factor * 1e-2, gamma_factor * 1e-3],
          'C': [1, 10, 50, 100, 500, 1000]
          },
         {'kernel': ['linear'],
          'C': [1, 10, 50, 100, 500, 1000]
         }
        ]

    if for_randomized:
        x = x[0]

    return x

def grid_dicts_RandomForest(for_randomized: bool = False):
    '''
    Creates dictonary ready to use in paramg_grid of model_selection.GridSearchCV for a random forest classifier
    :return: list if for_randomized True, else dict
    '''
    x = [{'bootstrap': [True, False],
          'max_depth': [10, 20, 30, 40, 50, 60, 70, 80, 90, 100, None],
          'max_features': ['auto', 'sqrt'],
          'min_samples_leaf': [1, 2, 4],
          'min_samples_split': [2, 5, 10],
          'n_estimators': [200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000]}
        ]

    if for_randomized:
        x = x[0]


    return x


def grid_dicts_LogisticRegression(for_randomized: bool = False):
    '''
    Creates dictonary ready to use in paramg_grid of model_selection.GridSearchCV for a random forest classifier
    :return: list if for_randomized True, else dict
    '''
    x = [{'penalty': ['elasticnet'],
          'tol': [1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1],
          'C': [1, 10, 20, 50, 70, 100, 200, 500, 700, 1000, 2000],
          'fit_intercept': [True, False],
          'solver': ['saga'],
          'multi_class': ['multinomial'],
          'l1_ratio': [0.00, 0.25, 0.50, 0.75, 1.00]
          }
         ]

    if for_randomized:
        x = x[0]


    return x

def grid_dicts_KNeighbors(for_randomized: bool = False):
    '''
    Creates dictonary ready to use in paramg_grid of model_selection.GridSearchCV for a kNN classifier
    :return: list if for_randomized True, else dict
    '''

    x = [{'n_neighbors': [2, 5, 10, 25, 50, 100, 200],
          'weights': ['uniform', 'distance'],
          'p': [1, 2],
          }
         ]

    if for_randomized:
        x = x[0]


    return x


def grid_dicts_MLPClassifier(for_randomized: bool = False):
    '''
    Creates dictonary ready to use in paramg_grid of model_selection.GridSearchCV for a Multi-layer Perceptron classifier
    :return: list if for_randomized True, else dict
    '''

    x = [{'hidden_layer_sizes': [(100), (100, 100), (100, 100, 100),
                                 (500), (500, 500), (500, 500, 500),
                                 (1000), (1000, 1000), (1000, 1000, 1000)],
          'activation': ['logistic', 'relu'],
          'solver': ['adam', 'lbfgs'],
          'learning_rate': ['constant', 'adaptive']
          }
         ]

    if for_randomized:
        x = x[0]


    return x

def grid_dicts_GaussianProcess(for_randomized: bool = False):
    '''
    Creates dictonary ready to use in paramg_grid of model_selection.GridSearchCV for a Multi-layer Perceptron classifier
    :return: list if for_randomized True, else dict
    '''

    x = [{'n_restarts_optimizer': [0, 1, 2, 5, 10],
          'max_iter_predict': [20, 50, 100, 200, 500],
          'warm_start': [True, False],
          'multi_class': ['one_vs_rest', 'one_vs_one']
          }
         ]

    if for_randomized:
        x = x[0]


    return x
