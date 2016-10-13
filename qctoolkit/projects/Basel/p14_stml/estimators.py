from __future__ import absolute_import, division, print_function, unicode_literals
from future_builtins import *
# code from michael.eickenberg@gmail.com 
# meant to be sklearn module

import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)

from sklearn.base import BaseEstimator, TransformerMixin
import numpy as np
import scipy

class OrthogonalLeastSquares(BaseEstimator, TransformerMixin):

    def __init__(self, max_num_atoms=None, tol=1e-6):
        self.max_num_atoms = max_num_atoms
        self.tol = tol

    def _normalize_dictionary(self):
        norms = np.sqrt(
            np.einsum('ij, ij -> j', self.dictionary, self.dictionary))
        norms[norms == 0] = 1.
        self.dictionary /= norms
        #print(norms)

    def _orthogonalize_dictionary(self, wrt):
        # assumes dictionary and wrt are of unit norm
        similarities = wrt.dot(self.dictionary)
        self.dictionary -= similarities * wrt[:, np.newaxis]
        # assert((wrt.dot(self.dictionary) == 0).all())
        

    def _fit(self, X, y):
        self.dictionary = X.copy()
        self._normalize_dictionary()

        n, p = X.shape
        max_num_atoms = self.max_num_atoms
        if max_num_atoms is None:
            max_num_atoms = np.min(X.shape)

        self.component_list_ = []
        self.components = np.zeros([n, max_num_atoms], dtype=X.dtype)

        residual = y.copy().ravel().astype('float')
        residual /= np.sqrt(np.dot(residual, residual))

        self.residuals = np.zeros_like(self.components)
        for m in range(max_num_atoms):
            if scipy.linalg.norm(residual) < self.tol:
                #stop
                break
            self.residuals[:, m] = residual

            similarities = residual.dot(self.dictionary)
            strongest = np.abs(similarities).argmax()
            component = self.dictionary[:, strongest].copy()
            self.components[:, m] = component
            residual -= similarities[strongest] * component
            self.component_list_.append(strongest)
            self.dictionary[:, strongest] = 0.
            self._orthogonalize_dictionary(component)
            # assert (component == self.components[:, m]).all()
            # assert((np.abs(component.dot(self.dictionary)) < 1e-10).all())
            # assert(
                # (np.abs(self.components.T.dot(self.dictionary)) < 1e-10).all())
            if m < max_num_atoms - 1:
                self._normalize_dictionary()
            # assert(
            #     (np.abs(self.components.T.dot(self.dictionary)) < 1e-10).all())

        return self


from sklearn.linear_model import Ridge
#def ols_fit_predict(X_train, y_train, X_test,
#                    max_num_atoms=None,
#                    nums_atoms=None,
#                    ridge_penalty=0.):
#    if max_num_atoms is None:
#        if nums_atoms is not None:
#            max_num_atoms = max(nums_atoms)
#        else:
#            max_num_atoms = min(X_train.shape)
#    if nums_atoms is None:
#        nums_atoms = range(1, max_num_atoms + 1)
#
#    ols = OrthogonalLeastSquares(max_num_atoms=max_num_atoms)._fit(
#            X_train, y_train)
#    components = ols.component_list_
#    
#    ridge = Ridge(alpha=ridge_penalty)
#    predictions = list()
#    for num_atoms in nums_atoms:
#        X = X_train[:, components[:num_atoms]]
#        Xt = X_test[:, components[:num_atoms]]
#        y_pred = ridge.fit(X, y_train).predict(Xt)
#        predictions.append(y_pred)
#
#    return np.array(predictions).T
#

from sklearn.cross_validation import check_cv, cross_val_score, ShuffleSplit
from sklearn.externals.joblib import Parallel, delayed
from sklearn.utils import check_random_state
from numbers import Number
class OrthogonalLeastSquaresCV(BaseEstimator, TransformerMixin):

    def __init__(self, max_num_atoms=None, nums_atoms=None,
                 ridge_penalty=0., split=.5, ridge_cv=5,
                 scoring='mean_squared_error',
                 n_jobs=1, random_state=42):
        self.max_num_atoms = max_num_atoms
        self.nums_atoms = nums_atoms
        self.ridge_cv = ridge_cv
        self.scoring = scoring
        self.ridge_penalty = ridge_penalty
        self.split = split
        self.random_state = random_state


    def fit(self, X, y):
        if self.max_num_atoms is None:
            if self.nums_atoms is not None:
                max_num_atoms = max(self.nums_atoms)
            else:
                max_num_atoms = min(X.shape)
        else:
            max_num_atoms = self.max_num_atoms
        if self.nums_atoms is None:
            nums_atoms = range(1, max_num_atoms + 1)
        else:
            nums_atoms = self.nums_atoms
        if isinstance(self.split, Number):
            rng = check_random_state(self.random_state)
            split = next(iter(
                ShuffleSplit(X.shape[0], n_iter=1, test_size=self.split,
                             random_state=rng)))

        ols_train_set, ridge_cv_set = split
        self.ols = OrthogonalLeastSquares(
                    max_num_atoms=max_num_atoms)._fit(X[ols_train_set],
                                                      y[ols_train_set])
        components = self.ols.component_list_
        all_scores = []
        for num_atoms in nums_atoms:
            X_selected = X[:, components[:num_atoms]][ridge_cv_set]
            scores = cross_val_score(
                Ridge(alpha=self.ridge_penalty), X_selected, y[ridge_cv_set],
                      cv=self.ridge_cv, scoring=self.scoring)
            all_scores.append(scores)
        self.all_scores = np.array(all_scores)
        mean_scores = self.all_scores.mean(axis=1)
        best_model = mean_scores.argmax()
        self.best_model = components[:nums_atoms[best_model]]
        self.estimator_ = Ridge(alpha=self.ridge_penalty).fit(
            X[:, self.best_model], y)

    def predict(self, X):
        X_selected = X[:, self.best_model]
        return self.estimator_.predict(X_selected)


from sklearn.preprocessing import StandardScaler
from sklearn.metrics import mean_squared_error, mean_absolute_error
from sklearn.linear_model import Ridge

def ols_path(Xtrain, ytrain, Xtest, ytest, ridge_penalty=.1,
             standardize=True, max_num_atoms=None, graph=None,
             return_components=False, return_lr_coefficients=False,
             return_ols_coefficients=False):

    if return_lr_coefficients:
        raise NotImplementedError("Didn't implement lr coefs yet")

    if graph is None:
        ols = OrthogonalLeastSquares(max_num_atoms=max_num_atoms)
        ols._fit(Xtrain, ytrain)
        components = ols.component_list_
        if return_ols_coefficients:
            ytrain_norm = ytrain / np.linalg.norm(ytrain)
            ols_coefficients = ols.components.T.dot(ytrain_norm)
    else:
        components = dag_constrained_ols(Xtrain, ytrain,
                                         max_num_atoms=max_num_atoms,
                                         graph=graph)

    if standardize:
        scaler = StandardScaler()
        Xtrain = scaler.fit_transform(Xtrain)
        Xtest = scaler.transform(Xtest)

    ridge = Ridge(alpha=ridge_penalty)

    mse_path, mae_path = [], []
    for i in range(len(components)):
        X = Xtrain[:, components[:i + 1]]
        Xt = Xtest[:, components[:i + 1]]
        ypred = ridge.fit(X, ytrain).predict(Xt)
        mse_path.append(mean_squared_error(ytest, ypred))
        mae_path.append(mean_absolute_error(ytest, ypred))

    output = np.array(mse_path), np.array(mae_path)
    if return_components:
        output = output + (components,)
    if return_ols_coefficients:
        output = output + (ols_coefficients,)

    return output


from sklearn.cross_validation import check_cv
def ols_path_cv(X, y, cv, ridge_penalty=.1, standardize=True, 
                max_num_atoms=None, graph=None,
                return_components=False,
                return_ols_coefficients=False):

    mse_paths = []
    mae_paths = []
    if return_components:
        componentss = []
    if return_ols_coefficients:
        ols_coefficientss = []

    for train, test in cv:
        output = ols_path(X[train], y[train], X[test], y[test],
                              ridge_penalty=ridge_penalty, 
                              standardize=standardize,
                              max_num_atoms=max_num_atoms,
                              graph=graph,
                              return_components=return_components,
                              return_ols_coefficients=return_ols_coefficients)
        mse_path, mae_path = output[:2]    
        mse_paths.append(mse_path)
        mae_paths.append(mae_path)
        if return_components:
            components = output[2]
            componentss.append(components)
        if return_ols_coefficients:
            ols_coefficients = output[-1]
            ols_coefficientss.append(ols_coefficients)
    all_output = (mse_paths, mae_paths)
    if return_components:
        all_output = all_output + (componentss,)
    if return_ols_coefficients:
        all_output = all_output + (ols_coefficientss,)
    return all_output




"""
    def _normalize_dictionary(self):
        norms = np.sqrt(
            np.einsum('ij, ij -> j', self.dictionary, self.dictionary))
        norms[norms == 0] = 1.
        self.dictionary /= norms
        #print(norms)

    def _orthogonalize_dictionary(self, wrt):
        # assumes dictionary and wrt are of unit norm
        similarities = wrt.dot(self.dictionary)
        self.dictionary -= similarities * wrt[:, np.newaxis]
        # assert((wrt.dot(self.dictionary) == 0).all())

"""


def _normalize_dictionary(dictionary, inplace=True):
    norms = np.linalg.norm(dictionary, axis=0)
    norms[norms == 0.] = 1.
    if not inplace:
        dictionary = dictionary.copy()
    dictionary /= norms
    return dictionary


def _orthogonalize_dictionary(dictionary, wrt, inplace=True):
    # assumes that dictionary and wrt are of unit norm

    similarities = wrt.dot(dictionary)
    if not inplace:
        dictionary = dictionary.copy()
    dictionary -= similarities * wrt[:, np.newaxis]
    return dictionary


def _full_graph(p):
    return dict(start=np.arange(p).astype('int'))


def dag_constrained_ols(X, y, graph=None, max_num_atoms=None, tol=1e-6):

    n, p = X.shape

    if max_num_atoms is None:
        max_num_atoms = p
    if graph is None:
        graph = _full_graph(p)

    authorized_atoms = np.unique(graph['start'])
    max_iterations = min(n, p, max_num_atoms)
    component_list = []
    component_values = np.zeros((n, max_iterations))
    residual = y.copy().ravel().astype('float')
    residual /= scipy.linalg.norm(residual)
    dictionary = _normalize_dictionary(X, inplace=False)


    for i in range(max_iterations):
        if scipy.linalg.norm(residual) < tol:
            break
        similarities = residual.dot(dictionary[:, authorized_atoms])
        strongest_index_ = np.abs(similarities).argmax()
        strongest_index = authorized_atoms[strongest_index_]
        component_list.append(strongest_index)
        component = dictionary[:, strongest_index].copy()
        component_values[:, i] = component
        residual -= similarities[strongest_index_] * component
        dictionary[:, strongest_index] = 0.  # not really necessary
        _orthogonalize_dictionary(dictionary, wrt=component)
        if np.isnan(dictionary).any(): stop
        if i < max_iterations - 1:
            _normalize_dictionary(dictionary)
        if np.isnan(dictionary).any(): stop
        new_authorized_atoms = np.array(graph.get(strongest_index, []),
                                        dtype='int')
        authorized_atoms = np.unique(
            np.concatenate([authorized_atoms,
                            new_authorized_atoms]))
    return component_list


# Maybe these functions should go somewhere else, but right now they are
# clearly helpers for the DAG ols

def order2_index(j1, r1, j2, r2, J, L):
    if j2 <= j1:
        return None
    offset = L ** 2 * (j1 * J - (j1 * (j1 + 1) // 2)) + L * (J - j1 - 1) * r1
    index = L * (j2 - j1 - 1) + r2

    return offset + index


def order2_index_snlight(j1, r1, j2, r2, J, L):
    # this gives indices for the ordering by layer 2
    # as regrouped by scatnet light

    if j2 <= j1:
        return None

    offset = (((j2 - 1) * j2) // 2) * L ** 2
    index = j1 * L ** 2 + r1 * L + r2

    return offset + index


#def order2_indices(j1s, r1s, j2s, r2s, J, L):
#
#    j1s, r1s, j2s, r2s = map(np.atleast_1d, map(np.array, (j1s, r1s, j2s, r2s)))
#
#    for j1, r1, j2, r2 in np.broadcast(j1s, r1s, j2s, r2s):
#        pass
#

def scattering_graph(n_scales, n_orientations, offset=0, scatnet_light=True):
    """Assuming the scattering coefficients are ordered as j1, r1, j2, r2,
       create a graph freeing second order coefs after selection of first."""
    
    # add layer 0
    graph = dict(start=np.array([0]) + offset)

    # add layer 1
    graph[0 + offset] = 1 + np.arange(n_scales * n_orientations) + offset

    if scatnet_light:
        indexer = order2_index_snlight
    else:
        indexer = order2_index

    # add layer 2
    l1_counter = 0
    for j1 in range(n_scales):
        for r1 in range(n_orientations):
            i = j1 * n_orientations + r1
            assert i == l1_counter
            index_list = []
            for j2 in range(j1 + 1, n_scales):
                for r2 in range(n_orientations):
                    index_list.append(indexer(j1, r1, j2, r2,
                                              n_scales, n_orientations))
            graph[i + offset] = np.array(index_list).astype('int') + offset
            l1_counter += 1
    return graph


def multi_feature_graph(n_scales, n_orientations, n_features,
                        scatnet_light=True): 
    # if several features are extracted per scattering coefficient:
    # matrices are hstacked, 0-levels and first levels are taken together

    scattering_size = (n_scales * (n_scales - 1) // 2 * n_orientations ** 2
                       + n_scales * n_orientations + 1)

    graphs = [scattering_graph(n_scales, n_orientations,
                               offset=i * scattering_size,
                               scatnet_light=scatnet_light)
              for i in range(n_features)]

    # now mix the graphs
    graph = dict()
    for g in graphs:
        for k, v in g.items():
            if k != 'start':
                graph[k] = v

    layers0 = np.concatenate([g['start'] for g in graphs])
    layers1 = np.concatenate([g[g['start'][0]] for g in graphs])

    graph['start'] = layers0
    for k in layers0:
        graph[k] = layers1

    return graph

    




"""
        self.component_list_ = []
        self.components = np.zeros([n, max_num_atoms], dtype=X.dtype)

        residual = y.copy().ravel().astype('float')
        residual /= np.sqrt(np.dot(residual, residual))

        self.residuals = np.zeros_like(self.components)
        for m in range(max_num_atoms):
            if scipy.linalg.norm(residual) < self.tol:
                #stop
                break
            self.residuals[:, m] = residual

            similarities = residual.dot(self.dictionary)
            strongest = np.abs(similarities).argmax()
            component = self.dictionary[:, strongest].copy()
            self.components[:, m] = component
            residual -= similarities[strongest] * component
            self.component_list_.append(strongest)
            self.dictionary[:, strongest] = 0.
            self._orthogonalize_dictionary(component)
            # assert (component == self.components[:, m]).all()
            # assert((np.abs(component.dot(self.dictionary)) < 1e-10).all())
            # assert(
                # (np.abs(self.components.T.dot(self.dictionary)) < 1e-10).all())
            if m < max_num_atoms - 1:
                self._normalize_dictionary()
            # assert(
            #     (np.abs(self.components.T.dot(self.dictionary)) < 1e-10).all())



"""


    # make a strongly constraining graph
def _linked_list(p, backwards=False):
    points = np.arange(p)
    if backwards:
        points = points[::-1]

    linked_list = dict(start=points[0])
    for start, end in zip(points[:-1], points[1:]):
        linked_list[start] = np.array([end])
    return linked_list


from scipy.linalg import lstsq
def bilinear_least_squares(X, y, b0=None, n_iter=10, fit_intercept=True):
    """assumes X.shape = n_samples, n_matrices, h, wi
       and does linear regression as a sum of rank 1 matrices"""

    if X.ndim == 3:
        X = X[:, np.newaxis]
    n_samples, n_matrices, n_feat_a, n_feat_b = X.shape
    if b0 is None:
        b0 = np.ones((n_matrices, n_feat_b)) / n_feat_b
    b = b0.copy()

    if fit_intercept:
        X_mean, y_mean = X.mean(0), y.mean()
        X = X - X_mean
        y = y - y_mean

    for i in range(n_iter):
        a_estimation_matrix = np.einsum(
                                 "ijkl, jl -> ijk", X, b).reshape(n_samples, -1)
        a = lstsq(a_estimation_matrix, y)[0].reshape(n_matrices, n_feat_a)
        b_estimation_matrix = np.einsum(
                                 "ijkl, jk -> ijl", X, a).reshape(n_samples, -1)
        b = lstsq(b_estimation_matrix, y)[0].reshape(n_matrices, n_feat_b)


    if fit_intercept:
        intercept = y_mean - np.einsum("jkl, jk, jl", X_mean, a, b)
        return a, b, intercept

    return a, b


class BilinearRegression(BaseEstimator):

    def __init__(self, fit_intercept=True):
        self.fit_intercept = fit_intercept

    def fit(self, X, y):
        output = bilinear_least_squares(X, y, fit_intercept=self.fit_intercept)

        if self.fit_intercept:
            self.coef_a_, self.coef_b_, self.intercept_ = output
        else:
            self.coef_a_, self.coef_b_ = output
            self.intercept_ = 0.
        return self

    def predict(self, X):
        return np.einsum('ijkl, jk, jl -> i', X, self.coef_a_, self.coef_b_) + self.intercept_

def test_bilinear_regression():
    from numpy.testing import assert_array_almost_equal, assert_almost_equal
    a = np.arange(1, 5)
    b = np.arange(6, 10)
    abmat = a[:, np.newaxis] * b
    X = np.random.randn(100, abmat.size)
    y = X.dot(abmat.ravel())
    X.shape = X.shape[0], 1, len(a), len(b)

    a_, b_ = bilinear_least_squares(X, y, fit_intercept=False)
    a_b_mat = a_[0][:, np.newaxis] * b_[0]
    assert_array_almost_equal(a_b_mat, abmat)
    offset = 2
    a_, b_, c_ = bilinear_least_squares(X, y + offset, fit_intercept=True)
    assert_almost_equal(offset, c_)
    #return a, b, offset, a_, b_, c_


def test_bilinear_regression_estimator():
    from numpy.testing import assert_array_almost_equal, assert_almost_equal
    a = np.arange(1, 5)
    b = np.arange(6, 10)
    abmat = a[:, np.newaxis] * b
    X = np.random.randn(100, abmat.size)
    y = X.dot(abmat.ravel())

    from sklearn.linear_model import LinearRegression
    lr = LinearRegression(fit_intercept=False)
    br = BilinearRegression(fit_intercept=False)

    lr.fit(X, y)
    br.fit(X.reshape(X.shape[0], 1, 1, -1), y)
    br_coef = br.coef_b_.ravel() * br.coef_a_.ravel()[0]
    assert_array_almost_equal(lr.coef_, br_coef)

    br.fit(X.reshape(X.shape[0], 1, -1, 1), y)
    br_coef = br.coef_a_.ravel() * br.coef_b_.ravel()[0]
    assert_array_almost_equal(lr.coef_, br_coef)

    offset = 2
    lr.fit_intercept = True
    br.fit_intercept = True
    lr.fit(X, y + offset)
    X.shape = X.shape[0], 1, len(a), len(b)
    br.fit(X, y + offset)

    assert_almost_equal(br.intercept_, lr.intercept_)



    # first test wether it falls back to linear regression

if __name__ == "__main__":
#    est = OrthogonalLeastSquares(max_num_atoms=2)
#
#    X = np.random.randn(10, 100)
#    y = np.random.randn(10)
#
#    est._fit(X, y)
#
#     # test graph ols
#    print('dag constrained ols')
#    components = dag_constrained_ols(X, y, max_num_atoms=2)
#
#    print(est.component_list_)
#    print(components)
#
#    from sklearn.datasets import make_regression
#    X, y = make_regression(n_samples=200, n_features=500, n_informative=5, noise=1.)
#
#    # Watch out! Fitting OLS and cross-validating for correct number of atoms
#    # on the same set leads to wrong model selection!
#    est_cv = OrthogonalLeastSquaresCV()
#    est_cv.fit(X, y)
#    from sklearn.cross_validation import cross_val_score
##    scores = cross_val_score(est_cv, X, y, cv=10, scoring='mean_squared_error')
#
#    # test graph ols
#    print('dag constrained ols')
#    components = dag_constrained_ols(X, y)
#    est.max_num_atoms = None
#    est._fit(X, y)
#
#    print(est.component_list_)
#    print(components)
#
#
#    components2 = dag_constrained_ols(X, y, graph=_linked_list(X.shape[1]))
#    components3 = dag_constrained_ols(X, y, graph=_linked_list(X.shape[1], backwards=True))
#    print(components2)
#    print(components3)
#
#    metas = []
#    indices = []
#    for j1 in range(6):
#        for r1 in range(8):
#            for j2 in range(j1 + 1, 6):
#                for r2 in range(8):
#                    metas.append((j1, r1, j2, r2))
#                    indices.append(order2_index(j1, r1, j2, r2, 6, 8))
#
#
#    metas_snl = []
#    indices_snl = []
#    for j2 in range(1, 6):
#        for j1 in range(j2):
#            for r1 in range(8):
#                for r2 in range(8):
#                    metas_snl.append((j1, r1, j2, r2))
#                    indices_snl.append(
#                        order2_index_snlight(j1, r1, j2, r2, 6, 8))

    test_bilinear_regression()
    test_bilinear_regression_estimator()

