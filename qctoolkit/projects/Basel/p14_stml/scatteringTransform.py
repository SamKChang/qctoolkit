import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)
import qctoolkit as qtk
from numbers import Number
import pickle
import gzip
import copy
import numpy as np
from numpy.polynomial.hermite_e import hermeval
import matplotlib.pyplot as plt
import pkgutil
from scipy.special import binom
fft_eggs_loader = pkgutil.find_loader('pyfftw')
fft_found = fft_eggs_loader is not None
if fft_found:
    import pyfftw.interfaces.numpy_fft as fft
skl_eggs_loader = pkgutil.find_loader('pyfftw')
skl_found = skl_eggs_loader is not None
if skl_found:
    from sklearn.cross_validation import ShuffleSplit
    from sklearn.cross_validation import cross_val_score
    from sklearn.linear_model import Ridge
    from sklearn.kernel_ridge import KernelRidge
    from sklearn.externals.joblib import Parallel, delayed

def make_dirac_densities(x, z, grid_step=.01,
                         left_bound=None,
                         right_bound=None,
                         padding=0.,
                         dirac_width=0.001,
                         keep_integral='l1'):
    if left_bound is None:
        left_bound = x.min() - padding
    if right_bound is None:
        right_bound = x.max() + padding
    
    grid = np.mgrid[left_bound:right_bound + grid_step:grid_step]
    distances = (x[..., np.newaxis] - grid[np.newaxis, np.newaxis]) ** 2
    gaussians = np.exp(-.5 * distances / dirac_width ** 2) / np.sqrt(2 * np.pi * dirac_width ** 2)
    if keep_integral == 'l1':
        gaussian_integrals = gaussians.sum(axis=2)
        gaussians /= gaussian_integrals[..., np.newaxis]
    elif keep_integral == 'l2':
        norm = np.sqrt((gaussians ** 2).sum(axis=2))
        gaussians /= norm[..., np.newaxis]
    densities = (gaussians * z[..., np.newaxis]).sum(axis=1)
    
    return densities

def get_valence(Z):
    shells = np.array([0, 2] + [8.] * 10)
    shellsc = np.cumsum(shells)
    Z_ = np.atleast_1d(Z).ravel()
    diffs = Z_ - shellsc[:, np.newaxis]
    #print(diffs)
    diffs[diffs < 0] = np.nan
    values = np.nanmin(diffs, axis=0).astype(int)
    if isinstance(Z, Number):
        return values[0]
    return values.reshape(Z.shape)

def gabor_k(k, grid_x, xi0=np.pi * 3 / 4, sigma0=.5, scale=0):
    xi = xi0 / 2 ** scale
    sigma = sigma0 * 2 ** scale
    k_choose_j = binom(k, np.arange(k + 1))
    minus_one_power = (-1) ** np.arange(k + 1)
    i_xi_power = (1j * xi * sigma) ** (k - np.arange(k + 1))
    herm_coef = k_choose_j * minus_one_power * i_xi_power
    hermite_eval = hermeval(grid_x / sigma, herm_coef)
    gauss = np.exp(-.5 * (grid_x/sigma) ** 2) / np.sqrt(2 * np.pi * sigma ** 2) / sigma ** k
    wave = np.exp(1j * xi * grid_x)
    return hermite_eval * gauss * wave

def trimData(data, start, end):
    data_all = copy.deepcopy(data)
    for k in data_all.keys():
        data[k] = np.array(data_all[k])[start:end]
    return data

def loadData(fname, batch = (1, 1), return_range = False, size = None):
#    # hack from http://stackoverflow.com/questions/11305790
#    # /pickle-incompatability-of-numpy-arrays-between-python-2-and-3
#    with open(fname, 'rb') as f:
#        u = pickle._Unpickler(f)
#        u.encoding = 'latin1'
#        data = u.load()
#        #data = p
    with open(fname) as pkl:
        data = pickle.load(pkl)

    data['E'] = np.array(data['E'])
    
    positions = data['xyz'].sum(axis=2)
    p_min = positions.min()
    p_max = positions.max()
    
    if size is not None:
        data = trimData(data, 0, size)
    elif batch[0] != 1:
        N = len(data['E'])
        start = int(np.round((batch[1] - 1) * N / batch[0]))
        end = int(np.round((batch[1]) * N / batch[0]))
        data = trimData(data, start, end)
    
    if return_range:
        return data, [p_min, p_max]
    else:
        return data

def getSignal_1d(data, padding=10, p_range = None, grid_step = 0.01, components = ['all', 'valence', 'core']):
    
    global make_dirac_densities, get_valence
                
    val = get_valence(data['Z'])
    full = data['Z']
    cor = full - val
    
    positions = data['xyz'].sum(axis=2)
    if p_range is None:
        left_bound = positions.min()
        right_bound = positions.max()
    else:
        left_bound = p_range[0]
        right_bound = p_range[1]

    signal = []
    if 'all' in components or 'full' in components:
        full_densities = make_dirac_densities(positions, full, grid_step = grid_step,
                                              left_bound=left_bound - padding,
                                              right_bound=right_bound + padding)
        signal.append(full_densities)

    if 'val' in components or 'valence' in components:
        val_densities = make_dirac_densities(positions, val, grid_step = grid_step,
                                             left_bound=left_bound - padding,
                                             right_bound=right_bound + padding)
        signal.append(val_densities)

    if 'core' in components or 'cor' in components:
        cor_densities = make_dirac_densities(positions, cor, grid_step = grid_step,
                                             left_bound=left_bound - padding,
                                             right_bound=right_bound + padding)
        signal.append(cor_densities)
    
    # endpoint=False gives exactly the same results as np.arange
    x_coord = np.linspace(left_bound - padding, right_bound + padding, full_densities.shape[-1], endpoint=False)
    
            
    return np.stack(signal, axis = 1), x_coord

def getFilter_1d(signal, derivatives = np.array([0, 1, 2]), scales = np.arange(3, 10), wavelet = gabor_k):
    
    grid = np.arange(signal.shape[-1])
    grid -= len(grid) // 2
    
    wavelets = np.array([[wavelet(derivative, grid, scale=scale)
                          for scale in scales] for derivative in derivatives])
    wavelets /= np.abs(wavelets).sum(-1)[..., np.newaxis]
    return wavelets

def ifft(x):
    return fft.ifft(x)

delayed_ifft = delayed(ifft)

def stCoefs_1d(signals, filters, second_layer=False, parallel = False, block_size = 20):
    """
    compute 1D signals scattering coefficents
    input
      signals: 2D numpy array (signal_dimensions, signal_length)
      filters: 2D numpy array (filter_dinensions, signal_length)
    output
      filtered_signal: 
      3D numpy array (signal_dinensions, filter_dinensions, signal_length)
    for second layer, assuming scale dimension is the second last
    that is filters[..., scale, :]
    """
    
    filters = filters.view()
    signals = signals.view()
    
    if not second_layer:
        
        signal_length = signals.shape[-1]
        signal_shape = signals.shape[:-1]
        signal_size = signals[..., 0].size # same as np.prod(signal_shape)
        signals.shape = (signal_size, signal_length)
        filter_length = filters.shape[-1]
        filter_shape = filters.shape[:-1]
        filter_size = np.prod(filter_shape)
        #filter_size = functools.reduce(lambda x, y: x * y, filter_shape)
        filters.shape = (filter_size, filter_length)
        print("%d features" % (signal_size * filter_size / signal_shape[0]))

        # for first layer filter and signal are always 2D array
        f_signals = fft.fft(signals, axis = -1)
        f_filters = fft.fft(fft.ifftshift(filters, axes=(1,)), axis = -1)
        f_conv = f_signals[:, np.newaxis, :] * f_filters[np.newaxis, :, :]
        
    else:
        n_samples = signals.shape[0]
        signal_length = signals.shape[-1]
        n_scale = signals.shape[-2]
        signal_shape = list(signals.shape[1:-2])
        filter_shape = list(filters.shape[:-2])
        
        new_signal_shape = [n_samples]
        new_signal_shape.extend(signal_shape)
        new_signal_shape.extend(filter_shape)
        n_scale_new = n_scale * (n_scale - 1) // 2
        new_signal_shape.extend([n_scale_new, signal_length])
        feature_size = np.prod(new_signal_shape[1:-1])
        #conv_size = np.prod(new_signal_shape[:-1])
        new_singal_shape = tuple(new_signal_shape)
        print("%d features" % feature_size)
        
        f_signals = fft.fft(np.abs(signals), axis = -1)
        f_filters = fft.fft(fft.ifftshift(filters, axes=(-1,)), axis = -1)
        #f_conv = np.zeros(new_signal_shape, dtype=np.complex64)
        f_conv = np.zeros(new_signal_shape, dtype=np.complex)
        
        counter = 0
        for i in range(n_scale - 1):
            # scale is always at second last dimension
            f_filter = f_filters[..., i + 1:, :].view()
            f_signal = f_signals[..., i, :].view()
            n_extra_scales = f_filter.shape[-2]
            
            # convert signal and filter to same size for broadcasting
            # filter: [x1, x2, ..., scale, length] -> [1, 1, ..., x1, x2, ..., new_scale, length]
            # signal: [y1, y2, ..., scale, length] -> [y1, y2, ..., 1, 1, ..., new_scale, length]
            # result: [y1, y2, ..., x1, x2, ..., new_scale, length]
            filter_range = len(f_signal.shape) - 1
            signal_range = len(f_filter.shape) - 1
            for j in range(filter_range):
                f_filter = np.expand_dims(f_filter, axis = 0)
            for j in range(signal_range):
                f_signal = np.expand_dims(f_signal, axis = -2)
                
            f_conv_i = f_signal * f_filter
            
            # update every batch on second last dimension
            f_conv[..., counter:counter+n_extra_scales, :] = f_conv_i
            
            counter += n_extra_scales
        
        signal_shape = signal_shape[:-1]
        filter_shape = filter_shape[:-1]
        #f_conv.shape = (conv_size, signal_length)
        #f_conv = np.zeros(signal_size, filter_size, signal_length)

    if parallel:
        filtered = Parallel(n_jobs=-1)(delayed_ifft(f_conv[i:i + block_size]) 
                                       for i in range(0, len(f_conv), block_size))
        filtered = np.concatenate(filtered, axis = 0)
    else:
        filtered = fft.ifft(f_conv)
    
    #del f_conv
    if not second_layer:
        filtered.shape = signal_shape + filter_shape + (signal_length,)

    return filtered
    #return np.abs(filtered).astype('float32')
    #return np.abs(filtered)

def regressionMatrix(*filtered_list, **kwargs):
    if 'normalize' in kwargs:
        normalize = kwargs['normalize']
    else:
        normalize = False
    features = []
    for filtered in filtered_list:
        filtered = filtered.view()
        data_size = filtered.shape[-1]
        sample_size = filtered.shape[0]
        shape = filtered.size // data_size // sample_size
        filtered_reshaped = filtered.reshape(sample_size, shape, data_size)

        l1 = np.abs(filtered_reshaped).sum(-1)
        l2 = (np.abs(filtered_reshaped) ** 2).sum(-1)
        if normalize:
            l1 /= np.sqrt((l1 ** 2).sum(0))
            l2 /= np.sqrt((l2 ** 2).sum(0))
        features.append(np.hstack([l1, l2]))

#        l1_list = np.atleast_2d()
#        l2_list = np.atleast_2d()
#        for i in range(shape):
#            # integrate over signal length
#            l1 = np.abs(filtered_reshaped[:, i, :]).sum(-1)
#            # normalize over all samples
#            l2 = (np.abs(filtered_reshaped[:, i, :]) ** 2).sum(-1)
#            if normalize:
#                l1 /= np.sqrt((l1 ** 2).sum(0))
#                l2 /= np.sqrt((l2 ** 2).sum(0))
#            try:
#                l1_list = np.hstack([l1_list, l1[:, np.newaxis]])
#                l2_list = np.hstack([l2_list, l2[:, np.newaxis]])
#            except:
#                l1_list = l1[:, np.newaxis]
#                l2_list = l2[:, np.newaxis]
#        features.append(np.hstack([l1_list, l2_list]))

    return np.hstack(features)

def stModel_1d(fname, batch = 1, signal_setting = {}, filter_setting = {}):
    st_matrix_chunks = []
    E = []
    for i in range(1, batch+1):
        print("processing %d out of %d" % (i, batch))
        data, p_range = loadData(fname, 
                                 batch = (batch, i), 
                                 return_range = True)
        rho, x = getSignal_1d(data, p_range = p_range, **signal_setting)
        wlt = getFilter_1d(rho, **filter_setting)
        st_1st = stCoefs_1d(rho, wlt)
        st_2nd = stCoefs_1d(st_1st, wlt, 
                            second_layer=True, 
                            parallel=True)
        chunk = regressionMatrix(st_1st, st_2nd, normalize = False)
        st_matrix_chunks.append(chunk)
        E.extend(data['E'])
        del data, rho, wlt, st_1st, st_2nd
    st_matrix = np.vstack(st_matrix_chunks)
    norm = np.sqrt(np.sum(st_matrix ** 2, axis = -1))
    return np.vstack(st_matrix_chunks), E

def _get_best_components_from_folds(n_components, best_components):
    selected_best_components = np.array(best_components)[:, :n_components]
    selected_components = np.argsort(np.bincount(selected_best_components.ravel()))[::-1][:n_components]
    return selected_components

def stScore(data,
            n_samples_list = None,
            alphas = [1e-11],
            n_components_list = None,
            ols_components = None,
            regression_matrix = None,
            cv = None,
            threads = 1,
            st_setting = {},
           ):

    vec = data['E']

    qtk.report("ML.tools.stScores setting", "\n",
               "alphas:", alphas, "\n", 
               "n_components_list:", n_components_list, "\n",
               "ols_components:", ols_components, "\n",
               "n_samples_list:", n_samples_list, "\n",
               "cross_validation:", cv, "\n",
               "cv_threads:", threads, "\n",
               "final score format: [alphas, gammas, samples, cv]")

    selected_components_list = [
      _get_best_components_from_folds(n_components, ols_components)\
      for n_components in n_components_list]
    all_st_scores = []
    for alpha in alphas:
        alpha_scores = []
        all_st_scores.append(alpha_scores)
        for selected_components in selected_components_list:
            component_scores = []
            alpha_scores.append(component_scores)
            reg = regression_matrix[:, selected_components]
            for n_samples in n_samples_list:
                #print((len(selected_components), n_samples), end=" ")
                #sys.stdout.flush()
                cv_ = [(train[:n_samples], test) for train, test in cv]
                scores = cross_val_score(Ridge(alpha=1e-8), 
                         reg, vec, cv=cv_, n_jobs=threads, 
                         scoring='mean_absolute_error')
                component_scores.append(scores)
    return -np.array(all_st_scores)
    
