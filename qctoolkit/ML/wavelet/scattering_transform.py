import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)

import numpy as np
import pkgutil

skl_eggs_loader = pkgutil.find_loader('sklearn')
skl_found = skl_eggs_loader is not None
fftw_eggs_loader = pkgutil.find_loader('pyfftw')
fftw_found = fftw_eggs_loader is not None
if skl_found:
  from sklearn.externals.joblib import Parallel, delayed
if fftw_found:
  import pyfftw.interfaces.numpy_fft as fft

def stCoefs_1d(signals, filters, block_size=100, second_order=False):
  """
  compute 1D signals scattering coefficents
  input
    signals: 2D numpy array (# samples, signal_length)
    filters: 2D numpy array (# filters, signal_length)
  output
    filtered_signal: 
    3D numpy array (# samples, # filters, signal_length)
  """

  signal_length = signals.shape[-1]
  signal_shape = signals.shape[:-1]
  signal_size = reduce(lambda x, y: x * y, signal_shape)
  signals.shape = (signal_size, signal_length)
  filter_length = filters.shape[-1]
  filter_shape = filters.shape[:-1]
  filter_size = reduce(lambda x, y: x * y, filter_shape)
  filters.shape = (filter_size, filter_length)

  f_signals = np.fft.fft(signals, axis=1)
  f_filters = np.fft.fft(fft.ifftshift(filters, axes=(1,)), axis=1)
  if not second_order:
    f_conv = f_signals[:, np.newaxis, :] * f_filters[np.newaxis, :, :]
  else:
    f_conv = np.zeros(signal_size, filter_size, signal_length)
    
  filtered = fft.ifft(f_conv, axis=2)

  return filtered
