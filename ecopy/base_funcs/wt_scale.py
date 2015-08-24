import numpy as np
from ..base_funcs import wt_mean, wt_var

def wt_scale(x, wt, bias=0):
	if wt is None:
		wt = np.array([1]*len(x))
	x = np.array(x, 'float')
	wt_array = np.array(wt, 'float')
	if np.isnan(np.sum(x)):
		msg = 'vector contains null values'
		raise ValueError(msg)
	if x.shape != wt_array.shape:
		msg = 'weight vector must have equal dimensions as observations'
		raise ValueError(msg)
	wt_array = wt_array / wt_array.sum()
	w_mu = wt_mean(x, wt_array)
	w_var = wt_var(x, wt_array, bias)
	w_scale = (x - w_mu) / np.sqrt(w_var)
	return w_scale