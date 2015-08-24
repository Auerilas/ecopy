import numpy as np
from ..base_funcs import wt_mean

def wt_var(x, wt, bias=0):
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
	if bias==0:
		p1 = np.sum(wt_array*(x-w_mu)**2)
		p2 = 1/(1 - (wt_array**2).sum())
		w_var = p2*p1
	if bias==1:
		w_var = np.sum(wt_array*(x - w_mu)**2)
	return w_var
