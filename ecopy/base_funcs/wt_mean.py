import numpy as np

def wt_mean(x, wt=None):
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
	w_mu = np.sum(x*wt_array)
	return w_mu