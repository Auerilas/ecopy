import pandas as pd

def load_data(x):
	"""Loads data from online repository. Requires internet."""
	path = 'https://github.com/Auerilas/ecopy-data/raw/master/data/{0}.csv'
	downloadPath = path.format(x)
	loaded = pd.read_csv(downloadPath)
	return loaded