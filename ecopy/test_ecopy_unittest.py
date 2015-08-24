import unittest
import numpy as np
from ecopy import load_data, wt_mean, wt_var, wt_scale, diversity, rarefy

class TestECOPY(unittest.TestCase):

	def setUp(self):
		pass

	def test_wt_mean(self):
		x = [2, 3, 4, 5]
		w = [0.1, 0.1, 0.5, 0.3]
		res = wt_mean(x, w)
		self.assertEqual(res, 4.0)

	def test_wt_var(self):
		x = [2, 3, 4, 5]
		w = [0.1, 0.1, 0.5, 0.3]
		res = wt_var(x, w)
		self.assertEqual(res, 1.25)	

	def test_wt_scale(self):
		x = [2, 3, 4, 5]
		w = [0.1, 0.1, 0.5, 0.3]
		res = np.round(wt_scale(x,w), 2)
		truth = res == [-1.79, -0.89, 0, 0.89]
		self.assertEqual(truth.sum(), 4)

	def test_diversity(self):
		sp = np.array([0, 1, 2, 3, 0]).reshape(1,5)
		div = np.round(diversity(sp))
		self.assertEqual(div, 1)

	def test_rarefy(self):
		BCI = load_data('BCI')
		rareRich = np.round(rarefy(BCI, 'rarefy'))
		self.assertEqual(rareRich[1], 77)

if __name__ == '__main__':
	unittest.main()