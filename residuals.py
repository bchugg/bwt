from operator import itemgetter

"""
Module implementing residuals class to obtain residual subsequences
after computing k optimal local alignments
"""

class residuals:
	
	def __init__(self, x):
		self.x = x
		self.n = len(x)


	def getResiduals(self,inds):
		"""
		Compute residual subsequences of x which do not appear
		in any alignment, given list of indices of subsequences
		"""
		xrep = [ "free" for _ in range(len(self.x)) ]
		#Mark xrep
		for ind in inds:
			if xrep[ind[0]] == "free":
				xrep[ind[0]] = ind[1]
			else:
				if xrep[ind[0]] < ind[1]:
					xrep[ind[0]] = ind[1]
		residuals = []
		left = 0
		while left < self.n:
			partial = ""
			while xrep[left] == "free":
				
				partial += self.x[left]
				left += 1
				if left == self.n : break
			else:
				if len(partial)>0 : 
					residuals += [[partial]]
					partial = ""
				left = left + 1 if left == xrep[left] else xrep[left]
		if len(partial) >0 : residuals += [[partial]] 
		self.residuals = residuals






