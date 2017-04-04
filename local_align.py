import numpy as np
import sys
from operator import itemgetter
from params import score_matrix, alphabet
from queue import Queue


"""
Implements class to perform local alignment
between an RNA sequence and itself
"""


""" 
Global Vars
"""
# minimum score
ms = -1*sys.maxsize
# score matrix
sm = score_matrix


class localAlignRNA:
	"""
	A class to solve local Alignment. Computes k-optimal disjoint
	local alignments 
	"""
	
	def __init__(self, x, offset=0):
		"""
		x with length n, y with legth m are two sequences to compare
		"""
		self.x = x
		self.y = x[::-1]
		self.n = len(x)+1
		self.offset = offset
		

	def computeAlignment(self):
		"""
		Fill out alignment matrices for local alignment between
		sequences x and y. Return alignment, indices, and score
		"""
		# Initialize 
		n = self.n  
		self.M = np.zeros((n,n))
		for i in range(1, n):
			b = self.findIndex(i-1,self.y)
			for j in range(1, n):
				a = self.findIndex(j-1,self.x)
				choice = [self.M[i-1,j-1]+sm[a,b],0]
				self.M[i,j] = max(choice)

		self.traceback()

	def findIndex(self, index, sequence):
		""" 
		Find index in alphabet of sequence[index]
		"""
		try:
			return alphabet.index(sequence[index])
		except ValueError:
			return alphabet.index('N')
	
	
	def traceback(self):
		"""
		Traceback procedure for computing alignment
		
		Returns [align, i, j, score]
		align[0] is subsequence of x
		align[1] is subsequence of y
		j is index in x where local alignment begins (x[i])
		i is index of y where local alignments begins (y[j])
		"""
		align = ["", ""]
		score = np.max(self.M)
		am = np.argmax(self.M)
		j = am % self.M.shape[1]
		i = int(am / self.M.shape[1])
		while self.M[i,j] != 0:
			align[0] = self.x[j-1] + align[0]
			align[1] = self.y[i-1] + align[1]
			self.M[i,j] = 0
			i -= 1
			j -= 1

		#self.align_info = [[align[0], align[1]], j,i, score]
		self.assemble([align[0], align[1]], j,i, score)

	def assemble(self, align, xind, yind, score):
		"""
		Map alignments back onto original sequence x. Return alignment 
		and index info
		"""
		lrx = yind
		rrx = yind+len(align[1])
		m1_ind = [self.offset+xind, self.offset+xind+len(align[0])]
		m2_ind = [self.offset+self.n-rrx-1, self.offset+self.n-lrx-1]
		self.align_info = [m1_ind, m2_ind, score, align]

	def emptyAlign(self):
		""" 
		Return true if empty alignment
		"""
		m1_length = self.align_info[0][0] - self.align_info[0][1]
		m2_length = self.align_info[1][0] - self.align_info[1][1]
		return (m1_length == 0 or m2_length == 0)

	
	@staticmethod
	def printAlignment(alignment):
		"""
		Print single alignment
		"""
		print("SCORE: ", alignment[2])
		print('Match', alignment[0][0], '-', alignment[0][1])
		print('With', alignment[1][0], '-', alignment[1][1])
		print("Alignment: ")
		print(alignment[3][0][0:60])
		print(alignment[3][1][0:60])
		print('\n')




