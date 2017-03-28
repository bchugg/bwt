import numpy as np
import sys
from params import score_matrix, alphabet
from queue import Queue
from local_align import localAlignRNA
from functools import reduce


"""
Implements the two classes disjointAlignments
and localAlign to solve the k-disjoint local
alignment problem
"""


""" 
Global Vars
"""
# minimum score
ms = -1*sys.maxsize
# score matrix
sm = score_matrix


class disjointAlignments:
	""" 
	A class to compute the k-disjoint optimal local
	alignment between two input sequences, for any k
	"""
	
	def __init__(self, x):
		"""
		x is RNA strand
		"""
		self.x = x

	def kAlignments(self, k):
		""" 
		Compute k Optimal disjoint local alignments on RNA sequence x
		""" 
		self.k = k
		self.alignments = [[] for _ in range(k)]
		
		q = Queue()
		q.put([self.x, 0])
		
		while k > 0 and not q.empty():
			seq_and_off = q.get()
			align = localAlignRNA(seq_and_off[0], seq_and_off[1])
			align.computeAlignment()
			if  not align.emptyAlign():
				self.alignments[self.k-k] = align.align_info
				self.split(align.align_info, q, seq_and_off[0], seq_and_off[1])
				k -= 1

		if k != 0:
			opt = self.k-k
			print("\nWarning: Unable to find",self.k,"optimal alignments")
			print("Number of alignments found:",opt)
			print("Resetting Optimal k to",opt, end=' ')
			self.k = opt
			print("... ... DONE\n")

		
	def split(self, info, q, seq, offset):
		"""
		Map aligned regions onto self.x. Split into subsequences for recursion
		"""
		m1l, m1r, m2l, m2r = info[0][0],info[0][1],info[1][0],info[1][1]
		if m1l > m2l:
			m1l, m2l = m2l, m1l
			m1r, m2r = m2r, m1r

		# Left subsequence
		sl = self.x[offset:m1l]
		if len(sl) >0 : q.put([sl, offset])

		# Middle subsequence
		if m1r < m2l:
			sm = self.x[m1r:m2l]
			q.put([sm,m1r])

		# Right subsequence
		sr = self.x[m2r:offset+len(seq)]
		if len(sr) > 0 : q.put([sr,m2r])


	def printAlignments(self):
		""" 
		Print Alignment. If index argument specified, print only 
		that alignment. Otherwise print all. 
		"""
		for i in range(self.k):
				localAlignRNA.printAlignment(self.alignments[i])


	@staticmethod
	def scoreAlignments(alignments):
		"""
		Score self-alignments based on RNA scoring model
		"""
		score = 0
		for a in alignments:
			s, x, y = 0, a[3][0], a[3][1]
			for i in range(np.min([len(x),len(y)])):
				s += disjointAlignments.complement(x[i],y[i])
			score += s
		return score


	@staticmethod
	def complement(b1, b2):
		"""
		Return true iff b1 is complement base to b2
		"""
		if b1 == "A" and b2 == "U" : return 1
		if b1 == "U" and b2 == "A" : return 1
		if b1 == "C" and b2 == "G" : return 1
		if b1 == "G" and b2 == "C" : return 1
		return 0
		



	
