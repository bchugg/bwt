import numpy as np
import sys
from params import score_matrix, alphabet, d, e
from queue import Queue
from local_align import localAlignRNA


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
		x and y are two input strings of DNA
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
		# indices = Queue()
		# indices.put([0,0])
		
		
		while k > 0 and not q.empty():
			seq_and_off = q.get()
			align = localAlignRNA(seq_and_off[0], seq_and_off[1])
			#inds = indices.get()
			align.computeAlignment()
			if  not align.emptyAlign():
				self.alignments[self.k-k] = align.align_info
				self.split(align.align_info, q, seq_and_off[0], seq_and_off[1])
				#self.splitAlignments(align.align_info, subs[0], subs[1], q, indices, inds)
				k -= 1

		if k != 0:
			opt = self.k-k
			print("\nWarning: Unable to find",self.k,"optimal alignments")
			print("Number of alignments found:",opt)
			print("Resetting Optimal k to",opt, end=' ')
			self.k = opt
			print("... ... DONE\n")

		#Sort by index of x
		#self.alignments = sorted(self.alignments,key=itemgetter(1))

	def split(self, info, q, seq, offset):
		m1l, m1r, m2l, m2r = info[0][0],info[0][1],info[1][0],info[1][1]
		if m1l > m2l:
			m1l, m2l = m2l, m1l
			m1r, m2r = m2r, m1r

		# Left subsequence
		sl = self.x[offset:m1l]
		if len(sl) >0 : q.put([sl, offset])
		#print('Putting', offset, '-', m1l, 'on queue')

		# Middle subsequence
		if m1r < m2l:
			sm = self.x[m1r:m2l]
			q.put([sm,m1r])
			#print('Putting', m1r, '-', m1l, 'on queue')
		# Right subsequence
		sr = self.x[m2r:offset+len(seq)]
		if len(sr) > 0 : q.put([sr,m2r])
		#print('Putting', m2r, '-', offset + len(seq), 'on queue')


	def printAlignments(self):
		""" 
		Print Alignment. If index argument specified, print only 
		that alignment. Otherwise print all. 
		"""
		for i in range(self.k):
				localAlignRNA.printAlignment(self.alignments[i])
		



	
