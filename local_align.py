import numpy as np
import sys
from operator import itemgetter
from params import score_matrix, alphabet, d, e
from queue import Queue


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
	
	def __init__(self, x, y):
		"""
		x and y are two input strings of DNA
		"""
		self.x = x
		self.y = y

	def kAlignments(self, k):
		""" 
		Compute k Optimal disjoint local alignments between sequences x and y
		""" 
		self.k = k
		self.alignments = [[] for _ in range(k)]
		
		q = Queue()
		q.put([self.x, self.y])
		indices = Queue()
		indices.put([0,0])
		
		
		while k > 0 and not q.empty():
			subs = q.get()
			align = localAlign(subs[0], subs[1])
			inds = indices.get()
			align.computeAlignment(inds[0], inds[1])
			if  not align.emptyAlign():
				self.alignments[self.k-k] = align.align_info
				self.splitAlignments(align.align_info, subs[0], subs[1], q, indices, inds)
				k -= 1

		if k != 0:
			opt = self.k-k
			print("\nWarning: Unable to find",self.k,"optimal alignments")
			print("Number of alignments found:",opt)
			print("Resetting Optimal k to",opt, end=' ')
			self.k = opt
			print("... ... DONE\n")

		#Sort by index of x
		self.alignments = sorted(self.alignments,key=itemgetter(1))



	def splitAlignments(self, align, s1, s2, q, indices, inds):
		"""
		Split s1, s2 into subsequences before and after
		sequence in align, place on queue q
		align[0][0] is subsequence of s1, 
		align[0][1] is subsequence of s2
		"""
		sub_s1 = align[0][0]
		sub_s2 = align[0][1]
		# print("Just computed alignment:")
		# print(sub_s1)
		# print(sub_s2)
		# print("At indices:")
		# print('Index of x', align[1], 'Index of y:', align[2])

		# # Beginning sequences
		start_s1 = s1[0:s1.find(sub_s1)]
		start_s2 = s2[0:s2.find(sub_s2)]
		if len(start_s1) > 0 and len(start_s2) > 0:
			q.put([start_s1, start_s2])
			indices.put([inds[0], inds[1]])
			# print("Putting sequences on Stack:")
			# print(start_s1)
			# print(start_s2)
			# print('At indices:')
			# print('Index of x', inds[0], 'Index of y:', inds[1])

		# Ending sequences
		end_s1 = s1[s1.find(sub_s1)+len(sub_s1):]
		end_s2 = s2[s2.find(sub_s1)+len(sub_s2):]
		if len(end_s1) > 0 and len(end_s2) > 0:
			q.put([end_s1, end_s2])
			indices.put([inds[0]+len(sub_s1)+len(start_s1), inds[1]+len(sub_s2)+len(start_s2)])
			# print("Putting sequences on Stack:")
			# print(end_s1)
			# print(end_s2)
			# print('At indices:')
			# print('Index of x', inds[0]+len(sub_s1)+len(start_s1), 'Index of y:', inds[1]+len(sub_s2)+len(start_s2))



	def printAlignments(self):
		""" 
		Print Alignment. If index argument specified, print only 
		that alignment. Otherwise print all. 
		"""
		for i in range(self.k):
				localAlign.printAlignment(self.alignments[i])
		



	

class localAlign:
	"""
	A class to solve local Alignment. Computes k-optimal disjoint
	local alignments 
	"""
	
	def __init__(self, x, y):
		"""
		x with length n, y with legth m are two sequences to compare
		"""
		self.x = x
		self.y = y
		self.n = len(x)+1
		self.m = len(y)+1
	

	def computeAlignment(self, xoff, yoff):
		"""
		Fill out alignment matrices for local alignment between
		sequences x and y. Return alignment, indices, and score
		"""
		# Initialize 
		n = len(self.x) + 1
		m = len(self.y) + 1
		self.M = np.zeros((m,n))
		for i in range(1, m):
			b = self.findIndex(i-1,self.y)
			for j in range(1, n):
				a = self.findIndex(j-1,self.x)
				choice = [self.M[i-1,j-1]+sm[a,b],0]
				self.M[i,j] = max(choice)

		self.traceback(xoff, yoff)

	def findIndex(self, index, sequence):
		""" 
		Find index in alphabet of sequence[index]
		"""
		letter = sequence[index]
		if letter == 'U':
			letter = 'C'
		return alphabet.index(letter)

	
	
	def traceback(self, xoff, yoff):
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
		self.align_info = [[align[0], align[1]], j+xoff,i+yoff, score]

	def emptyAlign(self):
		""" 
		Return true if empty alignment
		"""
		a1 = self.align_info[0][0]
		a2 = self.align_info[0][1]
		return (len(a1) == 0 or len(a2) == 0)

	
	@staticmethod
	def printAlignment(alignment):
		"""
		Print single alignment
		"""
		print("SCORE: ", alignment[3])
		print('LENGTH: ', len(alignment[0][0]))
		print("Index of x:", alignment[1]) 
		print("Index of y:", alignment[2]) 
		print("Alignment: ")
		print("x: ", end='')
		print(alignment[0][0][0:60])
		print("y: ", end='')
		print(alignment[0][1][0:60])
		print('\n')




