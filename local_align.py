import numpy as np
import sys
from params import score_matrix, alphabet, d, e
from Queue import Queue


"""
Implements the two classes disjointAlignments
and localAlign to solve the k-disjoint local
alignment problem
"""


""" 
Global Vars
"""
# minimum score
ms = -1*sys.maxint
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
		self.alignments = ["" for _ in range(k)]
		q = Queue()
		q.put([self.x, self.y])
		
		while k > 0 and not q.empty():
			subs = q.get()
			align = localAlign(subs[0], subs[1])
			align.computeAlignment()
			self.alignments[self.k-k] = align.align_info
			self.splitAlignments(align.align_info, subs[0], subs[1], q)
			k -= 1

		if k != 0:
			opt = self.k-k
			print "\nWarning: Unable to find",self.k,"optimal alignments"
			print "Number of alignments found:",opt
			print "Resetting Optimal k to",opt,
			self.k = opt
			print "... ... DONE\n"



	def splitAlignments(self, align, s1, s2, q):
		"""
		Split s1, s2 into subsequences before and after
		sequence in align, place on queue q
		align[0][0] is subsequence of s1, 
		align[0][1] is subsequence of s2
		"""
		sub_s1 = align[0][0]
		sub_s2 = align[0][1]

		# Beginning sequences
		start_s1 = s1[0:align[1]]
		start_s2 = s2[0:align[2]]
		if len(start_s1) > 0 and len(start_s2) > 0:
			q.put([start_s1, start_s2])

		# Ending sequences
		end_s1 = s1[align[1]+len(sub_s1):len(s1)]
		end_s2 = s1[align[2]+len(sub_s2):len(s2)]
		if len(end_s1) > 0 and len(end_s2) > 0:
			q.put([start_s1, start_s2])


	
	def printAlignments(self, index=""):
		""" 
		Print Alignment. If index argument specified, print only 
		that alignment. Otherwise print all. 
		"""
		if index == "":
			for i in range(self.k):
				localAlign.printAlignment(self.alignments[i])
		else:
			localAlign.printAlignment(self.alignments[int(index)])




	

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


	def initializeMatrices(self,n,m):
		""" 
		Initialize matrices and pointer matrices for use in recurrence
		"""
		self.M, self.GX, self.GY = np.zeros((m,n)), np.zeros((m,n)), np.zeros((m,n))
		self.Mp = [ [ [] for _ in range(n) ] for _ in range(m) ]
		self.GXp = [ [ [] for _ in range(n) ] for _ in range(m) ]
		self.GYp = [ [ [] for _ in range(n) ] for _ in range(m) ]

		# Base Cases
		self.GX[0,0], self.GY[0,0] = ms, ms
		self.Mp[0][0], self.GXp[0][0], self.GYp[0][0] = 'STOP', 'STOP', 'STOP'
		for i in range(1,m):
			self.GY[i,0] = ms
			self.GX[i,0] = -d - (i-1)*e
			self.GXp[i][0], self.GYp[i][0], self.Mp[i][0] = 'STOP', 'STOP', 'STOP'
		for j in range(1,n):
			self.GX[0,j] = ms
			self.GY[0,j] = -d - (j-1)*e
			self.GXp[0][j], self.GYp[0][j], self.Mp[0][j] = 'STOP', 'STOP', 'STOP'
		

	
	def computeAlignment(self):
		"""
		Fill out alignment matrices for local alignment between
		sequences x and y. Return alignment, indices, and score
		"""
		# Initialize 
		n = len(self.x) + 1
		m = len(self.y) + 1
		self.initializeMatrices(n, m)

		for i in range(1, m):
			b = self.findIndex(i-1,self.y)
			for j in range(1, n):
				a = self.findIndex(j-1,self.x)
				# Choices in recurrence
				M_choice = np.append(sm[a,b]+np.array([self.M[i-1,j-1], self.GX[i-1,j-1], self.GY[i-1,j-1]]),[0])
				GX_choice = [ self.M[i-1,j]-d, self.GX[i-1,j]-e ] 
				GY_choice = [ self.M[i,j-1]-d, self.GY[i,j-1]-e ] 
				# Max arguments
				M_am = np.argmax(M_choice)
				GX_am = np.argmax(GX_choice)
				GY_am = np.argmax(GY_choice)
				# Compute new values and update pointer matrices
				self.M[i,j] = max(M_choice)
				if self.M[i,j] == 0:
					self.Mp[i][j] = 'STOP'
				else:
					self.Mp[i][j] = 'M' if M_am == 0 else ('GX' if M_am == 1 else 'GY')
				self.GX[i,j] = max(GX_choice)
				self.GXp[i][j] = 'M' if GX_am == 0 else 'GX'
				self.GY[i,j] = max(GY_choice)
				self.GYp[i][j] = 'M' if GY_am == 0 else 'GY'

		return self.traceback()
		


	def findIndex(self, index, sequence):
		""" 
		Find index in alphabet of sequence[index]
		"""
		letter = sequence[index]
		if letter == 'U':
			letter = 'C'
		return alphabet.index(letter)

	
	
	def traceback(self):
		"""
		Traceback procedure for computing alignment
		M is scoring matrix, [M,GX,GY]p are pointer matrices

		Returns [align, i, j]
		align[1] is subsequence of x
		align[0] is subsequence of y
		j is index in x where local alignment begins (x[i])
		i is index of y where local alignments begins (y[j])
		"""
		align = ["", ""]
		table = self.Mp
		score = np.max(self.M)
		am = np.argmax(self.M)
		j = am % self.M.shape[1]
		i = am / self.M.shape[1]
		while True:
			if table == self.Mp:
				if table[i][j] == 'GX':
					table = self.GXp
				elif table[i][j] == 'GY':
					table = self.GYp
				align[0] = self.y[i-1] + align[0]
				align[1] = self.x[j-1] + align[1]
				i -= 1
				j -= 1
			elif table == self.GXp:
				if table[i][j] == 'M':
					table = self.Mp
				align[0] = self.y[i-1] + align[0]
				align[1] = "-" + align[1]
				i -= 1
			else:
				if table[i][j] == 'M':
					table = self.Mp
				align[0] = "-" + align[0]
				align[1] = self.x[j-1] + align[1]
				j -= 1
			if table[i][j] == 'STOP':
				break

		self.align_info = [[align[1], align[0]], j, i, score]

	
	@staticmethod
	def printAlignment(alignment):
		"""
		Print single alignment
		"""
		print "SCORE: ", alignment[3]
		print "Index of x:", alignment[1], 
		print "Index of y:", alignment[2]
		print "Alignment: "
		print "x: ",
		print alignment[0][0][0:60]
		print "y: ",
		print alignment[0][1][0:60]




