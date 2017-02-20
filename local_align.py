import numpy as np
import sys
from params import blosum, score_index, seq, seq_name, d, e

# GLOBAL VARS
min_best_score = -1*sys.maxint

		 
# Parse each protein in file and compute their score
# 	File: in fasta format			 
# 	maximum (optional): maximum number of proteins to compare
def parseAndComputeScore(file, maximum=sys.maxint):
	num_proteins = 0
	# Keep track of best matches: best_matches[i] = [name, index, score]
	best_matches = [[[] for _ in range(7)] for _ in range(3)]
	for i in range(0,3):
		best_matches[i][2] = -1*sys.maxint
	print best_matches
	for line in file:
		if line[0] == '>':
			if num_proteins != 0:
				align_details = computeScoreAndAlignment(amino_acids)
				updateBestMatches(best_matches, align_details, name, index)
			num_proteins += 1
			if num_proteins >= maximum+1:
				break
			amino_acids = ''
			# grab relevant protein info
			index = num_proteins-1
			name = line.split("|")[2].split(" ")[0]
		else: 
			amino_acids += line.strip('\n')

	updateBestMatches(best_matches, computeScoreAndAlignment(amino_acids), name, index)
	printMatches(best_matches)

	

# Determines if current protein is among the top 3 matches so far
#	best_matches is array of top 3 matches so far
#	score, name and index are of current protein
	# alignment = [score, X, Y, X_index, Y_index, length]	
def updateBestMatches(best_matches, alignment, name, index):
	min_index = np.argmin([best_matches[0][2], best_matches[1][2], best_matches[2][2]])
	if best_matches[min_index][2] < alignment[0]:
		best_matches[min_index][0] = name
		best_matches[min_index][1] = index
		best_matches[min_index][2] = alignment[0]
		best_matches[min_index][3] = alignment[3]
		best_matches[min_index][4] = alignment[4]
		best_matches[min_index][5] = alignment[1]
		best_matches[min_index][6] = alignment[2]
	# Update max score
	min_best_score = min(best_matches[0][2], best_matches[1][2], best_matches[2][2])	
	print best_matches[0][0:3],
	print best_matches[1][0:3],
	print best_matches[2][0:3] 

# Print names, scores, indices and alignments of top scores
def printMatches(best_matches):
	for i in range(3):
		print "Index=", best_matches[i][1], 
		print "Name=", best_matches[i][0],
		print "Score=", best_matches[i][2]
		print "START POS in", best_matches[i][0], ":", best_matches[i][3],
		print "START POS in", seq_name, ": ", best_matches[i][4],
		print "LENGTH: ", max(len(best_matches[i][5]), len(best_matches[i][6]))
		print best_matches[i][5]
		print best_matches[i][6]

# Compute alignment using Smith Waterman for local alignment
#	aa is string of amino acids
#   Return array of matching sc
def computeScoreAndAlignment(aa):
	n = len(seq) + 1
	m = len(aa)+1
	# Initialize matrices and pointer matrices
	M, GX, GY = np.zeros((m,n)), np.zeros((m,n)), np.zeros((m,n))
	Mpointer = [ [ [] for _ in range(n) ] for _ in range(m) ]
	GXpointer = [ [ [] for _ in range(n) ] for _ in range(m) ]
	GYpointer = [ [ [] for _ in range(n) ] for _ in range(m) ]

	# Initialization
	GX[0,0], GY[0,0] = -1*sys.maxint, -1*sys.maxint
	Mpointer[0][0], GXpointer[0][0], GYpointer[0][0] = 'STOP', 'STOP', 'STOP'
	for i in range(1,m):
		GY[i,0] = -d
		GX[i,0] = -d - (i-1)*e
		GXpointer[i][0], GYpointer[i][0], Mpointer[i][0] = 'STOP', 'STOP', 'STOP'
	for j in range(1,n):
		GX[0,j] = -1*sys.maxint
		GY[0,j] = -d 
		GXpointer[0][j], GYpointer[0][j], Mpointer[0][j] = 'STOP', 'STOP', 'STOP'
	
	# Compute recurrence
	for i in range(1,m):
		b = findAAIndex(i-1,aa)
		for j in range(1,n):
			a = findAAIndex(j-1,seq)
			# Choices in recurrence
			M_choice = np.append(blosum[a,b]+np.array([M[i-1,j-1], GX[i-1,j-1], GY[i-1,j-1]]),[0])
			GX_choice = [ M[i-1,j]-d, GX[i-1,j]-e ] 
			GY_choice = [ M[i,j-1]-d, GY[i,j-1]-e ] 
			# Max arguments
			M_am = np.argmax(M_choice)
			GX_am = np.argmax(GX_choice)
			GY_am = np.argmax(GY_choice)
			# Compute new values and update pointer matrices
			M[i,j] = max(M_choice)
			if M[i,j] == 0:
				Mpointer[i][j] = 'STOP'
			else:
				Mpointer[i][j] = 'M' if M_am == 0 else ('GX' if M_am == 1 else 'GY')
			GX[i,j] = max(GX_choice)
			GXpointer[i][j] = 'M' if GX_am == 0 else 'GX'
			GY[i,j] = max(GY_choice)
			GYpointer[i][j] = 'M' if GY_am == 0 else 'GY'
	
	score = np.max(M)
	if score > min_best_score:
		align = traceback(M,Mpointer,GXpointer,GYpointer,aa)
		return [score, align[0][0], align[0][1], align[1], align[2]]
	else:
		return [score, "", "", "", ""]

# Traceback procedure for computing an alignment
#	M is score matrix, Mp, GXp, GYp are respective pointer matrices
#   Return 2d-array align, align[0] is subsequence of x, 
#	align[1] is subsequence of y
def traceback(M,Mp,GXp,GYp,aa):
	align = ["", ""]
	table = Mp
	am = np.argmax(M)
	j = am % M.shape[1]
	i = am / M.shape[1]
	while True:
		if table == Mp:
			if table[i][j] == 'GX':
				table = GXp
			elif table[i][j] == 'GY':
				table = GYp
			align[0] = aa[i-1] + align[0]
			align[1] = seq[j-1] + align[1]
			i -= 1
			j -= 1
		elif table == GXp:
			if table[i][j] == 'M':
				table = Mp
			align[0] = aa[i-1] + align[0]
			align[1] = "-" + align[1]
			i -= 1
		else:
			if table[i][j] == 'M':
				table = Mp
			align[0] = "-" + align[0]
			align[1] = seq[j-1] + align[1]
			j -= 1
		if table[i][j] == 'STOP':
			break

	return [align, i+1, j+1]

	
# Find amino acid corresponding to index in given sequence
# If letter is U, switch it for C
def findAAIndex(index, sequence):
	letter = sequence[index]
	if letter == 'U':
		letter = 'C'
	return score_index.index(letter)
	

def main(argv):
	if len(argv) < 2:
		file = open('uniprot.fasta', 'r')
	else: 
		file = open(argv[1], 'r')

	parseAndComputeScore(file)

main(sys.argv)