import sys
import numpy as np
from disjoint_alignments import disjointAlignments as DA
from r_fold import RNA_Fold as RF
from residuals import residuals
from functools import reduce

"""
Main File for running local alignment RNA folding

Written jointly by:
Coulter Beeson, Ben Chugg, Kenny Drabble, Jeffrey Jeyachandren


STEPS:

1) Call k-Optimal Local Alignment with specific k

2) Call RNA Folding on leftover subsequences

3) Compare scores between k-local-folding and typical 
folding

"""


def kLocalFold(seq,k, params=[np.inf, 0,0]):
	"""
	Carry out k-local-folding and typical rna folding on sequence and 
	compare score. Return k-local score and score from typical RNA folding
	Restrict k-local alignments to have max_knots pseudoknots
	Optional extra parameters: 
	- params[0] = maximum number of pseudoknots allowed
	- params[1] = 
		- 1 if only want to perform k local folding, 
		- perform typical rna folding as well otherwise
	- params[2] = verbosity
		- 0 for no output
		- 1 for printing only scores
		- 2 for obscene verbosity
	"""
	sep = '*'*60 + '\n'

	# Find k optimal alignments between seq and itself. Score 
	# alignments. 
	align = DA(seq)
	align.kAlignments(k)
	try:
		if params[0] < np.inf: align.removePseudoknots(params[0])
	except IndexError:
		print('Not removing pseudoknots')

	local_score = DA.scoreAlignments(align.alignments)
	
	# Pass leftover sequences into RNA folding. 
	inds1 = list(map(lambda x: x[0], align.alignments))
	inds2 = list(map(lambda x: x[1], align.alignments))
	inds = inds1 + inds2

	res = residuals(seq)
	res.getResiduals(inds)
	res_score = 0
	for r in res.residuals:
		rf = RF()
		rf.fold(r[0])
		res_score += rf.F[0][len(rf.F[0])-1]

	# Score comparison
	total_score = res_score + local_score
	rna_score = 0
	try: 
		if params[1] != 0:
			rf = RF()
			rf.fold(seq)
			rna_score = rf.F[0][len(rf.F[0])-1]
	except IndexError:
		print("Performing regular RNA folding as comparison")
	
	try: 
		verbose = params[2]
	except IndexError:
		verbose = 0

	if verbose == 2:
		print(sep)
		print(' '*12, 'RNA FOLDING USING LOCAL ALIGNMENT\n')
		print(sep)
		print('Step 1: K-OPTIMAL LOCAL ALIGNMENTS\n')
		align.printAlignments()
		print('RNA FOLDING SCORE:', local_score)
		print(sep)
		print('Step 2: PASS UNMATCHED SEQUENCES TO RNA FOLDER\n')
		print('Residual subsequences:')
		print(res.residuals, '\n')
		print('Score from RNA folding on residuals:', res_score)
		print(sep)
		print('STEP 3: COMPARE SCORES FROM BOTH METHODS\n')
	
	if verbose >= 1:	
		print('Score of this method:        ', total_score)
		print('Score of typical RNA folding:', rna_score, '\n')

	return[total_score, rna_score]
