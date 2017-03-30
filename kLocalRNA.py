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


def kLocalFold(seq,k, verbose=1):
	"""
	Carry out k-local-folding and typical rna folding on sequence and 
	compare score. Return k-local score and score from typical RNA folding
	Print progress if verbose.
	"""
	sep = '*'*60 + '\n'

	# Find k optimal alignments between seq and itself. Score 
	# alignments. 
	align = DA(seq)
	align.kAlignments(k)
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
	rf = RF()
	rf.fold(seq)
	rna_score = rf.F[0][len(rf.F[0])-1]

	
	if verbose:
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
		
	print('Score of this method:        ', total_score)
	print('Score of typical RNA folding:', rna_score, '\n')

	return[total_score, rna_score]
