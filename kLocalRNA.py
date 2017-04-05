import sys
import numpy as np
from disjoint_alignments import disjointAlignments as DA
from r_fold import RNA_Fold as RF
from residuals import residuals
from functools import reduce
import time

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


def kLocalFold(seq,k, rna_fold=1, verbose=1):
	"""
	Carry out k-local-folding and typical rna folding on sequence and 
	compare score. Return k-local score and score from typical RNA folding
	Restrict k-local alignments to have max_knots pseudoknots
	Optional extra parameters: 
	- rna_fold:
		- 1 to do regular rna folding as well
		- 0 to do only k-local folding
	- Verbosity:
		- 0 for no output
		- 1 for printing only scores
		- 2 for obscene verbosity
	"""
	sep = '*'*60 + '\n'

	# Find k optimal alignments between seq and itself. Score 
	# alignments. 
	align = DA(seq)
	res = residuals(seq)

	start_klocal = time.time()
	align.kAlignments(k)
	local_score = DA.scoreAlignments(align.alignments)
	# Pass leftover sequences into RNA folding. 
	inds1 = list(map(lambda x: x[0], align.alignments))
	inds2 = list(map(lambda x: x[1], align.alignments))
	inds = inds1 + inds2

	res.getResiduals(inds)
	res_score = 0
	for r in res.residuals:
		rf = RF()
		rf.fold(r[0])
		res_score += rf.F[0][len(rf.F[0])-1]

	total_score = res_score + local_score
	end_klocal = time.time()

	# Score comparison
	rna_score = 0

	if rna_fold == 1:
		start_rna = time.time()
		rf = RF()
		rf.fold(seq)
		rna_score = rf.F[0][len(rf.F[0])-1]
		end_rna = time.time()
	else:
		start_rna, end_rna, rna_score = 0,0,0

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

	return [total_score, end_klocal-start_klocal, rna_score, end_rna-start_rna]
