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

2) Perform any caretaking on leftover subsequences needed 
to pass into RNA folding alg

3) Call RNA Folding on leftover subsequences

4) Append local alignment results and RNA folding results
together. 

5) Output in informative manner. 

"""

seq = "UAUAUAUAUAGCAUGAGUGUCGCGCGCGC"
k = 2
sep = '*'*60 + '\n'

print(sep)
print(' '*12, 'RNA FOLDING USING LOCAL ALIGNMENT\n')

# Find k optimal alignments between seq and itself
align = DA(seq)
align.kAlignments(k)

print(sep)
print('Step 1: K-OPTIMAL LOCAL ALIGNMENTS\n')
align.printAlignments()
local_score = DA.scoreAlignments(align.alignments)
print('RNA FOLDING SCORE:', local_score)
print(sep)


# Pass leftover sequences into RNA folding. 
print('Step 2: PASS UNMATCHED SEQUENCES TO RNA FOLDER\n')
inds1 = list(map(lambda x: x[0], align.alignments))
inds2 = list(map(lambda x: x[1], align.alignments))
inds = inds1 + inds2

res = residuals(seq)
res.getResiduals(inds)
print('Residual subsequences:')
print(res.residuals, '\n')

res_score = 0
for r in res.residuals:
	rf = RF()
	rf.fold(r[0])
	res_score += rf.F[0][len(rf.F[0])-1]

print('Score from RNA folding on residuals:', res_score)
print(sep)

# Score comparison
print('STEP 3: COMPARE SCORES FROM BOTH METHODS\n')
total_score = res_score + local_score
rf = RF()
rf.fold(seq)
rna_score = rf.F[0][len(rf.F[0])-1]

print('Score of this method:', total_score)
print('Score of typical RNA folding', rna_score)












	



