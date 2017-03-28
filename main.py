# Dependencies
import sys
import numpy as np
from local_align import disjointAlignments as DA
from local_align import localAlign as LA
from r_fold import RNA_Fold as RF


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

# Find k optimal alignments between seq and itself
rev = seq[::-1]
align = DA(seq,rev)
align.kAlignments(k)

print(sep)
print('Step 1: K-OPTIMAL LOCAL ALIGNMENTS\n')
align.printAlignments()
print(sep)

# Pass leftover sequences into RNA folding. 
print('Step 2: PASS UNMATCHED SEQUENCES TO RNA FOLDER\n')

leftovers = [[]]
newseq = seq
cumul_ind = 0
for i in range(align.k):
	cut = align.alignments[i][0][0]
	if len(cut) > 0:
		ind = newseq.find(cut)
		leftovers += [[newseq[:ind], cumul_ind]]
		cumul_ind += ind + len(cut)
		newseq = newseq[ind+len(cut):]
if len(newseq) != 0:
	leftovers = np.append(leftovers, [])

if np.size(leftovers) >0:
	print('Leftover Sequences to be passed into RNA Folding:')
	for i in range(1,np.size(leftovers)):
		s = leftovers[i][0]
		l = leftovers[i][1]
		print('Sequence:')
		print(s)
		print('Index: ', l, '-', l+len(s))
else:
	print('No unmatched regions found.')


# Run RNA folding on Leftovers 
rf = RF()
for i in range(1,np.size(leftovers)):
	if len(leftovers[i][0]) != 0:
		rf.fold(leftovers[i][0])
		rf.print_F()






	



