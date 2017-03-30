import sys
import numpy as np
from kLocalRNA import kLocalFold
from Bio import SeqIO
from random import random
from math import ceil


"""
Main File calling and testing k-local RNA folding. 

"""

def randomRNA(length):
	"""
	Generate Random RNA sequences drawn from uniform 
	distribution
	"""
	alph = ['A','C','U','G']
	seq = ""
	for i in range(length):
		r = int(ceil(4*random()))
		seq += alph[r-1]
	return seq


def main():
	fasta_seq = SeqIO.parse(open("rna.fasta"), "fasta")
	for fasta in fasta_seq:
		kLocalFold(fasta.seq[0:250],10,0)
	#for i in range(10):
	# 	kLocalFold(randomRNA(200), 5,0)


main()






	



