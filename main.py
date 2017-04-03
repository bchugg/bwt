import sys
import numpy as np
from kLocalRNA import kLocalFold
from Bio import SeqIO
from random import random
from math import ceil
import matplotlib.pyplot as plt
from functools import reduce


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

def testAndPlot(sequences, krange, params):
	"""
	Run k-local folding on sequences with for all k in krange. 
	Plot results using parameters in params 
	- params[0], ylabel
	- params[1], xlabel
	- params[2], title
	- params[3], filename (optional)
	"""
	n = len(sequences)
	averages = []
	for k in krange:
		k_fold_average, rna_fold_average = 0,0
		for seq in sequences:
			scores = kLocalFold(seq,k,0)
			k_fold_average += scores[0]
			rna_fold_average += scores[1]
		k_fold_average = k_fold_average/n
		rna_fold_average = rna_fold_average/n
		averages += [[k_fold_average,rna_fold_average]]

	k_local_ave = list(map(lambda x: x[0], averages))
	rna_ave = list(map(lambda x: x[1], averages))
	plt.plot(krange, k_local_ave, 'go', label = 'k-Local Folding')
	plt.plot(krange, rna_ave, 'bo', label ='Standard Folding')
	plt.ylabel(params[0])
	plt.xlabel(params[1])
	plt.title(params[2])
	plt.legend()
	if len(params) >= 4:
		path = 'figures/' + params[3]
		plt.savefig(path)
	plt.show()


def main():
	"""
	Call tests for k-local folding.

	Parameters:
	"""
	num_sequences = 10
	krange = [19,20,21,22,23,24,25]
	title = 'k-Local Folding vs Standard Folding'
	filename = 'k1-10_20seqs'
	params = ["Average Score", "k", title, filename]

	fasta_seq = list(SeqIO.parse(open("rna.fasta"), "fasta"))
	inds = np.floor(np.random.rand(num_sequences)*len(fasta_seq)).astype(int)
	sequences = list(map(lambda x: x.seq[0:250],[fasta_seq[i] for i in inds]))
	
	
	testAndPlot(sequences, krange, params)


main()






	



