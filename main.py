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

def testAndPlot(fasta_file, num_sequences, max_knots, krange, params):
	"""
	Run k-local folding on sequences with for all k in krange. 
	Plot results using parameters in params 
	- params[0], ylabel
	- params[1], xlabel
	- params[2], title
	- params[3], 1 if figure is to be plotted
	- params[4], filename (optional)

	"""

	fasta_seq = list(SeqIO.parse(open(fasta_file), "fasta"))
	inds = np.floor(np.random.rand(num_sequences)*len(fasta_seq)).astype(int)
	sequences = list(map(lambda x: x.seq[0:200],[fasta_seq[i] for i in inds]))
		
	n = len(sequences)
	averages = []
	rna_fold_average = 0
	for k in krange:
		k_fold_average = 0
		p = 1 if k == krange[0] else 0 
		for seq in sequences:
			scores = kLocalFold(seq,k,[max_knots,p,0])
			k_fold_average += scores[0]
			if p == 1 : rna_fold_average += scores[1]
		averages += [[k_fold_average/n, rna_fold_average/n]]

	k_local_ave = list(map(lambda x: x[0], averages))
	rna_ave = list(map(lambda x: x[1], averages))
	plt.plot(krange, k_local_ave, 'go', label = 'k-Local Folding')
	plt.plot(krange, rna_ave, 'b--', label ='Standard Folding')
	plt.ylabel(params[0])
	plt.xlabel(params[1])
	plt.title(params[2])
	plt.legend()
	if len(params) > 4:
		path = 'figures/' + params[4]
		plt.savefig(path)
	if params[3] == 1:
		plt.show()
	plt.clf()
		

def main():
	"""
	Call tests for k-local folding.

	"""
	#Parameters:
	files = ["data/rna.fasta", 'data/5sRNA.fasta', 'data/ciliate.fasta', 'vRna.fasta']
	num_seqs = 20
	k_ranges = [[0,1,2,3,4,5,6,7,8], [15,16,17,18,19,20]]
	max_knots_range = [0,3,5,7]
	global_title = 'k-Local Folding vs Standard Folding'
	plot = 0 # Change to 1 to have results plotted each iteration
	
	for f in files:
		for r in k_ranges:
			for m in max_knots_range:
				filename = f.split('data/')[1].split('.fasta')[0]
				title = global_title+' with '+str(num_seqs)+' sequences'
				saveas = 'k'+str(r[0])+'_-'+str(r[len(r)-1])+'_'+str(num_seqs)+'seqs_'
				saveas += str(m)+'knots_250length_'+filename
				params = ["Average Score", "k", title, plot, saveas]
				testAndPlot(f, num_seqs, m, r, params)


main()






	



