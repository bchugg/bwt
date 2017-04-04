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

def testAndPlot(fasta_file, num_sequences, krange, params):
	"""
	Run k-local folding on sequences with for all k in krange.
	Plot results using parameters in params
	- params[0], 1 if figure is to be plotted
	- params[1], filename (optional)

	"""

	fasta_seq = list(SeqIO.parse(open(fasta_file), "fasta"))
	inds = np.floor(np.random.rand(num_sequences)*len(fasta_seq)).astype(int)
	sequences = list(map(lambda x: x.seq[0:200],[fasta_seq[i] for i in inds]))
		
	n = len(sequences)
	averages = []
	times = []
	results = [0,0,0,0]
	for k in krange:
		p = 1 if k == krange[0] else 0
		results[0], results[1] = 0,0
		for seq in sequences:
			fold_results = kLocalFold(seq,k,p,0)
			results[0] += fold_results[0]
			results[1] += fold_results[1]
			if p == 1 : 
				results[2] += fold_results[2]
				results[3] += fold_results[3]

		averages += [[results[0]/n, results[2]/n]]
		times += [[results[1]/n, results[3]/n]]

	title = 'TIME_k'+str(krange[0])+'-'+str(krange[len(krange)-1])
	title += '_'+str(n)+'seqs_'+params[1]
	k_local_ave = list(map(lambda x: x[0], averages))
	rna_ave = list(map(lambda x: x[1], averages))
	plot_scores(krange, k_local_ave,rna_ave, title, params)
	time_ratio = list(map(lambda x: x[0]/float(x[1]), times))
	plot_times(krange, time_ratio, title, params)

def plot_scores(krange, k_local_ave, rna_ave, title, params):
	plt.plot(krange, k_local_ave, 'go', label = 'k-Local Folding')
	plt.plot(krange, rna_ave, 'b--', label = 'Standard Folding')
	plt.ylabel('Average Scores')
	plt.xlabel('k')
	plt.title('Scores of k-Local Folding vs Nussinov')
	plt.legend()
	if len(params)>1:
		name = 'SCORE_' + title
		path = 'figures/' + name
		plt.savefig(path)
	if params[0]:
		plt.show()
	plt.clf()
	

def plot_times(krange, time_ratio, title, params):
	plt.plot(krange, time_ratio)
	plt.ylabel('Time Ratio: k-Local to Nussinov')
	plt.xlabel('k')
	plt.title('Time Comparison of k-Local Folding to Nussinov')
	if len(params)>1:
		name = 'TIME_' + title
		path = 'figures/' + name
		plt.savefig(path)
	if params[0]:
		plt.show()
	plt.clf()
		

def main():
	"""
	Call tests for k-local folding.

	"""
	#Parameters:
	files = ["data/rna.fasta", 'data/5sRNA.fasta', 'data/ciliateRna.fasta', 'data/vRna.fasta']
	num_seqs = 20
	k_ranges = [[0,1,2,3,4,5,6,7,8], [15,16,17,18,19]]
	plot = 0 # Change to 1 to have results plotted each iteration
	save = 1
	
	for f in files:
		for r in k_ranges:
			filename = f.split('data/')[1].split('.fasta')[0]
			if save == 1:
				params = [plot,filename]
			else: 
				params = [plot]
			testAndPlot(f, num_seqs, r, params)



main()
