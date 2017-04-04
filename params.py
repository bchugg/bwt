import numpy as np

#Parameters file

alphabet = ['C', 'G', 'U', 'A', 'N']
score_matrix = np.matrix([
	[-1,1,-1,-1, 1],
	[1,-1,-1,-1, 1],
	[-1,-1,-1,1, 1],
	[-1,-1,1,-1, 1],
	[1, 1, 1, 1, 1]])


