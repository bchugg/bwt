# Dependencies
import sys
from r_fold.py import RNA_FOLD as RNA_FOLD
from local_align.py import localAlign as LA

"""
Main File for running BWT_RNA

Written jointly by:
Coulter Beeson, Ben Chugg, Kenny Drabble, Jeffrey Jeyachandren

* Run Instructions *

"""

"""
STEPS:

1) Call k-Optimal Local Alignment with specific k

2) Perform any caretaking on leftover subsequences needed 
to pass into RNA folding alg

3) Call RNA Folding on leftover subsequences

4) Append local alignment results and RNA folding results
together. 

5) Output in informative manner. 


"""