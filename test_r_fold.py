from r_fold import RNA_Fold as RF
from random import random
from math import *
import sys


rf = RF()


def gen_rand_seq(n):

	alph = ['A','C','U','G']
	seq = ""
	for i in xrange(n):
		r = int(ceil(4*random()))
		seq += alph[r-1]
	return seq


#-------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------

#sys.setrecursionlimit(1200)

#test 1
#Expect 3
seq = "GGGAAAUCC"
print("case1: " + str(rf.fold(seq)))

#test2
#Expect 10 
seq = "GGGUUUGGGCCCUUUGGUGGCCCCUUUCCC"
print("case2: " + str(rf.fold(seq)))


#test3
#Expect 16
seq = "CGUCACAGUAGAAGCCUAUAUGUGAGGCUACUCGGC"
print("case3: " + str(rf.fold(seq)))


#test4
#Expect ?
seq = "ACACGACCUCAUAUAAUCUUGGGAAUAUGGCCCAUAAGUUUCUACCCGGCAACCGUAAAUUGCCGGACUAUGCAGGGAAGUGA"
print("case4: " + str(rf.fold(seq)))

#test5
#Expect random 
seq = gen_rand_seq(200)
print(seq)
print("case5: " + str(rf.fold(seq)))

#test6
#Expect random 
seq = gen_rand_seq(200)
print(seq)
print("case6: " + str(rf.fold(seq)))

#test7
#Expect random 
seq = gen_rand_seq(200)
print(seq)
print("case7: " + str(rf.fold(seq)))

#test8
#Expect random 
seq = gen_rand_seq(200)
print(seq)
print("case8: " + str(rf.fold(seq)))






