import sys 
sys.path.append('../')
from disjoint_alignments import disjointAlignments as DA

"""
Tests for disjoint_alignments module 
"""
sep = "*"*50 + '\n'

"""
Test 1
Random
"""
print ('TEST 1\n')
seq = "CGUAGCUAGGCGGCGAGAGGAUCGAUAUAUAGCCGCGCCCUAUUUUUUU"
kaligns = DA(seq)
kaligns.kAlignments(10)
kaligns.printAlignments()
print(sep)

"""
Test 2
Forced hairpin at ends
"""
print('TEST 2\n')
seq = "CCCCCCCCCCAAAAAAAAAAACCCCGGGGAAAAAAAAAAGGGGGGGGGG"
kaligns = DA(seq)
kaligns.kAlignments(2)
kaligns.printAlignments()

