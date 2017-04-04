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

"""
Test 2
Forced hairpin at ends
"""
print('TEST 2\n')
seq = "CCCCCCCCCCAAAAAAAAAAACCCCGGGGAAAAAAAAAAGGGGGGGGGG"
kaligns = DA(seq)
kaligns.kAlignments(2)
kaligns.printAlignments()

"""
Test 3
Test scoring 
Expected: 6
"""
print('TEST 3')
alignments = [[0,0,0,['AAAAAA', 'UUUUUU']]]
print(DA.scoreAlignments(alignments), 'Expected:', 6)

"""
Test 4
Test scoring 
Expected: 9
"""
print('TEST 4')
alignments = [[0,0,0,['CGAACCCGUUA', 'GCUUGAACAAU']]]
print(DA.scoreAlignments(alignments), 'Expected:', 9)

"""
Test 5
Test scoring 
Expected: 11
"""
print('Test 5')
alignments = [[0,0,0,['CGAACGUUA', 'GCUUGCAAU']], [0,0,0,['AACG','UUCG']]]
print(DA.scoreAlignments(alignments), 'Expected:', 11)
print('\n')


"""
Test 6
Removing Pseudoknots
"""

print('Test 6')
align1 = [[3,6], [12,15], 20]
align2 = [[6,8], [16,18], 10]
kaligns = DA(seq)
kaligns.alignments = [align1, align2]
kaligns.k = 2
kaligns.removePseudoknots(0)
print(kaligns.alignments, 'k=', kaligns.k)

"""
Test 6
Removing Pseudoknots
"""

print('Test 6')
align1 = [[3,6], [12,15], 20]
align2 = [[6,8], [16,18], 10]
align3 = [[8,10], [18,19], 15]
kaligns = DA(seq)
kaligns.alignments = [align1, align2, align3]
kaligns.k = 2
kaligns.removePseudoknots(1)
print(kaligns.alignments, 'k=', kaligns.k)














