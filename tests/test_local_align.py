import sys
sys.path.append('../')
from local_align import localAlignRNA as LA

""" 
Test File for local_align module 
""" 

"""
Test 1 - Perfect Complimentary palindromes
"""
seq = "AAAGUCCCGGGACUUU"
align = LA(seq)
align.computeAlignment()
align.printAlignment(align.align_info)

seq = "AUGAGUGACGUAGG"
seq += seq[::-1]
align = LA(seq)
align.computeAlignment()
align.printAlignment(align.align_info)


""" 
Test 2
"""
seq = "CGUAGCUAGGCGGCGAGAGGAUCGAUAUAUAGCCGCGCCCUAUUUUUUU"
align = LA(seq)
align.computeAlignment()
align.printAlignment(align.align_info)





