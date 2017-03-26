import sys
sys.path.append('../')
from local_align import disjointAlignments as DA
from local_align import localAlign as LA

""" 
Test File for local_align module
""" 

"""
Tests for localAlign class
"""

""" 
Perfect matches
"""
seq1 = "AAAARRRRNNNNNDDDDDD"
seq2 = "AAAARRRRNNNNNDDDDDD"

align = DA(seq1, seq2)
align.kAlignments(2)
#align.printAlignments()


"""
Matches with gaps
"""
seq1 = "AAAADDDDDDRRRRDDDDDNNNNNFFFFFFDDDDDD"
seq2 = "FFFFFFRRRRFFFFFNNNNNDDDDDDDDDDDDAAAAA"
align = DA(seq1, seq2)
align.kAlignments(7)
align.printAlignments()



