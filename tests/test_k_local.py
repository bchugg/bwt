import sys
sys.path.append('../')
from kLocalRNA import kLocalFold

seq = "AAAAAAAACCCCCCCCCCCCGGGGGGGGGGGGUUUUUUUU"
print(kLocalFold(seq,2,[0,0,2]))
