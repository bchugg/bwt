import sys
sys.path.append('../')
from kLocalRNA import kLocalFold

seq = "AAAAAAAACCCCCCCCCCCCGGGGGGGGGGGGUUUUUUUU"
print(kLocalFold(seq,0))
