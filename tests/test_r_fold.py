import sys
sys.path.append('../')
from r_fold import RNA_Fold as RF

seq = "GGGAAAUCC"

rf = RF()

print("case1: " + str(rf.fold(seq)))

seq = "CGUCACAGUAGAAGCCUAUAUGUGAGGCUACUCGGC"

print("case2: " + str(rf.fold(seq)))
rf.print_F()




