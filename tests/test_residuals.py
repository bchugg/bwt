import sys
sys.path.append('../')
from residuals import residuals


"""
Test class for residuals module
"""


"""
Test 1
"""
print('Test 1')
# Should return only Y's
x = "XXXXXYYYYXXXXYYYYYY"
inds = [[0,5], [9,13]]
res = residuals(x)
res.getResiduals(inds)
print(res.residuals)

inds = [[0,4], [0,5], [1,3], [9,9], [10,10], [11,13]]
res = residuals(x)
res.getResiduals(inds)
print(res.residuals)

res = residuals(x)
res.getResiduals(inds[::-1])
print(res.residuals)


"""
Test 2
"""
print('Test 2')
x = "YYYYXYYYY"
inds = [[4,4], [4,4]]
res = residuals(x)
res.getResiduals(inds)
print(res.residuals)

"""
Test 3
"""
print('Test 3')
x = "YXXXXXXXXY"
inds = [[1,1], [1,7], [7,9]]
res = residuals(x)
res.getResiduals(inds)
print(res.residuals)

"""
Test 4
"""
print('Test 4')
x = "XXXXXYYYYYYYXXXXXXXXXX"
inds = [[0,5], [10,len(x)]]
res = residuals(x)
res.getResiduals(inds)
print(res.residuals)



