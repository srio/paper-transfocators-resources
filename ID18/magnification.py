

import numpy
from srxraylib.plot.gol import plot

p1 = 65.0
f1 = 28.2
f2 = 39.7
# D = 188.0 - 65.0
D = 164.0 - 65.0
L = 200.0




q1 = 1/(1/f1 - 1/p1)
p2 = D - q1 # 1/(1/f2 - 1/q2)

q2 = L - p1 - D
# q2 = 1/(1/f2 - 1/p2)

M = (1 - q2 / f2) / (1 - p1 / f1)

print("D: %g, q1+p2: %g" % (D, q1+p2))
print("p1: %g" % p1)
print("q1: %g" % q1)
print("p2: %g" % p2)
print("q2: %g" % q2)
print("L: %g, Sum: %g" % (L, p1+q1+p2+q2))

M1 = q1 / p1
M2 = q2 / p2

print("M: %g, M1: %g, M2: %g,  M1 M2: %g" % (M, M1, M2, M1*M2))

s1 = 70.
print("s1: %g, s1*M: %g, s1*M1: %g" % (s1, s1*M, s1*M1))