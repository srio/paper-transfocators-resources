


import numpy

# SRW thickness 30 um, diameter 1mm
case3_srw = [1.646e15, 2.399e14, 2.312e14, 2.078e14, 2.056e14]
case3_srw_absorbed = []

for i in range(len(case3_srw)):
    if i == 0:
        case3_srw_absorbed.append(0)
    else:
        case3_srw_absorbed.append( 100 * (case3_srw[i-1] - case3_srw[i]) / case3_srw[i-1])

print(case3_srw_absorbed)