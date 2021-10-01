


import numpy

# SRW before slit, after slit, after TF1, after TF2
case1_srw = [1.644e15, 3.966e13, 3.665e13, 1.794e13]
case2_srw = [1.644e15, 3.966e13, 3.703e13, 3.566e13]
case3_srw = [1.644e15, 1.605e14, 1.507e14, 1.177e14]
case4_srw = [1.644e15, 1.605e14, 1.481e14, 1.427e14]



for case in range(4):
    case_srw_absorbed = []
    if case == 0:
        case_srw = case1_srw
    elif case == 1:
        case_srw = case2_srw
    elif case == 2:
        case_srw = case3_srw
    elif case == 3:
        case_srw = case4_srw

    for i in range(4):
        if i == 0:
            case_srw_absorbed.append(0)
        else:
            case_srw_absorbed.append( 100 * (case_srw[i-1] - case_srw[i]) / case_srw[i-1])

    print(case_srw_absorbed)