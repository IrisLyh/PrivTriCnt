import numpy as np
from time import *

for dataset in range(3,4):
    for eps in range (1,11):
        read_start = time()
        print("Reading...")
        path_name = "./" + str(dataset) + "_" + str(eps) + "_" + "ptrMat.txt"
        matrix = []
        f=open(path_name,'r')
        lines = f.readlines()
        f.close()
        for l in lines:
            l1 = l.strip('\n')
            l2 = l1.split(' ')
            il = list(map(int,l2))
            matrix.append(il)
        read_end = time()
        read_time = read_end - read_start
        print("Read time:",read_time)
        cal_start_1 = time()
        print("Calculating...")
        noisyMat = np.array(matrix).astype(dtype=int)
        print(noisyMat.shape)
        print(type(noisyMat[0][0]))
        midMat = np.matmul(noisyMat,noisyMat)
        cal_end_1 = time()
        print("Calculation period 1:",(cal_end_1-cal_start_1)/3600)
        finalMat = np.matmul(midMat,noisyMat)
        cal_end_2 = time()
        print("Calculation period 2:",(cal_end_2 - cal_start_2)/3600)
        with open("noisyCnt_"+path_name) as of:
            for user in range(0,len(matrix)):
                print(finalMat[user][user])
