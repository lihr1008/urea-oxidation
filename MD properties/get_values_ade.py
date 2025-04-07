import numpy as np
import pandas as pd
import sys
import re

# start, stop = int(sys.argv[1]), int(sys.argv[2])
start, stop = 1,2
excel=f'../MD preprocess/results/output_{start}-{stop}.xlsx'
all = pd.read_excel(excel,sheet_name='groups')
all_num=all.values

ade_path=r'data/ade_final.xlsx'
ade_data = pd.read_excel(ade_path)
ade_data=ade_data.values
metals=ade_data[:,:6]
e_urea=ade_data[:,-2]
e_CO=ade_data[:,-1]


names=all.keys()
values=names.values

final_ade_urea=[]
final_ade_CO=[]

for num in range(all_num.shape[0]):
    one_ade_urea=0
    one_ade_CO=0
    for i in range(len(values)):
        name=values[i]
        six=[]
        m = 0
        for j in range(6):
            six.append(name[m:m+2])
            m+=2

        five=six[1:]
        five.sort()
        ture_six=[six[0]]
        for k in five:
            ture_six.append(k)

        for n in range(len(ade_data)):
            a = list(ade_data[n][1:6])
            a.sort()
            ture_a = [ade_data[n][0]]
            for k in a:
                ture_a.append(k)

            if ture_six == ture_a and all_num[num][i] !=0 :
                one_ade_urea=one_ade_urea+e_urea[n]*all_num[num][i]
                one_ade_CO= one_ade_CO + e_CO[n] * all_num[num][i]

    final_ade_urea.append(one_ade_urea/6000)
    final_ade_CO.append(one_ade_CO/6000)

print(final_ade_urea,final_ade_CO)
with open(f'results/ade/{start}-{stop}-ade.txt','w+') as f:
    for i in range(start,stop+1):
        f.write(str(i)+'\t'+str(final_ade_urea[i-start])+'\t'+str(final_ade_CO[i-start])+'\n')

