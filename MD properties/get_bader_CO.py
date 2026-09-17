import numpy as np
import pandas as pd
import sys
import re

# start, stop = int(sys.argv[1]), int(sys.argv[2])
start, stop = 1,2
excel=f'../MD preprocess/results/output_{start}-{stop}.xlsx'
all = pd.read_excel(excel,sheet_name='groups')
all_num=all.values

bader_path=r'data/bader_CO_final.xlsx'
bader_data = pd.read_excel(bader_path)
bader_data=bader_data.values
metals=bader_data[:,:6]
bader_CO=bader_data[:,-1]

names=all.keys()
values=names.values

final_bader_CO=[]

for num in range(all_num.shape[0]):
    one_bader_CO = 0
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

        for n in range(len(bader_data)):
            a = list(bader_data[n][1:6])
            a.sort()
            ture_a = [bader_data[n][0]]
            for k in a:
                ture_a.append(k)

            if ture_six == ture_a and all_num[num][i] !=0 :
                one_bader_CO= one_bader_CO + bader_CO[n] * all_num[num][i]

    final_bader_CO.append(one_bader_CO/6000)

with open(f'results/bader/processed_bader_CO_{start}-{stop}.txt','w+') as f:
    for i in range(start,stop+1):
        f.write(str(i)+'\t'+str(final_bader_CO[i-start])+'\n')


