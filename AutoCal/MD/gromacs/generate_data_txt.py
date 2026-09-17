import sys

import numpy as np

# ['Mn', 'Fe', 'Co', 'Ni', 'Zn', 'CN']

total = int(sys.argv[1])  # python generate_data_txt.py 20000
filename = 'data.txt'
lines = []
rng = np.random.RandomState(0)

for i in range(total):

    div = rng.randint(0, 56, 4)
    div = np.sort(div)

    Mn = div[0] + 1
    Fe = div[1] - div[0] + 1
    Co = div[2] - div[1] + 1
    Ni = div[3] - div[2] + 1
    Zn = 56 - div[3]
    CN = 40

    line = [str(x) for x in [Mn, Fe, Co, Ni, Zn, CN]]
    lines.append('\t'.join(line) + '\n')

with open(filename, 'a') as file:
    file.writelines(lines)
