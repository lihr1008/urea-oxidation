# -*- coding: utf-8 -*-

import os
import re
import sys

import numpy as np
import pandas as pd


def get_atom_dict(data):
    with open(data) as file:
        lines = file.read()
    atom = np.array(
        lines[re.search('Masses', lines).span()[1]:re.search('Bond Coeffs', lines).span()[0]].split()).reshape(-1, 4)
    num = atom[:, 0]
    atom = [re.findall('[a-zA-Z]*', i[-1])[0] for i in atom]
    atom = dict(zip(num, atom))
    return atom


def get_distance(a, b, size):
    distance = np.abs(a - b)
    distance = np.min([distance, size - distance], axis=0)
    distance = (distance ** 2).sum(axis=-1) ** 0.5
    return distance


def get_structs(direct, group):
    print(direct)
    if not os.path.exists(f'{direct}/movie.data'):
        metal_num.append([])
        group_num.append([])
        return
    atom_dict = get_atom_dict(f'{direct}/data.box')

    with open(f'{direct}/movie.data') as file:
        lines = file.readlines()
    group = {i: 0 for i in group}
    num = 580
    elements = []
    flag = 0
    for i, line in enumerate(lines):
        if line == 'ITEM: ATOMS id type x y z\n':
            time = int(lines[i - 7])
            if time < 10000:
                continue
            if time == 1000000:
                flag = 1

            size = lines[i - 3].split()
            size = float(size[1]) - float(size[0])

            atoms = lines[i + 1:i + num + 1]
            metals, elements = [], []
            for atom in atoms:
                atom = atom.split()
                if atom[1]== '8':
                    continue
                else:
                    element = atom_dict[atom[1]]
                    if element in Metal:
                        elements.append(element)
                        metals.append(atom[-3:])

            metals = np.array(metals).astype(float)
            elements = np.array(elements)

            for index,metal in enumerate(metals):
                distances=get_distance(metal, metals, size)
                distances=distances.tolist()
                for id,dis in enumerate(distances):
                    distances[id]=[dis,id]
                distances=np.array(distances)
                neighbor = distances[np.argsort(distances[:, 0])]
                group[elements[int(neighbor[0][1])] + elements[int(neighbor[1][1])] + elements[int(neighbor[2][1])] + elements[int(neighbor[3][1])]+elements[int(neighbor[4][1])]+elements[int(neighbor[5][1])]]+=1

    if flag == 0:
        metal_num.append([])
        group_num.append([])
        return

    metal_num.append([sum(elements == m) for m in Metal])
    group_num.append(group.values())


    return


Metal = ['Mn', 'Fe', 'Co', 'Ni', 'Zn']

group = []
for i in range(5):
    for j in range(5):
        for k in range(5):
            for m in range(5):
                for n in range(5):
                    for l in range(5):
                        group.append(Metal[i] + Metal[j] + Metal[k]+Metal[m]+Metal[n]+Metal[l])

metal_num = []
group_num = []
add=[]
start=1
stop=2
for i in range(start, stop+1):
    get_structs(f'data example/{i}',group)
metal_num = pd.DataFrame(metal_num, columns=Metal)
group_num = pd.DataFrame(group_num, columns=group)

for key,value in group_num.items():
    for num in value.values:
        if num != 0:
            add.append(key)
            break

for key, value in group_num.items():
    if key not in add:
        group_num.pop(key)

with pd.ExcelWriter(f'results/output_{start}-{stop}.xlsx') as writer:
    metal_num.to_excel(writer, sheet_name='metals', index=False)
    group_num.to_excel(writer, sheet_name='groups', index=False)
