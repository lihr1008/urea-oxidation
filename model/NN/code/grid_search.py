from train import premodel,model
import numpy as np
from sklearn.preprocessing import StandardScaler
import pandas as pd
from sklearn.model_selection import train_test_split
import torch

excel = '../../data/experiment.xlsx'
x = pd.read_excel(excel, sheet_name='metals')
y = pd.read_excel(excel, sheet_name='overpotential')
x_train, x_val, y_train, y_val = train_test_split(x, y, test_size=0.2,random_state=106)


norm_x = StandardScaler().fit(x_train)
norm_y = StandardScaler().fit(y_train)
txt = open(f'../../data/choose_proportion.txt')
lines = txt.readlines()
for i in range(len(lines)):
    lines[i]=lines[i].split()
    if i !=0:
        for j in range(len(lines[i])):
            lines[i][j]=float(lines[i][j])/100

x_val=np.array(lines)
x_val=np.delete(x_val,0,0)
norm_x = StandardScaler().fit(x_train)
norm_y = StandardScaler().fit(y_train)
x_val_ = norm_x.transform(x_val)
x_val_=torch.Tensor(x_val_)
val_pred_medium = premodel.forward(x_val_)
val_pred_final = model.forward(val_pred_medium)
val_pred_final = val_pred_final.detach().numpy()
val_pred_final = norm_y.inverse_transform(val_pred_final)
x_val=pd.DataFrame(x_val,columns=['Mn','Fe','Co','Ni','Zn'])
val_pred_final=pd.DataFrame(val_pred_final,columns=['overpotential'])
with pd.ExcelWriter(f'../results/space_one_percent.xlsx') as writer:
    x_val.to_excel(writer, sheet_name='ratio', index=False)
    val_pred_final.to_excel(writer, sheet_name='value', index=False)
