from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
RESULT_DIR = SCRIPT_DIR.parent / "results"
DATA_DIR = SCRIPT_DIR.parents[3] / "data"

import numpy as np
import pandas as pd
import pickle
import torch
import torch.nn as nn
from torch.nn import functional as F
from torch.nn import Module

class Model(Module):
    def __init__(self):
        super(Model,self).__init__()
        self.Layer1 = nn.Linear(in_features=5,out_features=512)
        self.Layer2 = nn.Linear(in_features=512,out_features=512)
        self.Layer3 = nn.Linear(in_features=512,out_features=4)

    def forward(self,x):
        x = F.relu(self.Layer1(x))
        x = F.relu(self.Layer2(x))
        x = self.Layer3(x)
        return x

premodel = Model()
premodel.load_state_dict(torch.load(RESULT_DIR / 'PBA-ml-pretrain_model_weights.pth'))
premodel.eval()


class Model(Module):
    def __init__(self):
        super(Model,self).__init__()
        self.Layer1 = nn.Linear(in_features=4,out_features=512)
        self.Layer2 = nn.Linear(in_features=512,out_features=512)
        self.Layer3 = nn.Linear(in_features=512,out_features=512)
        self.Layer4 = nn.Linear(in_features=512,out_features=1)
        self.dropout = nn.Dropout(0.4)

    def forward(self,x):
        x = self.dropout(F.relu(self.Layer1(x)))
        x = self.dropout(F.relu(self.Layer2(x)))
        x = self.dropout(F.relu(self.Layer3(x)))
        x = self.Layer4(x)

        return x

model = Model()
model.load_state_dict(torch.load(RESULT_DIR / 'train_model_weights.pth'))
model.eval()

with open(RESULT_DIR / 'pretrain_norm_x.pkl', 'rb') as f:
    pretrain_norm_x = pickle.load(f)
with open(RESULT_DIR / 'train_norm_y.pkl', 'rb') as f:
    train_norm_y = pickle.load(f)
txt = open(DATA_DIR / "choose_proportion.txt")
lines = txt.readlines()
for i in range(len(lines)):
    lines[i]=lines[i].split()
    if i !=0:
        for j in range(len(lines[i])):
            lines[i][j]=float(lines[i][j])/100

x_val=np.array(lines)
x_val=np.delete(x_val,0,0)
x_val_ = pretrain_norm_x.transform(x_val)
x_val_=torch.Tensor(x_val_)
val_pred_medium = premodel.forward(x_val_)
val_pred_final = model.forward(val_pred_medium)
val_pred_final = val_pred_final.detach().numpy()
val_pred_final = train_norm_y.inverse_transform(val_pred_final)
x_val=pd.DataFrame(x_val,columns=['Mn','Fe','Co','Ni','Zn'])
val_pred_final=pd.DataFrame(val_pred_final,columns=['overpotential'])
with pd.ExcelWriter(RESULT_DIR / 'space_one_percent.xlsx') as writer:
    x_val.to_excel(writer, sheet_name='ratio', index=False)
    val_pred_final.to_excel(writer, sheet_name='value', index=False)
