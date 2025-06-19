import numpy as np
from sklearn.preprocessing import StandardScaler
import pandas as pd
from sklearn.model_selection import train_test_split
import torch
import torch.nn as nn
from torch.nn import functional as F
from torch.nn import Module
import torch.optim as optim

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

model = torch.load('PBA-ml-pretrain_model.pkl')
state_dict = model.state_dict()
torch.save(state_dict, 'PBA-ml-pretrain_model_weights.pth')
premodel = Model()
premodel.load_state_dict(torch.load('PBA-ml-pretrain_model_weights.pth'))
premodel.eval()


class Model(Module):
    def __init__(self):
        super(Model,self).__init__()
        self.Layer1 = nn.Linear(in_features=4,out_features=512)
        self.Layer2 = nn.Linear(in_features=512,out_features=512)
        self.Layer3 = nn.Linear(in_features=512,out_features=512)
        self.Layer4 = nn.Linear(in_features=512,out_features=1)

    def forward(self,x):
        x = F.relu(self.Layer1(x))
        x = F.relu(self.Layer2(x))
        x = F.relu(self.Layer3(x))
        x = self.Layer4(x)

        return x

trained_model = torch.load('PBA-ml-train_model.pth')
state_dict = trained_model.state_dict()
model = Model()
torch.save(state_dict, 'train_model_weights.pth')
model.load_state_dict(torch.load('train_model_weights.pth'))
model.eval()

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
