import joblib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
from torch.nn import functional as F
from torch.nn import Module
import torch.optim as optim
from torch.utils.data import Dataset,DataLoader
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler

import pickle
import nni

class PreModel(nn.Module):
    def __init__(self):
        super(PreModel,self).__init__()
        self.Layer1 = nn.Linear(in_features=5,out_features=512)
        self.Layer2 = nn.Linear(in_features=512,out_features=512)
        self.Layer3 = nn.Linear(in_features=512,out_features=4)

    def forward(self,x):
        x = F.relu(self.Layer1(x))
        x = F.relu(self.Layer2(x))
        x = self.Layer3(x)

        return x

premodel = PreModel()
premodel.load_state_dict(torch.load('PBA-ml-pretrain_model_weights.pth'))
for param in premodel.parameters():
    param.requires_grad = False
premodel.eval()

def get_default_parameters():
    params = {
        'lr':0.001,
        'batch_size': 8,
        'units': 512,
        'epochs':300,
        'seed':106
    }
    return params

received_PARAMS = nni.get_next_parameter()
config = get_default_parameters()

excel = r'../../data/experiment.xlsx'
x = pd.read_excel(excel, sheet_name='metals')
y = pd.read_excel(excel, sheet_name='overpotential')
x_train, x_val, y_train, y_val = train_test_split(x, y, test_size=0.2, random_state=int(config['seed']))

norm_x = StandardScaler().fit(x_train)
norm_y = StandardScaler().fit(y_train)
x_train_ = norm_x.transform(x_train)
x_val_ = norm_x.transform(x_val)

y_train_ = norm_y.transform(y_train)
y_val_ = norm_y.transform(y_val)


class GetLoader(Dataset):

    def __init__(self, data_root, data_label):
        self.data = data_root
        self.label = data_label

    def __getitem__(self, index):
        datas = self.data[index]
        labels = self.label[index]
        return (datas, labels)

    def __len__(self):
        return len(self.data)

with torch.no_grad():
    x_train_tensor = torch.Tensor(x_train_)
    x_val_tensor = torch.Tensor(x_val_)
    train_pred_ = premodel(x_train_tensor)
    val_pred_ = premodel(x_val_tensor)

trainDataset = GetLoader(train_pred_,y_train_)
testDataset = GetLoader(val_pred_,y_val_)

train_loader = DataLoader(
    trainDataset,
    batch_size=config['batch_size'],
    shuffle=False
)

test_loader = DataLoader(
    testDataset,
    batch_size=config['batch_size'],
    shuffle=False
)

class Model(Module):
    def __init__(self):
        super(Model,self).__init__()
        self.Layer1 = nn.Linear(in_features=4,out_features=config['units'])
        self.Layer2 = nn.Linear(in_features=config['units'],out_features=config['units'])
        self.Layer3 = nn.Linear(in_features=config['units'],out_features=config['units'])
        self.Layer4 = nn.Linear(in_features=config['units'],out_features=1)

    def forward(self,x):
        x = F.relu(self.Layer1(x))
        x = F.relu(self.Layer2(x))
        x = F.relu(self.Layer3(x))
        x = self.Layer4(x)

        return x

model = Model()
model.load_state_dict(torch.load('train_model_weights.pth'))
model.eval()
with torch.no_grad():
    train_pred = model(train_pred_)
    val_pred = model(val_pred_)

train_pred = train_pred.detach().numpy()
train_pred = norm_y.inverse_transform(train_pred)
val_pred = val_pred.detach().numpy()
val_pred = norm_y.inverse_transform(val_pred)

train_corr = np.corrcoef(y_train.iloc[:, 0], train_pred[:, 0])[0, 1]
val_corr = np.corrcoef(y_val.iloc[:, 0], val_pred[:, 0])[0, 1]
print(train_corr,val_corr)


fig, axes = plt.subplots(1, 2, figsize=(12, 12), tight_layout=True)
axes[0].set_title('Train', fontsize=25)
axes[1].set_title('Test', fontsize=25)

for i in range(y.shape[1]):
    y_min = min(y.iloc[:, i].min(), train_pred[:, i].min(), val_pred[:, i].min())
    y_max = max(y.iloc[:, i].max(), train_pred[:, i].max(), val_pred[:, i].max())
    lim = [y_min - (y_max - y_min) * 0.05, y_max + (y_max - y_min) * 0.05]

    for j in range(2):
        axes[j].set_aspect('equal')
        axes[j].plot(lim, lim, lw=3)
        axes[j].set_xticks(axes[j].get_xticks())
        axes[j].set_yticks(axes[j].get_xticks())
        axes[j].set_xlim(lim)
        axes[j].tick_params(labelsize=25)
        axes[j].set_ylim(lim)

    axes[0].scatter(y_train.iloc[:, i], train_pred[:, i], s=10)
    axes[1].scatter(y_val.iloc[:, i], val_pred[:, i], s=10)

    axes[0].text(y_max, y_min, f'r={train_corr:.3f}', horizontalalignment='right', fontsize=40)
    axes[1].text(y_max, y_min, f'r={val_corr:.3f}', horizontalalignment='right', fontsize=40)

    axes[0].set_ylabel(y.columns[i], fontsize=35)

plt.savefig(f'train.png')

f = open('train_norm_y.pkl', 'wb')
pickle.dump(norm_y, f)
f.close()

f = open('train_norm_x.pkl', 'wb')
pickle.dump(norm_x, f)
f.close()

trian_pred=train_pred.tolist()
val_pred=val_pred.tolist()

train_pred=pd.DataFrame(train_pred,columns=['overpotential'])
val_pred=pd.DataFrame(val_pred,columns=['overpotential'])
with pd.ExcelWriter(f'train_result_all.xlsx') as writer:
    y_train.to_excel(writer, sheet_name='y_train', index=False)
    train_pred.to_excel(writer, sheet_name='train_pred', index=False)
    y_val.to_excel(writer, sheet_name='y_val', index=False)
    val_pred.to_excel(writer, sheet_name='val_pred', index=False)
