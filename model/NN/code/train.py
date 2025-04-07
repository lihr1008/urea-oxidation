# -*- coding: UTF-8 -*-
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
from pretrain import premodel
import pickle

excel = '../../data/experiment.xlsx'
x = pd.read_excel(excel, sheet_name='metals')
y = pd.read_excel(excel, sheet_name='overpotential')
x_train, x_val, y_train, y_val = train_test_split(x, y, test_size=0.2, random_state=106)

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


x_train_=torch.Tensor(x_train_)
x_val_=torch.Tensor(x_val_)
train_pred_ = premodel.forward(x_train_)
val_pred_ = premodel.forward(x_val_)

trainDataset = GetLoader(train_pred_,y_train_)
testDataset = GetLoader(val_pred_,y_val_)

train_loader = DataLoader(
    trainDataset,
    batch_size=8,
    shuffle=False
)

test_loader = DataLoader(
    testDataset,
    batch_size=8,
    shuffle=False
)
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

model = Model()

optimizer = optim.Adam(model.parameters(),lr=0.001)
scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(optimizer,'min',factor=0.3, min_lr=1e-5,patience=10)
criteon = nn.MSELoss()
epochs=300
for epoch in range(epochs):
    temp_loss = []
    for batch_idx,train_data in enumerate(train_loader):
        # model.train()
        (data, label) = train_data
        # data=data.to(torch.float32)
        label = label.to(torch.float32)
        logits = model(data)
        #logits = model.train(data)
        loss = criteon(logits,label)
        temp_loss.append(loss.data)
        optimizer.zero_grad()
        loss.backward(retain_graph=True)
        optimizer.step()

    temp_loss = np.array(temp_loss)
    train_loss = np.average(temp_loss)
    print('Train Epoch : {}\tLoss:{}'.format(
                epoch,train_loss
        ))

    temp_loss = []
    test_loss = 0.0
    correct = 0
    for index,test_data in enumerate(test_loader):
        model.eval()
        (data,label) = test_data
        data = data.to(torch.float32)
        label = label.to(torch.float32)
        logits = model(data)
        loss = criteon(logits, label)
        temp_loss.append(loss.data)
        test_loss += loss.item()
        pred = logits.data.argmax(0)
        correct += pred.eq(label.data.argmax(0)).sum().item()

    temp_loss = np.array(temp_loss)
    val_loss = np.average(temp_loss)
    scheduler.step(val_loss)

    print('\nTest epoch:{}, Average loss: {}\n'.format(
        epoch,val_loss))


print("完成")

train_pred = model.forward(train_pred_)
#train_pred = model.train(train_pred_)
val_pred = model.forward(val_pred_)
#val_pred = model.eval(val_pred_)

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

plt.savefig(f'../results/train.png')

torch.save(model, '../results/PBA-ml-train_model.pkl')


class FullModel(nn.Module):
    def __init__(self, premodel, formal_model):
        super().__init__()
        self.premodel = premodel
        self.formal_model = formal_model
    def forward(self, x):
        x = self.premodel(x)
        return self.formal_model(x)

full_model = FullModel(premodel, model)
torch.save(full_model, '../results/full_model.pkl')

#保存
f = open('../results/train_norm_y.pckl', 'wb')
pickle.dump(norm_y, f)
f.close()

f = open('../results/train_norm_x.pckl', 'wb')
pickle.dump(norm_x, f)
f.close()

trian_pred=train_pred.tolist()
val_pred=val_pred.tolist()

train_pred=pd.DataFrame(train_pred,columns=['overpotential'])
val_pred=pd.DataFrame(val_pred,columns=['overpotential'])
with pd.ExcelWriter(f'../results/train_result_all.xlsx') as writer:
    y_train.to_excel(writer, sheet_name='y_train', index=False)
    train_pred.to_excel(writer, sheet_name='train_pred', index=False)
    y_val.to_excel(writer, sheet_name='y_val', index=False)
    val_pred.to_excel(writer, sheet_name='val_pred', index=False)


import shap
def combined_model(x_np):
    x_norm = norm_x.transform(x_np)
    x_tensor = torch.Tensor(x_norm)
    with torch.no_grad():
        intermediate = premodel.forward(x_tensor)
        y_tensor = model.forward(intermediate)
    return y_tensor.detach().numpy()

background = x_train.values
explainer = shap.KernelExplainer(combined_model, background)
x_explain = x_val.values
shap_values_list = explainer.shap_values(x_explain)

if isinstance(shap_values_list, list):
    shap_values_arr = shap_values_list[0]
else:
    shap_values_arr = shap_values_list

if shap_values_arr.ndim == 3:
    shap_values_2d = np.squeeze(shap_values_arr, axis=2)
else:
    shap_values_2d = shap_values_arr

plt.figure()
shap.summary_plot(shap_values_2d, x_explain, feature_names=x.columns, show=False)
plt.savefig("../results/shap_summary_plot.png", dpi=300, bbox_inches='tight')
plt.close()

df_input = pd.DataFrame(x_explain, columns=x.columns)
df_shap = pd.DataFrame(shap_values_2d, columns=["shap_" + col for col in x.columns])
df_combined = pd.concat([df_input, df_shap], axis=1)
df_combined.to_excel("../results/shap_data.xlsx", index=False)


