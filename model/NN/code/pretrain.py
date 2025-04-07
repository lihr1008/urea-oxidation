import torch
import torch.nn as nn
from torch.nn import functional as F
from torch.nn import Module
import torch.optim as optim
from pytorchtools import EarlyStopping
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from torch.utils.data import Dataset,DataLoader
from matplotlib import pyplot as plt

excel = r'../../data/metals_activity.xlsx'
x = pd.read_excel(excel,sheet_name='metals')
y = pd.read_excel(excel,sheet_name='activity')
x_train,x_val,y_train,y_val = train_test_split(x,y,test_size=0.2,random_state=0)

norm_x = StandardScaler().fit(x_train)
norm_y = StandardScaler().fit(y_train)
x_train_ = norm_x.transform(x_train)
x_val_ = norm_x.transform(x_val)

y_train_ = norm_y.transform(y_train)
y_val_ = norm_y.transform(y_val)

class GetLoader(Dataset):
    # 初始化函数，得到数据
    def __init__(self, data_root, data_label):
        self.data = data_root
        self.label = data_label

    def __getitem__(self, index):
        datas = self.data[index]
        labels = self.label[index]
        return (datas, labels)

    def __len__(self):
        return len(self.data)


trainDataset = GetLoader(x_train_,y_train_)
testDataset = GetLoader(x_val_,y_val_)

train_loader = DataLoader(
    trainDataset,
    batch_size=256,
    shuffle=False
)

test_loader = DataLoader(
    testDataset,
    batch_size=256,
    shuffle=False
)

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

optimizer = optim.Adam(premodel.parameters(),lr=1e-3)
scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(optimizer,'min',factor=0.5, min_lr=1e-6,patience=10)
early_stopping = EarlyStopping(delta=1e-5, patience=15)
criteon = nn.MSELoss()
epochs = 500
for epoch in range(epochs):
    temp_loss = []
    for batch_idx,train_data in enumerate(train_loader):
        # premodel.train()
        (data, label) = train_data
        data=data.to(torch.float32)
        label = label.to(torch.float32)
        logits = premodel(data)
        #logits = model.train(data)
        loss = criteon(logits,label)
        temp_loss.append(loss.data)
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()
        #scheduler.step(loss)
    temp_loss = np.array(temp_loss)
    train_loss = np.average(temp_loss)
    print('Train Epoch : {}\tLoss:{}'.format(
                epoch,train_loss
        ))

    temp_loss = []
    test_loss = 0.0
    correct = 0
    for index,test_data in enumerate(test_loader):
        premodel.eval()
        (data,label) = test_data
        data = data.to(torch.float32)
        label = label.to(torch.float32)
        logits = premodel(data)
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

    early_stopping(val_loss,premodel)
    if early_stopping.early_stop:
        print("Early stopping")
        break


print("完成")


x_train_ = torch.tensor(x_train_).to(torch.float32)
x_val_= torch.tensor(x_val_).to(torch.float32)

train_pred = premodel.forward(x_train_)
val_pred = premodel.forward(x_val_)


train_pred = train_pred.detach().numpy()
train_pred = norm_y.inverse_transform(train_pred)
val_pred = val_pred.detach().numpy()
val_pred = norm_y.inverse_transform(val_pred)

train_corr = [np.corrcoef(y_train.iloc[:, i], train_pred[:, i])[0, 1] for i in range(4)]#返回预测值和真实值之间的相关系数
val_corr = [np.corrcoef(y_val.iloc[:, i], val_pred[:, i])[0, 1] for i in range(4)]
print(train_corr,val_corr)

fig, axes = plt.subplots(4, 2, figsize=(9, 12), tight_layout=True)#fig窗口为3*2，
axes[0, 0].set_title('Train', fontsize=16)
axes[0, 1].set_title('Validation', fontsize=16)

for i in range(y.shape[1]):
    y_min = min(y.iloc[:, i].min(), train_pred[:, i].min(), val_pred[:, i].min()) #返回三个数据集中的最小值
    y_max = max(y.iloc[:, i].max(), train_pred[:, i].max(), val_pred[:, i].max())
    lim = [y_min - (y_max - y_min) * 0.05, y_max + (y_max - y_min) * 0.05] #设置坐标轴的范围

    for j in [0, 1]:
        axes[i, j].set_aspect('equal')
        axes[i, j].plot(lim, lim, lw=1)
        axes[i, j].set_xticks(axes[i, j].get_xticks())
        axes[i, j].set_yticks(axes[i, j].get_xticks())
        axes[i, j].set_xlim(lim)
        axes[i, j].set_ylim(lim)

    axes[i, 0].scatter(y_train.iloc[:, i], train_pred[:, i], s=5)#散点图
    axes[i, 1].scatter(y_val.iloc[:, i], val_pred[:, i], s=5)

    axes[i, 0].text(y_max, y_min, f'r={train_corr[i]:.3f}', horizontalalignment='right', fontsize=12)
    axes[i, 1].text(y_max, y_min, f'r={val_corr[i]:.3f}', horizontalalignment='right', fontsize=12)

    axes[i, 0].set_ylabel(y.columns[i], fontsize=16)

plt.savefig(f'../results/pretrain.png')
plt.show()

torch.save(premodel, '../results/PBA-ml-pretrain_model.pkl')


