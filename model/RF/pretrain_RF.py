import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestRegressor
import shap
import matplotlib.pyplot as plt
import pickle
import nni


def get_default_parameters():
    # 随机森林模型的默认参数
    params = {
        'seed': 0,
        'n_estimators': 500,  # 树的数量
        'min_samples_split': 2,  # 内部节点再划分所需最小样本数
        'min_samples_leaf': 1,   # 叶子节点最少样本数
        'max_features': 1.0,  # 每次分割时考虑的特征数量
        'min_weight_fraction_leaf':0.0
    }
    return params

config = get_default_parameters()

# 读取训练数据
excel = r'D:\pythonProject\program\urea\add\data\metals_activity.xlsx'
x = pd.read_excel(excel,sheet_name='metals')
y = pd.read_excel(excel,sheet_name='activity')

# 数据划分
x_train, x_test, y_train, y_test = train_test_split(
    x, y, test_size=0.2, random_state=int(config['seed']), shuffle=True
)

# 数据标准化
norm_x = StandardScaler().fit(x_train)
norm_y = StandardScaler().fit(y_train)

x_train_ = norm_x.transform(x_train)
y_test_ = norm_y.transform(y_test)
x_test_ = norm_x.transform(x_test)
y_train_ = norm_y.transform(y_train)

# 创建随机森林模型 - 随机森林原生支持多输出回归
model = RandomForestRegressor(
    n_estimators=int(config['n_estimators']),
    min_samples_split=int(config['min_samples_split']),
    min_samples_leaf=int(config['min_samples_leaf']),
    max_features=config['max_features'],
    min_weight_fraction_leaf=config['min_weight_fraction_leaf'],
    n_jobs=-1  # 使用所有可用的CPU核心
)

# 训练模型
model.fit(x_train_, y_train_)

# 预测
train_pred = model.predict(x_train_)
test_pred = model.predict(x_test_)

train_pred = norm_y.inverse_transform(train_pred)
test_pred = norm_y.inverse_transform(test_pred)

# 反标准化预测结果
train_corr = [np.corrcoef(y_train.iloc[:, i], train_pred[:, i])[0, 1] for i in range(4)]#返回预测值和真实值之间的相关系数
test_corr = [np.corrcoef(y_test.iloc[:, i], test_pred[:, i])[0, 1] for i in range(4)]
train_rmse = [np.sqrt(np.mean((y_train.iloc[:, i] - train_pred[:, i]) ** 2))for i in range(4)]  # 计算训练集的 RMSE
test_rmse = [np.sqrt(np.mean((y_test.iloc[:, i] - test_pred[:, i]) ** 2))for i in range(4)] # 计算验证集的 RMSE
print(train_corr, test_corr, train_rmse, test_rmse)

fig, axes = plt.subplots(4, 2, figsize=(9, 12), tight_layout=True)#fig窗口为3*2，
axes[0, 0].set_title('Train', fontsize=16)
axes[0, 1].set_title('Test', fontsize=16)

for i in range(y.shape[1]):
    y_min = min(y.iloc[:, i].min(), train_pred[:, i].min(), test_pred[:, i].min()) #返回三个数据集中的最小值
    y_max = max(y.iloc[:, i].max(), train_pred[:, i].max(), test_pred[:, i].max())
    lim = [y_min - (y_max - y_min) * 0.05, y_max + (y_max - y_min) * 0.05] #设置坐标轴的范围

    for j in [0, 1]:
        axes[i, j].set_aspect('equal')
        axes[i, j].plot(lim, lim, lw=1)
        axes[i, j].set_xticks(axes[i, j].get_xticks())
        axes[i, j].set_yticks(axes[i, j].get_xticks())
        axes[i, j].set_xlim(lim)
        axes[i, j].set_ylim(lim)

    axes[i, 0].scatter(y_train.iloc[:, i], train_pred[:, i], s=5)#散点图
    axes[i, 1].scatter(y_test.iloc[:, i], test_pred[:, i], s=5)

    axes[i, 0].text(y_max, y_min, f'r={train_corr[i]:.4f}\nRMSE={train_rmse[i]:.4f}', horizontalalignment='right', fontsize=12)
    axes[i, 1].text(y_max, y_min, f'r={test_corr[i]:.4f}\nRMSE={test_rmse[i]:.4f}', horizontalalignment='right', fontsize=12)

    axes[i, 0].set_ylabel(y.columns[i], fontsize=16)

plt.savefig(f'pretrain_RF.png')
plt.show()

# 保存模型
with open('rf_premodel.pkl', 'wb') as f:
    pickle.dump(model, f)


# 保存结果到Excel
train_pred_df = pd.DataFrame(train_pred,columns=['ade_urea','ade_CO','e_urea','e_CO'])
test_pred_df = pd.DataFrame(test_pred,columns=['ade_urea','ade_CO','e_urea','e_CO'])
with pd.ExcelWriter(f'pretrain_RF_result_all.xlsx') as writer:
    y_train.to_excel(writer, sheet_name='y_train', index=False)
    train_pred_df.to_excel(writer, sheet_name='train_pred', index=False)
    y_test.to_excel(writer, sheet_name='y_test', index=False)
    test_pred_df.to_excel(writer, sheet_name='test_pred', index=False)

# 保存标准化器
with open('pretrain_RF_norm_y.pckl', 'wb') as f:
    pickle.dump(norm_y, f)

with open('pretrain_RF_norm_x.pckl', 'wb') as f:
    pickle.dump(norm_x, f)

print("随机森林模型训练和评估完成")