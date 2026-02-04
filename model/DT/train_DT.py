import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.tree import DecisionTreeRegressor
import matplotlib.pyplot as plt
import pickle
import nni
class Predictor:
    def __init__(self, model_path):
        """初始化预测器"""
        with open(model_path, 'rb') as f:
            self.model = pickle.load(f)

    def predict(self, x_data):
        y_scaled = self.model.predict(x_data)
        return y_scaled


def get_default_parameters():
    # 决策树模型的默认参数
    params = {
        'seed1': 106, #25,57,93,16,42,131,197,107
        'max_depth':6,
        'min_samples_split':4,
        'min_samples_leaf':1,
        'max_leaf_nodes':64,
        'max_features':1,
        'seed2':77
    }
    return params

# received_PARAMS = nni.get_next_parameter()
# config = get_default_parameters()
# config.update(received_PARAMS)

config = get_default_parameters()

# 读取训练数据
excel = 'D:\\pythonProject\\program\\urea\\add\\data\\experiment.xlsx'
x = pd.read_excel(excel, sheet_name='metals')
y = pd.read_excel(excel, sheet_name='potential')

# 数据划分
x_train, x_test, y_train, y_test = train_test_split(
    x, y, test_size=0.2, random_state=int(config['seed1']), shuffle=True
)

# 数据标准化
norm_x = StandardScaler().fit(x_train)
norm_y = StandardScaler().fit(y_train)

# 标准化数据
x_train_ = norm_x.transform(x_train)
x_test_ = norm_x.transform(x_test)
y_train_ = norm_y.transform(y_train).ravel()  # 转换为1D数组
y_test_ = norm_y.transform(y_test).ravel()    # 转换为1D数组

# 使用预训练模型获取特征
predictor = Predictor('dt_premodel.pkl')
train_pred_ = predictor.predict(x_train_)
test_pred_ = predictor.predict(x_test_)

# 确保特征维度正确
if len(train_pred_.shape) == 1:
    train_pred_ = train_pred_.reshape(-1, 1)
if len(test_pred_.shape) == 1:
    test_pred_ = test_pred_.reshape(-1, 1)

# 创建决策树模型
model = DecisionTreeRegressor(
    max_depth=config['max_depth'],
    min_samples_leaf=config['min_samples_leaf'],
    min_samples_split=config['min_samples_split'],
    max_leaf_nodes=config['max_leaf_nodes'],
    max_features=config['max_features'],
    random_state=int(config['seed2'])
)

# 训练模型 - 使用1D的y
model.fit(train_pred_, y_train_)

# 预测
train_pred_scaled = model.predict(train_pred_)
test_pred_scaled = model.predict(test_pred_)

# 将预测结果转换为2D数组以便反标准化
train_pred_scaled_2d = train_pred_scaled.reshape(-1, 1)
test_pred_scaled_2d = test_pred_scaled.reshape(-1, 1)

# 反标准化
train_pred = norm_y.inverse_transform(train_pred_scaled_2d)
test_pred = norm_y.inverse_transform(test_pred_scaled_2d)

# 计算指标
train_corr = np.corrcoef(y_train.iloc[:, 0], train_pred.ravel())[0, 1]  # 使用ravel()转换为1D
test_corr = np.corrcoef(y_test.iloc[:, 0], test_pred.ravel())[0, 1]
train_rmse = np.sqrt(np.mean((y_train.iloc[:, 0] - train_pred.ravel()) ** 2))
test_rmse = np.sqrt(np.mean((y_test.iloc[:, 0] - test_pred.ravel()) ** 2))
print(f"Train - Corr: {train_corr:.4f}, RMSE: {train_rmse:.4f}")
print(f"Test - Corr: {test_corr:.4f}, RMSE: {test_rmse:.4f}")

if train_corr > test_corr:
    nni.report_final_result(test_corr)
else:
    nni.report_final_result(0)

# 绘图
fig, axes = plt.subplots(1, 2, figsize=(12, 6), tight_layout=True)  # 调整图形大小
axes[0].set_title('Train', fontsize=20)
axes[1].set_title('Test', fontsize=20)

# 计算坐标轴范围
y_min = min(y.iloc[:, 0].min(), train_pred.min(), test_pred.min())
y_max = max(y.iloc[:, 0].max(), train_pred.max(), test_pred.max())
lim = [y_min - (y_max - y_min) * 0.05, y_max + (y_max - y_min) * 0.05]

for j in range(2):
    axes[j].set_aspect('equal')
    axes[j].plot(lim, lim, lw=2, color='red', alpha=0.7)
    axes[j].set_xlim(lim)
    axes[j].set_ylim(lim)
    axes[j].tick_params(labelsize=12)

# 绘制散点图
axes[0].scatter(y_train.iloc[:, 0], train_pred.ravel(), s=30, alpha=0.6)
axes[1].scatter(y_test.iloc[:, 0], test_pred.ravel(), s=30, alpha=0.6)

# 添加文本标签
axes[0].text(lim[0] + (lim[1]-lim[0])*0.05, lim[1] - (lim[1]-lim[0])*0.05,
             f'r={train_corr:.4f}\nRMSE={train_rmse:.4f}',
             horizontalalignment='left', verticalalignment='top',
             fontsize=14, bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))

axes[1].text(lim[0] + (lim[1]-lim[0])*0.05, lim[1] - (lim[1]-lim[0])*0.05,
             f'r={test_corr:.4f}\nRMSE={test_rmse:.4f}',
             horizontalalignment='left', verticalalignment='top',
             fontsize=14, bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))

# 修正标签设置
axes[0].set_ylabel('Predicted ' + y.columns[0], fontsize=14)
axes[0].set_xlabel('Actual', fontsize=14)
axes[1].set_xlabel('Actual', fontsize=14)

plt.savefig('train_DT.png', dpi=300, bbox_inches='tight')
plt.show()

# 保存结果到Excel
train_pred_df = pd.DataFrame(train_pred, columns=['potential'])
test_pred_df = pd.DataFrame(test_pred, columns=['potential'])

with pd.ExcelWriter('train_DT_result_all.xlsx') as writer:
    y_train.to_excel(writer, sheet_name='y_train', index=False)
    train_pred_df.to_excel(writer, sheet_name='train_pred', index=False)
    y_test.to_excel(writer, sheet_name='y_test', index=False)
    test_pred_df.to_excel(writer, sheet_name='test_pred', index=False)

# 保存模型
with open('dt_model.pkl', 'wb') as f:
    pickle.dump(model, f)

# 保存标准化器
with open('train_DT_norm_y.pckl', 'wb') as f:
    pickle.dump(norm_y, f)

with open('train_DT_norm_x.pckl', 'wb') as f:
    pickle.dump(norm_x, f)

print("决策树模型训练和评估完成")

txt = open(f'D:\\pythonProject\\program\\urea\\add\\data\\choose_proportion.txt')
lines = txt.readlines()
for i in range(len(lines)):
    lines[i]=lines[i].split()
    if i !=0:
        for j in range(len(lines[i])):
            lines[i][j]=float(lines[i][j])/100

x_val=np.array(lines)
x_val=np.delete(x_val,0,0)
x_val_ = norm_x.transform(x_val)
val_pred_medium = predictor.predict(x_val_)
val_pred_final = model.predict(val_pred_medium)

# 修复：将1D数组转换为2D数组
val_pred_final = val_pred_final.reshape(-1, 1)
val_pred_final = norm_y.inverse_transform(val_pred_final)

x_val=pd.DataFrame(x_val,columns=['Mn','Fe','Co','Ni','Zn'])
val_pred_final=pd.DataFrame(val_pred_final,columns=['potential'])
with pd.ExcelWriter(f'space_one_percent.xlsx') as writer:
    x_val.to_excel(writer, sheet_name='ratio', index=False)
    val_pred_final.to_excel(writer, sheet_name='value', index=False)