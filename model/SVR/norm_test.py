import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVR
import matplotlib.pyplot as plt
import pickle
import nni

class SVRPredictor:
    def __init__(self, model_path):
        """初始化预测器"""
        with open(model_path, 'rb') as f:
            self.model = pickle.load(f)

    def predict(self, x_data):

        y_scaled = self.model.predict(x_data)

        return y_scaled


# 读取训练数据
excel = 'D:\\pythonProject\\program\\urea\\add\\data\\experiment.xlsx'  # 使用双反斜杠或原始字符串
x = pd.read_excel(excel, sheet_name='metals')
y = pd.read_excel(excel, sheet_name='potential')

# 数据划分
x_train, x_test, y_train, y_test = train_test_split(
    x, y, test_size=0.2, random_state=106, shuffle=True
)

# 数据标准化
norm_x = StandardScaler().fit(x_train)

# 读取训练数据
excel = r'D:\pythonProject\program\urea\add\data\metals_activity.xlsx'
x = pd.read_excel(excel,sheet_name='metals')
y = pd.read_excel(excel,sheet_name='activity')

# 数据划分
x_train, x_test, y_train, y_test = train_test_split(
    x, y, test_size=0.2, random_state=0, shuffle=True
)

norm_y = StandardScaler().fit(y_train)

# 标准化数据
x_train_ = norm_x.transform(x_train)
x_test_ = norm_x.transform(x_test)
y_train_ = norm_y.transform(y_train).ravel()  # 转换为1D数组
y_test_ = norm_y.transform(y_test).ravel()    # 转换为1D数组

# 使用预训练模型获取特征
predictor = SVRPredictor('svr_premodel.pkl')
train_pred = predictor.predict(x_train_)
test_pred = predictor.predict(x_test_)

train_pred = norm_y.inverse_transform(train_pred)
test_pred = norm_y.inverse_transform(test_pred)

# 反标准化预测结果
train_corr = [np.corrcoef(y_train.iloc[:, i], train_pred[:, i])[0, 1] for i in range(4)]#返回预测值和真实值之间的相关系数
test_corr = [np.corrcoef(y_test.iloc[:, i], test_pred[:, i])[0, 1] for i in range(4)]
train_rmse = [np.sqrt(np.mean((y_train.iloc[:, i] - train_pred[:, i]) ** 2))for i in range(4)]  # 计算训练集的 RMSE
test_rmse = [np.sqrt(np.mean((y_test.iloc[:, i] - test_pred[:, i]) ** 2))for i in range(4)] # 计算验证集的 RMSE
print(train_corr, test_corr, train_rmse, test_rmse)

