# -*- coding: UTF-8 -*-

import matplotlib.pyplot as plt
import pandas as pd
import shap

excel = '../results/new_shap_data.xlsx'
shap_data = pd.read_excel(excel,sheet_name='shap')
shap_values = shap_data.values
x_val=pd.read_excel(excel,sheet_name='x_val')
x_explain = x_val.values      # 选择需要解释的样本（例如验证集全部）

# 绘制summary plot，展示各特征对预测结果的影响
plt.figure()
shap.summary_plot(shap_values, x_explain, feature_names=x_val.columns,show=False)
# 获取当前figure并保存为png格式，dpi可根据需要调整
plt.savefig("../results/new_shap_summary_plot.png", dpi=300, bbox_inches='tight')

