# -*- coding: UTF-8 -*-

import matplotlib.pyplot as plt
import pandas as pd
import shap

excel = '../results/new_shap_data.xlsx'
shap_data = pd.read_excel(excel,sheet_name='shap')
shap_values = shap_data.values
x_val=pd.read_excel(excel,sheet_name='x_val')
x_explain = x_val.values 

plt.figure()
shap.summary_plot(shap_values, x_explain, feature_names=x_val.columns,show=False)
plt.savefig("../results/new_shap_summary_plot.png", dpi=300, bbox_inches='tight')

