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
shap.summary_plot(
    shap_values,
    x_explain,
    feature_names=x_val.columns,
    show=False 
)

fig = plt.gcf()  
axes = fig.axes  
ax_main = axes[0]
ax_main.set_xlabel(
    ax_main.get_xlabel(),
    fontweight='bold',
    fontsize=12 
)

for tick in ax_main.get_yticklabels():
    tick.set_fontweight('bold')
    tick.set_fontsize(12)  

for tick in ax_main.get_xticklabels():
    tick.set_fontweight('bold')
    tick.set_fontsize(12) 

if len(axes) > 1: 
    ax_colorbar = axes[1]
    ax_colorbar.set_ylabel(
        ax_colorbar.get_ylabel(),
        fontweight='bold',
        fontsize=12
    )

    for tick in ax_colorbar.get_yticklabels():
        tick.set_fontweight('bold')
        tick.set_fontsize(12)


plt.savefig(
    "../results/new_shap_summary_plot.png",
    dpi=300,
    bbox_inches='tight' 
)
plt.close() 
