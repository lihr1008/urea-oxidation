# -*- coding: UTF-8 -*-
import shap
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.ticker as ticker
from matplotlib.cm import ScalarMappable
import matplotlib
matplotlib.use('TkAgg')
import pandas as pd
import numpy as np

plt.rcParams['font.serif'] = ['Arial']
plt.rcParams['font.weight'] = 'bold'
plt.rcParams['axes.labelweight'] = 'bold'
plt.rcParams['axes.titleweight'] = 'bold'
plt.rcParams['xtick.labelsize'] = 12
plt.rcParams['ytick.labelsize'] = 12
plt.rcParams['axes.linewidth'] = 2
plt.rcParams['xtick.major.width'] = 2
plt.rcParams['ytick.major.width'] = 2
plt.rcParams['xtick.minor.width'] = 2.0
plt.rcParams['ytick.minor.width'] = 1.0

excel = '../results/new_shap_data.xlsx'
shap_data = pd.read_excel(excel, sheet_name='shap')
shap_values = shap_data.values
x_val = pd.read_excel(excel, sheet_name='x_val')
x_explain = x_val.values

plt.figure()
shap.summary_plot(shap_values, x_explain, feature_names=x_val.columns, show=False)
plt.savefig("../results/new_shap_summary_plot.png", dpi=300, bbox_inches='tight')

data = pd.read_excel(excel, sheet_name='average')
feature_names = data.columns.tolist()
shap_series = data.values.flatten()
sorted_shap_values = pd.read_excel(excel, sheet_name='percent').values.flatten()
num_vars = len(sorted_shap_values)

pink = (253/255, 0/255, 89/255)
blue = (0/255, 134/255, 250/255)
cmap = mcolors.LinearSegmentedColormap.from_list('blue_pink_manual', [blue, pink], N=256)
color_norm = plt.Normalize(data.values.min(), data.values.max())
colors = cmap(color_norm(shap_series))

fig = plt.figure(figsize=(16, 8))
left_margin, right_margin, bottom_margin, top_margin = 0.08, 0.08, 0.12, 0.12
space_between_plots, colorbar_width = 0.04, 0.02
plot_bottom, plot_height = bottom_margin, 1.0 - bottom_margin - top_margin
cbar_left = left_margin
main_ax_left = cbar_left + colorbar_width + space_between_plots
main_ax_width = 1.0 - main_ax_left - right_margin

ax_cbar = fig.add_axes([cbar_left, plot_bottom, colorbar_width, plot_height])
ax_bar = fig.add_axes([main_ax_left, plot_bottom, main_ax_width, plot_height])

sm = ScalarMappable(cmap=cmap, norm=color_norm)
cbar = fig.colorbar(sm, cax=ax_cbar, orientation='vertical')
cbar.set_ticks([])
cbar.ax.yaxis.set_ticks_position('left')
ax_cbar.text(0.5, 1.01, 'High', transform=ax_cbar.transAxes, ha='center', va='bottom', fontsize=24, fontweight="bold")
ax_cbar.text(0.5, -0.01, 'Low', transform=ax_cbar.transAxes, ha='center', va='top', fontsize=24, fontweight="bold")
cbar.outline.set_visible(False)
ax_cbar.text(-1.2, 0.5, 'Mean |SHAP Value|', transform=ax_cbar.transAxes, fontsize=24, rotation=90, va='center', fontweight="bold")

ax_bar.xaxis.tick_bottom()
ax_bar.xaxis.set_label_position("bottom")
ax_bar.invert_xaxis()
ax_bar.barh(range(num_vars), shap_series, color=colors, height=0.6)
ax_bar.invert_yaxis()
ax_bar.set_xlabel('Mean |SHAP Value|', size=24, labelpad=20, fontweight="bold")
ax_bar.set_yticks([])
ax_bar.spines[['left', 'top']].set_visible(False)
ax_bar.spines['right'].set_position(('data', 0))
ax_bar.spines['right'].set_visible(True)
ax_bar.spines['bottom'].set_visible(True)
ax_bar.tick_params(axis='x', which='major', direction='out', labelsize=20, length=6, pad=8)
ax_bar.xaxis.set_minor_locator(ticker.AutoMinorLocator(2))
ax_bar.tick_params(axis='x', which='minor', direction='out', length=4)
for label in ax_bar.get_xticklabels():
    label.set_fontweight('bold')
for i, feature in enumerate(feature_names[:num_vars]):
    ax_bar.text(-0.002, i, feature, ha='left', va='center', color='black', fontsize=24, fontweight="bold")

inset_left = main_ax_left - 0.2
inset_bottom = plot_bottom - 0.08
inset_size = min(main_ax_width, plot_height) * 0.75
ax_radial_inset = fig.add_axes([inset_left, inset_bottom, inset_size, inset_size], projection='polar')
ax_radial_inset.patch.set_alpha(0)

percentages = (sorted_shap_values / sorted_shap_values.sum()) * 100
widths = (sorted_shap_values / sorted_shap_values.sum()) * 2 * np.pi
base_length, fixed_increment, colored_ring_width = 3.0, 0.5, 2.0
total_lengths = [base_length + i * fixed_increment for i in range(num_vars)]
inner_heights = [max(0, tl - colored_ring_width) for tl in total_lengths]
inner_colors = (['#EAEAEA', '#FFFFFF'] * (num_vars // 2 + 1))[:num_vars]
one_oclock_offset = np.pi / 21
thetas = np.cumsum([0] + widths[:-1].tolist()) - one_oclock_offset

ax_radial_inset.bar(thetas, inner_heights, width=widths, color=inner_colors, align='edge', edgecolor='white', linewidth=1.5)
ax_radial_inset.bar(thetas, [colored_ring_width]*num_vars, width=widths, bottom=inner_heights, color=colors, align='edge', edgecolor='white', linewidth=1.5)
for i in range(num_vars):
    label_angle_rad = thetas[i] + widths[i] / 1.8
    label_radius = total_lengths[i] + 1.2
    ax_radial_inset.text(label_angle_rad, label_radius, f'{percentages[i]:.1f}%', ha='center', va='center', fontsize=16, fontweight="bold")

ax_radial_inset.set_yticklabels([])
ax_radial_inset.set_xticklabels([])
ax_radial_inset.spines['polar'].set_visible(False)
ax_radial_inset.grid(False)
ax_radial_inset.set_theta_zero_location('N')
ax_radial_inset.set_theta_direction(-1)
ax_radial_inset.set_ylim(0, max(total_lengths) + 3)

plt.savefig(r'../results/shap_all.png', dpi=300, bbox_inches='tight')
