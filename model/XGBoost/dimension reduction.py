import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
import numpy as np

# 设置中文字体支持
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['axes.unicode_minus'] = False

data = pd.read_excel('analysis.xlsx')
features = data.iloc[:, :5]
target = data.iloc[:, -1]


pca = PCA(n_components=2, random_state=3)
pca_reduced = pca.fit_transform(features)

pca_df = pd.DataFrame({
    'PCA_Dimension1': pca_reduced[:, 0],
    'PCA_Dimension2': pca_reduced[:, 1],
    'potential': target
})

pca_df.to_csv('pca_results.csv', index=False)

plt.figure(figsize=(10, 8))
scatter = plt.scatter(pca_reduced[:, 0], pca_reduced[:, 1], c=target, cmap='jet')
colorbar = plt.colorbar(scatter)
colorbar.set_label('potential', fontsize=25)
colorbar.ax.tick_params(labelsize=20)
plt.title('PCA Dimensionality Reduction Visualization', fontdict={'fontsize': 25})
plt.xlabel('Dimension 1', fontdict={'fontsize': 25})
plt.ylabel('Dimension 2', fontdict={'fontsize': 25})
plt.tick_params(labelsize=20)
plt.savefig('pca.png')
plt.show()
