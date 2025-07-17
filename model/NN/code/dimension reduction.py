import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
import umap
import numpy as np

# 设置中文字体支持
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['axes.unicode_minus'] = False

# 读取数据
data = pd.read_excel('space_one_percent_.xlsx')

# 提取特征和目标变量
features = data.iloc[:, :5]  # 前5列为特征
target = data.iloc[:, -1]    # 最后一列为目标变量

# ----------------------------
# PCA 降维并保存结果
pca = PCA(n_components=2, random_state=3)
pca_reduced = pca.fit_transform(features)

# 创建包含PCA结果和目标变量的DataFrame
pca_df = pd.DataFrame({
    'PCA_Dimension1': pca_reduced[:, 0],
    'PCA_Dimension2': pca_reduced[:, 1],
    'Overpotential': target
})

# 保存为CSV文件
pca_df.to_csv('pca_results.csv', index=False)
print("PCA结果已保存到 pca_results.csv")

# PCA可视化
plt.figure(figsize=(10, 8))
scatter = plt.scatter(pca_reduced[:, 0], pca_reduced[:, 1], c=target, cmap='jet')
colorbar = plt.colorbar(scatter)
colorbar.set_label('overpotential', fontsize=25)
colorbar.ax.tick_params(labelsize=20)
plt.title('PCA Dimensionality Reduction Visualization', fontdict={'fontsize': 25})
plt.xlabel('Dimension 1', fontdict={'fontsize': 25})
plt.ylabel('Dimension 2', fontdict={'fontsize': 25})
plt.tick_params(labelsize=20)
plt.savefig('pca.png')
plt.show()

# # ----------------------------
# # t-SNE 降维并保存结果
# tsne = TSNE(n_components=2, random_state=15)
# tsne_reduced = tsne.fit_transform(features)
#
# # 创建包含t-SNE结果和目标变量的DataFrame
# tsne_df = pd.DataFrame({
#     'tSNE_Dimension1': tsne_reduced[:, 0],
#     'tSNE_Dimension2': tsne_reduced[:, 1],
#     'Overpotential': target
# })
#
# # 保存为CSV文件
# tsne_df.to_csv('tsne_results.csv', index=False)
# print("t-SNE结果已保存到 tsne_results.csv")
#
# # t-SNE可视化
# plt.figure(figsize=(10, 8))
# scatter = plt.scatter(tsne_reduced[:, 0], tsne_reduced[:, 1], c=target, cmap='jet')
# colorbar = plt.colorbar(scatter)
# colorbar.set_label('overpotential', fontsize=25)
# colorbar.ax.tick_params(labelsize=20)
# plt.title('t-SNE Dimensionality Reduction Visualization', fontdict={'fontsize': 25})
# plt.xlabel('Dimension 1', fontdict={'fontsize': 25})
# plt.ylabel('Dimension 2', fontdict={'fontsize': 25})
# plt.tick_params(labelsize=20)
# plt.savefig('tsne.png')
# plt.show()
#
# # ----------------------------
# # UMAP 降维并保存结果
# umap_model = umap.UMAP(n_components=2, random_state=15)
# umap_reduced = umap_model.fit_transform(features)
#
# # 创建包含UMAP结果和目标变量的DataFrame
# umap_df = pd.DataFrame({
#     'UMAP_Dimension1': umap_reduced[:, 0],
#     'UMAP_Dimension2': umap_reduced[:, 1],
#     'Overpotential': target
# })
#
# # 保存为CSV文件
# umap_df.to_csv('umap_results.csv', index=False)
# print("UMAP结果已保存到 umap_results.csv")
#
# # UMAP可视化
# plt.figure(figsize=(10, 8))
# scatter = plt.scatter(umap_reduced[:, 0], umap_reduced[:, 1], c=target, cmap='jet')
# colorbar = plt.colorbar(scatter)
# colorbar.set_label('overpotential', fontsize=25)
# colorbar.ax.tick_params(labelsize=20)
# plt.title('UMAP Dimensionality Reduction Visualization', fontdict={'fontsize': 25})
# plt.xlabel('Dimension 1', fontdict={'fontsize': 25})
# plt.ylabel('Dimension 2', fontdict={'fontsize': 25})
# plt.tick_params(labelsize=20)
# plt.savefig('umap.png')
# plt.show()

# 保存特征重要性(仅PCA有明确的特征重要性)
if hasattr(pca, 'components_'):
    pca_importance = pd.DataFrame(
        pca.components_,
        columns=features.columns,
        index=[f'PC{i+1}' for i in range(pca.n_components)]
    )
    pca_importance.to_csv('pca_feature_importance.csv')
    print("PCA特征重要性已保存到 pca_feature_importance.csv")