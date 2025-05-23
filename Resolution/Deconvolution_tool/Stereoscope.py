!pip install stereoscope
import stereoscope as st
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# --------------------------
# 1. load data
# --------------------------

# generate simulated single cell data (1000 cells, 200 genes)
sc_data = sc.AnnData(
    X=np.random.negative_binomial(n=5, p=0.1, size=(1000, 200)),  # 模拟UMI计数
    obs={"cell_type": np.random.choice(["Neuron", "Astrocyte", "Microglia"], size=1000)}  # 细胞类型标签
)

# 生成模拟空间数据 (500 spots, 200 genes)
spatial_data = sc.AnnData(
    X=np.random.negative_binomial(n=5, p=0.1, size=(500, 200)),  # 模拟UMI计数
    obsm={"spatial": np.random.rand(500, 2)}  # 空间坐标
)

# --------------------------
# 2. 数据预处理（关键步骤）
# --------------------------
# 确保单细胞和空间数据的基因一致
shared_genes = list(set(sc_data.var_names) & set(spatial_data.var_names))
sc_data = sc_data[:, shared_genes].copy()
spatial_data = spatial_data[:, shared_genes].copy()

# 标准化（Stereoscope要求原始计数，无需对数变换）
sc.pp.normalize_total(sc_data, target_sum=1e4)
sc.pp.normalize_total(spatial_data, target_sum=1e4)

# 将单细胞数据转换为细胞类型×基因的伪bulk矩阵
cell_type_matrix = pd.crosstab(
    index=sc_data.obs["cell_type"],
    columns=sc_data.var_names,
    values=sc_data.X.toarray().sum(axis=0),  # 假设X是稀疏矩阵
    aggfunc="sum"
)
sc_data_pseudo = sc.AnnData(cell_type_matrix)  # 伪bulk数据

# --------------------------
# 3. 训练 Stereoscope 模型
# --------------------------
model = st.Stereoscope(
    sc_data=sc_data_pseudo,      # 单细胞伪bulk数据（细胞类型×基因）
    st_data=spatial_data,        # 空间数据
    n_components=10,            # 潜在空间维度（默认10）
    use_gpu=True                # 是否使用GPU加速
)
model.fit(
    n_epochs=100,               # 训练轮次
    lr=0.01,                    # 学习率
    batch_size=64               # 批大小
)

# --------------------------
# 4. 预测 spot 细胞组成
# --------------------------
results = model.predict()
print(results.head())  # 输出每个spot的细胞类型比例（DataFrame）

# --------------------------
# 5. 可视化结果（可选）
# --------------------------
# 绘制空间分布（以Neuron比例为例）
plt.figure(figsize=(10, 6))
plt.scatter(
    spatial_data.obsm["spatial"][:, 0],
    spatial_data.obsm["spatial"][:, 1],
    c=results["Neuron"],        # 颜色映射到Neuron比例
    s=20, cmap="viridis", alpha=0.8
)
plt.colorbar(label="Neuron Proportion")
plt.title("Stereoscope Deconvolution Results")
plt.show()
