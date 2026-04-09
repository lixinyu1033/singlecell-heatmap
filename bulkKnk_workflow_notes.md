# bulkKnk 代码工作流与原理梳理

仓库：<https://github.com/fengsh27mail/bulkKnk>

## 总体定位
- `bulkKnk` 是一个 R 包，目标是基于 bulk RNA-seq 做 **in silico 基因敲除（virtual knockout）** 与 **表型逆转评估**。
- 方法链路按 README 的 Quick Start 可概括为：
  1) 特征筛选（HVG + 生物签名保送）
  2) WGCNA 模块/Hub 基因识别
  3) GENIE3 因果网络推断 + ARACNE/DPI 剪枝
  4) 基于网络传播/RWR 的虚拟敲除
  5) 机器学习与生存分析评估干预收益

## 关键依赖（对应算法栈）
- `DESCRIPTION` 中可见主要依赖：`WGCNA`、`GENIE3`、`randomForest`、`clusterProfiler`、`org.Hs.eg.db`、`ggplot2`、`patchwork`、`reshape2`。
- 这对应了“网络构建 + 干预传播 + 预测评估 + 富集可视化”的完整流水线。

## 代码模块与职责

### 0) 输入与数据对接
- `import_xena_cohort()`：用于导入 UCSC Xena 风格表达矩阵与生存文件。
- README 同时给了内置 demo（TCGA-STAD）的快速跑通路径。

### 1) 特征筛选（Module 0.1）
- `select_advanced_features(expr_matrix, n_features, fallback_threshold)`：
  - 先尝试保留内置签名基因（来自 `predict_ko_phenotype` 内部签名定义）；
  - 若匹配太少（低于阈值），自动降级为纯 HVG 模式（按标准差排序补齐）；
  - 用意是兼顾“人类签名先验”与“跨物种通用性”。

### 2) 枢纽模块与靶点候选
- `identify_hub_modules(sub_expr, trait_vec)`：
  - 用 WGCNA 在高维表达里先找与目标性状（如 OS_status）相关模块；
  - 再从模块内给出 hub 候选基因（README 示例直接取第一名作为 `top_target`）。

### 3) 因果网络推断与净化
- `infer_causal_network(focus_mat)`：
  - 用 GENIE3 从表达矩阵学习有向调控权重。
- `prune_network_dpi_fast()` / `prune_network_dpi()` / `refine_network_aracne()`：
  - 用 ARACNE 的 DPI 思路做边剪枝，减少间接边，得到更稀疏且更“直接”的网络。
- README 中随后会做列归一化，形成传播用转移矩阵 `W_norm`。

### 4) 虚拟敲除核心引擎
- `run_virtual_knockout(target_genes, E_init, adj_matrix, alpha, steps)`：
  - 若指定靶点：
    - 把靶点基线表达置 0；
    - 同时把靶点相关入边/出边置 0（相当于“物理断边”）；
  - 然后迭代：
    - `E_new = alpha * (W %*% E_current) + (1 - alpha) * E_anchor`
  - 本质上是带重启的随机游走/网络传播求稳态。
- `batch_virtual_knockout(expr_matrix, pheno_vec, target_pheno, target_genes, adj_matrix, method)`：
  - 先从队列中选中目标表型子集（例如 Tumor）；
  - 用该子集均值做 `E_baseline`；
  - 跑 `run_virtual_knockout` 得到群体 KO 稳态向量；
  - 再用 `build_ko_matrix` 把群体层面的偏移“桥接”回样本层（`shift` 或 `replace`）。

### 5) 表型逆转与临床收益评估
- `predict_ko_phenotype(WT_matrix, KO_matrix, ko_gene, clinical_data, ...)`：
  - 使用 AI 分类器（仓库介绍写明为 Random Forest）评估 KO 后样本向“健康/敏感/低风险”方向逆转的概率；
  - 结合生存信息输出风险变化、逆转率、获益概率等指标。
- 下游可视化函数：
  - `plot_km_comparison`
  - `plot_cox_forest`
  - `plot_risk_score_shift`
  - `plot_ko_delta`
  - `plot_dual_volcano`
- 还有 `run_vk_enrichment()` 用于 KO 前后差异相关通路富集解释。

## 一条“从数据到结论”的典型执行路径
1. 读取表达和临床数据（demo 或 Xena）。
2. `select_advanced_features` 压缩维度并保留关键生物签名。
3. `identify_hub_modules` 找到与表型最相关模块并提名 hub 靶点。
4. 取核心基因子集，`infer_causal_network` 推断网络，再 `prune_network_dpi_fast` 去间接边。
5. 构建归一化传播矩阵后，`batch_virtual_knockout` 生成 KO 队列表达。
6. `predict_ko_phenotype` + 生存分析/可视化判断“是否逆转、逆转多少、临床是否可能获益”。

## 方法原理（简化理解）
- **为什么先 WGCNA 再 GENIE3？**
  - WGCNA 先做结构化降维与候选缩圈；GENIE3 在更聚焦的基因集上学有向依赖，降低噪声与计算成本。
- **为什么要 DPI 剪枝？**
  - 机器学习网络常含大量间接边；DPI 近似剔除 A→C 这类可被 A→B→C 解释的弱直接关系，减少传播时的伪级联。
- **虚拟 KO 如何作用到全网？**
  - 不是只把一个基因设零，而是“断开网络连接 + 迭代传播到稳态”，因此可模拟全转录组级联反应。
- **为何还要样本级桥接？**
  - 核心传播通常在“群体均值状态”上完成；`build_ko_matrix` 把群体偏移映射回个体样本，方便后续分类、生存、可视化。

## 使用时的注意点（实操）
- 表达矩阵方向要统一（README 中不同函数有行列方向约定，调用前要核对）。
- 网络矩阵和表达向量的基因名必须严格对齐（`run_virtual_knockout` 会校验）。
- `target_pheno` 需要在 `pheno_vec` 中真实存在，否则会报错。
- 非人类或面板测序建议显式提供 `custom_signatures`，避免完全依赖默认内置签名。

## 我对仓库实现风格的判断
- 这是一个“方法整合型”工具包：把 WGCNA/GENIE3/ARACNE/RWR/RF 串成可执行工程流。
- README 的示例路径非常接近作者心中的“标准工作流”。
- 文档层面，`man/*.Rd` 可快速确认函数接口；源码里部分文件呈现压缩风格（阅读性一般），但不影响理解主算法脉络。
