# Spatial Expression Web Viewer

`web_viewer` 提供空间转录组表达量的本地网页查看器。当前版本使用真实细胞轮廓绘制，并按生物重复拆分展示样本。空间图下方显示当前基因在合并 `seurat_clusters` 上的单行 cluster dotplot。

## 当前显示方式

- S1 拆分为 7 个重复：
  - 第一行：`s1_rep1`、`s1_rep2`、`s1_rep3`、`s1_rep4`
  - 第二行：`s1_rep5`、`s1_rep6`、`s1_rep7`
- S2 拆分为 3 个重复：
  - 一行展示 `s2_rep1`、`s2_rep2`、`s2_rep3`
- 同一样本内所有重复统一显示尺寸：
  - S1 以 `s1_rep1` 的 bbox 尺寸为显示基准。
  - S2 以 `s2_rep1` 的 bbox 尺寸为显示基准。
- 每个细胞用真实轮廓填充表达量颜色，并用填充色的反色描边。
- 细胞轮廓线宽会随缩放动态变化：
  - 总览时轮廓变细，减少黑色/反色线条遮盖。
  - 放大后轮廓逐渐变粗，方便查看细胞边界。

## 数据文件

网页运行时使用服务端 SQLite 表达量索引和预处理轮廓 JSON：

```text
web_viewer/data/
  expression.sqlite
  S1_cells.json
  S2_cells.json
  genes.json
  replicates.json
  contours/
    S1/
      manifest.json
      tile_*.json
    S2/
      manifest.json
      tile_*.json
```

说明：

- `expression.sqlite`：运行时优先使用的表达量数据库，包含基因列表、每个样本/基因的表达范围、非零表达 cell，以及从 Rda 派生出的重复分组 JSON。也可包含 dotplot 表：`dotplot_clusters` 和 `dotplot_gene_cluster_stats`。
- `S1_cells.json` / `S2_cells.json`：每个 cell 的 bbox 和面积，是导出重复分组时使用的中间元数据。
- `replicates.json`：从 Seurat 对象的 `orig.ident` 导出，记录每个重复包含的 cell id、bbox 和需要加载的轮廓 tile；作为 SQLite 建库输入，数据库缺失时回退使用。
- `contours/*/manifest.json`：轮廓 tile 索引。
- `contours/*/tile_*.json`：真实细胞轮廓坐标。
- `genes.json`：旧版基因列表缓存；SQLite 缺失时回退使用。

服务端运行时不会读取 `.npy` 或 `.Rda`。
查询基因表达量时也不会扫描 `S1_all_genes.csv` / `S2_all_genes.csv`，除非 `expression.sqlite` 不存在或不可用。

## 离线导出

如果原始数据更新，按顺序重新生成预处理数据。

1. 导出 cell 空间元数据：

```bash
python3 web_viewer/export_cells.py \
  --mask npy/S1_cells.npy \
  --expr S1_all_genes.csv \
  --sample S1 \
  --out web_viewer/data/S1_cells.json

python3 web_viewer/export_cells.py \
  --mask npy/S2_cells.npy \
  --expr S2_all_genes.csv \
  --sample S2 \
  --out web_viewer/data/S2_cells.json
```

2. 导出真实细胞轮廓：

```bash
python3 web_viewer/export_contours.py \
  --mask npy/S1_cells.npy \
  --expr S1_all_genes.csv \
  --sample S1 \
  --out-dir web_viewer/data/contours/S1

python3 web_viewer/export_contours.py \
  --mask npy/S2_cells.npy \
  --expr S2_all_genes.csv \
  --sample S2 \
  --out-dir web_viewer/data/contours/S2
```

默认轮廓模式是 `cv2.CHAIN_APPROX_SIMPLE`。如果需要保留更多原始边界点，可以使用 `--chain none`，但前端绘制负担会增加。

3. 导出重复分组：

```bash
python3 web_viewer/export_replicates.py \
  --rda seurat_object.02-dims64.res0.6.Rda \
  --out web_viewer/data/replicates.json
```

该脚本读取 Seurat 对象 `st@meta.data$orig.ident`，生成 S1/S2 的重复分组。S2 中重复出现的数字 cell id 会按空间位置归属到一个重复，避免同一个 mask label 被重复绘制。

4. 构建 SQLite 表达量索引：

```bash
python3 web_viewer/build_expression_db.py \
  --s1-csv S1_all_genes.csv \
  --s2-csv S2_all_genes.csv \
  --replicates-json web_viewer/data/replicates.json \
  --out web_viewer/data/expression.sqlite \
  --replace
```

建库后，`/api/genes`、`/api/gene` 和 `/api/replicates` 会优先读取 SQLite。未建库时，服务端仍会回退到 `genes.json`、`replicates.json` 和 CSV/cache。

5. 导入 cluster dotplot 统计值：

```bash
python3 web_viewer/import_dotplot_stats.py \
  --rda seurat_object.02-dims64.res0.6.Rda \
  --db web_viewer/data/expression.sqlite \
  --cluster-column seurat_clusters \
  --expect-cells 60890 \
  --expect-clusters 23
```

该步骤会离线调用 R/Seurat 读取 `seurat_object.02-dims64.res0.6.Rda`，按 `st@meta.data$seurat_clusters` 分组，用当前 assay 的 `data` layer/slot 复刻 Seurat `DotPlot()` 的统计口径，然后写入 SQLite：

```sql
dotplot_clusters(
  cluster_id TEXT PRIMARY KEY,
  cluster_order INTEGER NOT NULL,
  cell_count INTEGER NOT NULL
)

dotplot_gene_cluster_stats(
  gene_id INTEGER NOT NULL,
  cluster_id TEXT NOT NULL,
  avg_expr REAL NOT NULL,
  pct_expr REAL NOT NULL,
  expressing_count INTEGER NOT NULL,
  cell_count INTEGER NOT NULL,
  PRIMARY KEY (gene_id, cluster_id)
)
```

导入完成后，网页运行时不会读取 `.Rda`，也不会调用 R/Rscript。

dotplot 统计口径：

- `avg_expr` 对应 Seurat `DotPlot()` 输出数据里的 `avg.exp`：
  `avg_expr = mean(expm1(data_slot_value))`，在当前 `gene × cluster` 的所有细胞内计算，包含 0 表达细胞。
- `pct_expr` 对应 Seurat `DotPlot()` 输出数据里的 `pct.exp * 100`：
  `pct_expr = count(data_slot_value > 0) / cell_count * 100`。
- `expressing_count` 是 `data_slot_value > 0` 的细胞数。
- 当前默认使用 Seurat 对象的默认 assay，本数据为 `SCT` assay 的 `data` layer/slot。

当前项目已完成建库：

- `web_viewer/data/expression.sqlite` 约 3.0 GB。
- 基因数：27006。
- S1：9993 个 cell，9,728,800 个非零表达值。
- S2：50897 个 cell，58,662,697 个非零表达值。

## 启动服务

```bash
python3 web_viewer/server.py --host 0.0.0.0 --port 8000
```

访问：

```text
http://127.0.0.1:8000/
```

## 查询与绘图逻辑

查询基因时，前端请求：

```text
/api/gene?gene=<gene_id>
/api/dotplot?gene=<gene_id>
```

后端会：

- 优先从 `web_viewer/data/expression.sqlite` 校验基因并读取 S1/S2 表达量。
- SQLite 中每个基因只保存非零表达 cell、表达范围和统计信息，查询时按 `(sample, gene)` 索引读取。
- 如果 SQLite 不存在或不可用，则回退到 CSV/cache 读取。
- 只返回非零表达 cell，未返回的 cell 在前端按 0 表达处理。
- 返回 S1/S2 合并后的表达范围，用于统一颜色归一化。
- CSV 回退路径会将查询结果缓存到 `web_viewer/cache/`。
- `/api/dotplot` 只读取 SQLite。缺少 dotplot 表时返回明确错误，不影响 `/api/gene`。
- `/api/dotplot` 返回 SQLite 中的 `avgExpr` 和 `pctExpr`，并额外返回 `avgExprScaled`。`avgExprScaled` 按 Seurat 默认 DotPlot 颜色口径计算：对当前基因各 cluster 的 `log1p(avgExpr)` 做 z-score，并 clamp 到 `[-2.5, 2.5]`。

前端会：

- 根据当前视口按需加载轮廓 tile。
- 根据 `/api/replicates` 返回的重复分组把 cell 放到对应重复面板。
- 将每个重复的原始 bbox 拉伸到同一样本的统一面板尺寸。
- 使用表达量填充真实细胞轮廓。
- 使用动态线宽绘制反色细胞轮廓。
- 在空间图下方绘制单行 cluster dotplot：X 轴为 `seurat_clusters`，Y 轴为当前查询基因。切换 S1/S2 只影响空间图，不会改变 dotplot。
- Dotplot 点颜色使用 `avgExprScaled`，对应 Seurat 默认的 scaled average expression。
- Dotplot 点大小使用当前基因的动态百分比范围映射，而不是固定 `0-100%`：
  当前基因所有 cluster 的最小 `pctExpr` 映射到最小半径，最大 `pctExpr` 映射到最大半径，中间值线性插值。这样低表达基因也能显示 cluster 间的相对差异。
- 如果当前基因所有 cluster 的 `pctExpr` 完全相同，点大小统一使用中等半径。

动态轮廓线宽的屏幕像素范围为：

```text
min: 0.15px
base at fitted view: 0.22px
max: 1.25px
```
