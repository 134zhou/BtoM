# 复现指南

## 1. Python 环境

| 用途 | 解释器 | 备注 |
|---|---|---|
| **反演 / 实验**（本仓库的算法部分） | `E:\Python\Miniforge\envs\discretize\python.exe` | 已装好全部依赖，**不需要新增任何包** |
| **站点构建**（`web/`） | `E:\Python\Miniforge\envs\HTML\python.exe` | 只有 mkdocs 相关依赖，见 `web/requirements.txt` |

算法环境（Python 3.14.5）：

| 包 | 版本 | 用途 |
|---|---|---|
| `numpy` | 2.4.6 | 全部数值 |
| `scipy` | 1.17.1 | `LinearOperator`、`lsqr`、`minimize(L-BFGS-B)`、`cKDTree`、`CubicSpline` |
| `discretize` | 0.12.0 | `TreeMesh` 八叉树自适应网格、面梯度、面对偶体积 |
| `numba` | 0.67.0 | `matmul` / `rmatvec` 并行核 |
| `matplotlib` | 3.10.9 | 出图 |
| `pyvista` | 0.48.4 | 读写 `.vtk` / `.vtu` |
| `pandas` | 3.0.3 | 读 CSV |

---

## 2. 目录结构

```
BtoM/
├── README.md                  仓库总览（本文件的精简版）
├── src/                       公用代码
│   ├── symmetry.py            对称面自动检测 + 镜像重合处理（库 + CLI + 自测）
│   ├── csv2vtk.py             CSV -> 结构化 VTK 点云
│   └── build_notebooks.py     由 .py（# %% 标记）生成 .ipynb
├── notebooks/                 全部 notebook 与对应的 .py 源
│   ├── B2M.ipynb              生产版：specialU 数据集（手工维护，无 .py 源）
│   ├── B2M_noreg.py/.ipynb    去掉正则化的对照
│   ├── B2M_newprob.py/.ipynb  新版小探头数据集 + 对称面自动检测
│   └── plot_2d_sections.py/.ipynb  文章配图
├── data/
│   ├── raw/                   ✅ 入库：原始测量数据
│   │   ├── specialU/  UshapeNormal/  UshapeNormal_New_prob/  legacy/
│   └── derived/               ❌ 不入库：反演产物（.vtu/.vtk），可由脚本重跑
├── figures/                   ❌ 不入库：论文/站点配图（由 notebooks 生成）
├── archive/                   历史文件：旧 notebook、早期试验、参数快照
├── docs/                      ⭐【构建产物】GitHub Pages 发布的就是这个目录
│                                 由 mkdocs build 生成，别手改，也别 .gitignore
└── web/                       站点工程
    ├── mkdocs.yml             站点配置（site_dir: ../docs）
    ├── requirements.txt
    └── docs/                  【站点源码】分页 markdown 平铺在这里
```

!!! note "路径无关"
    `notebooks/` 与 `src/` 里所有脚本/notebook 的开头都有一段**仓库路径引导**：
    从当前目录向上找到含 `.git` 的仓库根，之后所有数据/输出路径都相对它解析。
    所以**在仓库里的任何位置运行都能找到同一批文件**，不会再出现
    "换个目录跑就 FileNotFoundError"。

---

## 3. 跑反演

```bash
# 全部命令用 discretize 环境的解释器
PY=E:/Python/Miniforge/envs/discretize/python.exe

# ① 旧数据集（specialU）——生产版
$PY -m jupyter nbconvert --to notebook --execute notebooks/B2M.ipynb --inplace

# ② 新版小探头数据集（含对称面自动检测）
$PY notebooks/B2M_newprob.py          # 直接当脚本跑
# 或者：$PY -m jupyter nbconvert --to notebook --execute notebooks/B2M_newprob.ipynb --inplace

# ③ 无正则化对照
$PY notebooks/B2M_noreg.py

# ④ 文章配图
$PY notebooks/plot_2d_sections.py
```

产物：`data/derived/**/inverted_M_results_*.vtu`、`figures/*.png|pdf`。

> `notebooks/B2M_noreg.py`、`B2M_newprob.py`、`plot_2d_sections.py` 既是可直接执行的脚本，
> 也是对应 `.ipynb` 的**唯一事实来源**（`# %%` 分格）。改了 `.py` 之后，
> 用 `$PY src/build_notebooks.py` 重新生成 `.ipynb`（有 mtime 防呆，不会覆盖更新的 notebook）。
> 例外：`notebooks/B2M.ipynb` 是手工维护的，没有 `.py` 源。

---

## 4. 跑对称面检测

```bash
$PY src/symmetry.py --selftest                              # 合成数据自测，8/8 通过
$PY src/symmetry.py data/raw/UshapeNormal_New_prob           # 检测 + 报告 + 出图
$PY src/symmetry.py data/raw/UshapeNormal_New_prob --policy average
$PY src/symmetry.py data/raw/UshapeNormal_New_prob --regression   # 与旧公式回归对照
```

期望输出（关键几行）：`h_raw = 7.698 mm` → 吸附到 `h = 8.000 mm`，
`accepted = True`，重合 `18405`、补点 `30379`、最终 `79163` 点。

---

---

## 5. 构建与发布站点

```bash
# 依赖（首次/换机）——已装好则跳过
E:/Python/Miniforge/envs/HTML/python.exe -m pip install -r web/requirements.txt

# 本地预览
cd web
E:/Python/Miniforge/envs/HTML/python.exe -m mkdocs serve
# 浏览器打开 http://127.0.0.1:8000/

# 构建（--strict 会把断链/坏引用当错误）
# 产物直接写到【仓库根的 docs/】，因为 mkdocs.yml 里设了 site_dir: ../docs
E:/Python/Miniforge/envs/HTML/python.exe -m mkdocs build --strict
cd ..

# 发布：把源码和产物一起提交到 main（没有同步脚本，改了 md 直接构建）
git add web docs
git commit -m "docs: ..."
git push
```

GitHub 侧只需要设置一次：仓库 **Settings → Pages → Source** →
**Deploy from a branch** → 分支 **`main`** / 目录 **`/docs`**。
站点地址 `https://134zhou.github.io/BtoM/`。

!!! note "为什么用仓库根的 docs/"
    这是 GitHub Pages "Deploy from a branch" 唯一能选的两个目录之一
    （另一个是 `/`）。把构建产物直接输出到那里，就**不需要 `gh-pages` 分支、
    也不需要 `mkdocs gh-deploy`** —— 一个仓库、一条 `main` 分支搞定。

    两个 `docs/` 别搞混：

    | 路径 | 是什么 | 谁维护 |
    |---|---|---|
    | `web/docs/` | markdown **源文件** | 人写 |
    | `docs/` | 构建出来的 **HTML 产物** | `mkdocs build` 生成，**别手改** |

    产物里的 `docs/.nojekyll` 是**手工放的一个空文件**（已提交进仓库），
    MkDocs 构建时不会删掉它；作用是不让 Pages 走 Jekyll。

!!! warning "发布是手动动作"
    本仓库**没有**配置 GitHub Actions 自动部署。改完文档要自己跑一次
    `mkdocs build --strict` + `git add docs` 再提交，否则线上站点不会更新。
