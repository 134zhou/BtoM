# BtoM 文档站（MkDocs + Material）

BtoM 项目的静态文档站，用 **MkDocs + Material for MkDocs** 构建，公式由 **MathJax 3**
在浏览器端渲染（支持 `\frac`、`\sum`、`align`、`cases`、`\tag{}` 编号）。

## 环境

- 解释器：`E:\Python\Miniforge\envs\HTML`（Python 3.13）
- 依赖：见 `requirements.txt`（mkdocs-material 9.7.7 / mkdocs 1.6.1 / pymdown-extensions 11.0.2）

```bash
# 首次或换机重建环境
E:/Python/Miniforge/envs/HTML/python.exe -m pip install -r requirements.txt
```

## 本地预览

```bash
cd web
E:/Python/Miniforge/envs/HTML/python.exe -m mkdocs serve
# 浏览器打开 http://127.0.0.1:8000/
```

## 构建

```bash
cd web
E:/Python/Miniforge/envs/HTML/python.exe -m mkdocs build --strict
# 产物在 web/site/，--strict 把断链/坏引用当错误
```

## 发布到 GitHub Pages

`mkdocs gh-deploy` 用 ghp-import 把 `site/` 的内容提交到仓库的 **`gh-pages` 分支**
（不污染 `main`），然后在 GitHub 仓库 Settings → Pages 里把 Source 设为
`gh-pages` / `(root)` 即可。

```bash
cd web
E:/Python/Miniforge/envs/HTML/python.exe -m mkdocs gh-deploy --force
```

站点地址：`https://134zhou.github.io/BtoM/`（仓库名不是 `<用户名>.github.io`，
所以路径里带仓库名，`mkdocs.yml` 里的 `site_url` 已按此填好）。

> 本仓库**没有**配置自动部署的 GitHub Actions，发布是手动动作。

## 内容来源：`meshtest/*.md` 是唯一事实来源

站点里 `理论推导` / `对照实验` 三页**不是**手写的，而是仓库根文档的自动同步副本：

| 事实来源（仓库根） | 站点页面 |
|---|---|
| `meshtest/THEORY.md` | `docs/theory/index.md` |
| `meshtest/README.md` | `docs/experiments/index.md` |
| `meshtest/results.md` | `docs/experiments/results.md` |
| `figures/*.png`、`meshtest/figures/*.png` | `docs/assets/*.png` |

```bash
# 改完 meshtest 的文档/图，同步到站点
E:/Python/Miniforge/envs/HTML/python.exe web/tools/sync_docs.py

# 只检查是否已同步（有漂移返回码 1，可以挂到提交前钩子）
E:/Python/Miniforge/envs/HTML/python.exe web/tools/sync_docs.py --check
```

**不要直接改 `docs/theory/` 与 `docs/experiments/` 下的文件** —— 下次同步会被覆盖。

其余页面（`index.md`、`methods/*`、`results/new-prob.md`、`reproduce.md`）是手写的，
直接改站点里的文件即可。

## 目录

```
web/
  mkdocs.yml                  站点配置（主题、导航、arithmatex、MathJax 脚本）
  requirements.txt            依赖固定版本
  README.md                   本文件
  tools/
    sync_docs.py              同步 meshtest 文档与图片到 docs/
  docs/
    index.md                  首页（结论速查 + 记号约定）
    javascripts/mathjax.js    MathJax 3 配置（Material + arithmatex 官方配方）
    methods/b2m.md            反演方法：正演、网格、目标函数、λ 标定
    methods/symmetry.md       对称面自动检测与镜像重合处理
    theory/index.md           ← meshtest/THEORY.md（自动同步）
    experiments/index.md      ← meshtest/README.md（自动同步）
    experiments/results.md    ← meshtest/results.md（自动同步）
    results/new-prob.md       数据集、几何、反演结果图
    reproduce.md              环境、目录、全部复现命令
    assets/                   配图（由 sync_docs.py 拷贝）
  site/                       构建产物（.gitignore 已忽略）
```

## 已知事项

- 公式在浏览器端由 MathJax 渲染，读者需要能访问 `unpkg.com`。要完全离线，
  把 `mkdocs.yml` 里 MathJax 的 URL 换成本地文件。
- `docs/assets/` 里的图是从 `figures/`（`.gitignore`）拷来的。这些图已经入库，
  所以 clone 下来就能构建；但 `sync_docs.py --check` 在缺图的机器上会报"需要拷贝"。
