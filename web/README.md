# BtoM 文档站（MkDocs + Material）

BtoM 项目的静态文档站，用 **MkDocs + Material for MkDocs** 构建，公式由 **MathJax 3**
在浏览器端渲染（支持 `\frac`、`\sum`、`align`、`cases`、`\tag{}` 编号）。

**所有分页 markdown 都平铺在 `docs/` 这一个文件夹里**，导航在 `mkdocs.yml` 的 `nav`
里手工列出（与 FPD_blog 的做法一致）。页面全是手写的，没有任何同步脚本。

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
# 产物直接写到【仓库根的 docs/】（mkdocs.yml 里 site_dir: ../docs）
# --strict 把断链/坏引用当错误
```

---

## 部署：一个仓库 + 仓库根的 `docs/` 文件夹

站点和代码在**同一个仓库**里。构建产物写在仓库根的 `docs/`，然后把它一起提交到 `main`，
GitHub Pages 直接从那里发布 —— 不用 `gh-pages` 分支，也不用 `mkdocs gh-deploy`。

| 方式 | 做法 | 本项目 |
|---|---|---|
| **① Deploy from a branch** | Settings → Pages → Source 选 `main` + `/docs` | ✅ **采用** |
| ② GitHub Actions | workflow 里构建 + 部署 | ⚪ 可选，没配 |
| ③ `gh-deploy` 推 `gh-pages` | `mkdocs gh-deploy` | ❌ 不用（会多出一个分支） |

完整流程：

```bash
# 1) 改文档：web/docs/*.md、mkdocs.yml
# 2) 重新构建（产物进仓库根的 docs/）
cd web
E:/Python/Miniforge/envs/HTML/python.exe -m mkdocs build --strict
cd ..

# 3) 把源码和产物一起提交
git add web docs
git commit -m "docs: ..."
git push
```

GitHub 侧设置一次即可：仓库 **Settings → Pages → Source** →
**Deploy from a branch** → 分支 **`main`** / 目录 **`/docs`**。

- 站点地址：`https://134zhou.github.io/BtoM/`
  （仓库名不是 `<用户名>.github.io`，所以路径里带仓库名；`mkdocs.yml` 的 `site_url`
  已按此填好，用于 sitemap/canonical）
- 仓库根的 `docs/` 是**构建产物**，**不要手改**，也不要加进 `.gitignore`
- 仓库根 `docs/.nojekyll` 是一个**手工放的空文件**（已提交）。它让 Pages 不走 Jekyll，
  否则下划线开头的文件会被跳过。MkDocs 构建时不会删掉它，不用管。

> ⚠️ 代价：`docs/` 是生成物却要入库，所以**改完文档必须记得重新 `mkdocs build`
> 再提交**，否则线上还是旧内容（`mkdocs build` 不会自动 commit）。

---

## 图片：手工放置

要加图就把文件丢进 `web/docs/`，然后在页面里按同目录相对路径引用：

```markdown
![说明文字](文件名.png)
```

`mkdocs build --strict` 会检查图片是否存在，链接写错会直接构建失败。

候选图（源文件在 `figures/`，是 `.gitignore` 的产物目录，由 notebook 脚本生成）
见站点里 `datasets.md` 的「配图由人工维护」一节。

## 目录

```
BtoM/
├── docs/                      ← 【构建产物】GitHub Pages 就是发布这个目录
│                                 自动生成，不要手改，也不要 .gitignore
└── web/                       站点工程
    ├── mkdocs.yml             站点配置（主题、导航、arithmatex、MathJax）
    ├── requirements.txt       依赖固定版本
    ├── README.md              本文件
    └── docs/                  ← 【站点源码】所有分页 markdown 都平铺在这里
        ├── index.md           首页（结论速查 + 记号约定）
        ├── b2m.md             反演方法：正演、网格、目标函数、λ 标定
        ├── symmetry.md        对称面自动检测与镜像重合处理
        ├── datasets.md        数据集、几何、反演结果
        ├── reproduce.md       环境、目录、全部复现命令
        └── javascripts/mathjax.js  MathJax 3 配置（Material + arithmatex 官方配方）
```

## 已知事项

- 公式在浏览器端由 MathJax 渲染，读者需要能访问 `unpkg.com`。要完全离线，
  把 `mkdocs.yml` 里 MathJax 的 URL 换成本地文件。
- **构建产物要提交**：`docs/` 是生成物却进了版本库，这是"用 main 分支的 `/docs`
  发布"这个方案的必要代价。
- `web/docs/` 下不放子文件夹（除 `javascripts/`）：`use_directory_urls: true` 时
  每页 URL 是 `/BtoM/<文件名>/`，同目录引用图片直接写文件名即可。
- 两个 `docs/` 别搞混：`web/docs/` 是 **markdown 源**，
  仓库根的 `docs/` 是 **HTML 产物**。
