# BtoM 文档站（MkDocs + Material）

BtoM 项目的静态文档站，用 **MkDocs + Material for MkDocs** 构建，公式由 **MathJax 3**
在浏览器端渲染（支持 `\frac`、`\sum`、`align`、`cases`、`\tag{}` 编号）。

**所有分页 markdown 都平铺在 `docs/` 这一个文件夹里**，导航在 `mkdocs.yml` 的 `nav`
里手工列出（与 FPD_blog 的做法一致）。

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

---

## 部署：一个仓库就够了

站点和代码在**同一个仓库**里，不需要另建仓库。三种发布方式里选的是**第三种**：

| 方式 | 做法 | 本项目 |
|---|---|---|
| ① Deploy from a branch | 在仓库 Settings → Pages 里选分支 + 目录，目录**只能选 `/` 或 `/docs`**（指**仓库根目录下**的 `docs/`） | ❌ 用不了：站点源码在 `web/docs/`，构建产物在 `web/site/`，都不是仓库根的 `docs/` |
| ② GitHub Actions | 写 workflow，检出仓库 → `mkdocs build` → 部署 | ⚪ 可选（本项目没配，按约定手动发布） |
| ③ **`gh-deploy` 推 `gh-pages` 分支** | `mkdocs gh-deploy` 用 ghp-import 把 `site/` 提交到**同一仓库**的 `gh-pages` 分支 | ✅ 采用 |

```bash
cd web
E:/Python/Miniforge/envs/HTML/python.exe -m mkdocs gh-deploy --force
```

然后在 GitHub 仓库 **Settings → Pages → Source** 选
**Deploy from a branch** → 分支 `gh-pages` / 目录 `(root)`。

- `gh-pages` 就是同一个仓库的一个分支，**不是第二个仓库**
- 站点地址：`https://134zhou.github.io/BtoM/`
  （仓库名不是 `<用户名>.github.io`，所以路径里带仓库名，`mkdocs.yml` 里的
  `site_url` 已按此填好，用于 sitemap/canonical）
- `gh-pages` 分支里的内容是构建产物，不用手改；改内容永远改 `main` 分支的 markdown

> 如果以后想省掉手动这一步，可以加一个 GitHub Actions workflow：
> 在 `web/` 下 `mkdocs build` 然后 `mkdocs gh-deploy`。目前按约定没配。

---

## 内容来源：`meshtest/*.md` 是唯一事实来源

站点里这三页**不是**手写的，而是仓库根文档的自动同步副本：

| 事实来源（仓库根） | 站点页面 |
|---|---|
| `meshtest/THEORY.md` | `docs/theory.md` |
| `meshtest/README.md` | `docs/experiments.md` |
| `meshtest/results.md` | `docs/experiment-results.md` |

```bash
# 改完 meshtest 的文档，同步到站点
E:/Python/Miniforge/envs/HTML/python.exe web/tools/sync_docs.py

# 只检查是否已同步（有漂移返回码 1，可以挂到提交前钩子）
E:/Python/Miniforge/envs/HTML/python.exe web/tools/sync_docs.py --check
```

**不要直接改 `docs/theory.md` / `docs/experiments.md` / `docs/experiment-results.md`** ——
下次同步会被覆盖。

其余页面（`index.md`、`b2m.md`、`symmetry.md`、`datasets.md`、`reproduce.md`）是手写的，
直接改站点里的文件即可。

## 图片：手工放置，不做同步

`sync_docs.py` **只管 markdown，不碰图片**。要加图就把文件丢进 `web/docs/`，
然后在页面里按同目录相对路径引用：

```markdown
![说明文字](文件名.png)
```

`mkdocs build --strict` 会检查图片是否存在，链接写错会直接构建失败。

候选图（源文件在 `figures/`、`meshtest/figures/`，都是 `.gitignore` 的产物目录，
由脚本生成）见站点里 `datasets.md` 的「配图由人工维护」一节。

## 目录

```
web/
  mkdocs.yml                  站点配置（主题、导航、arithmatex、MathJax 脚本）
  requirements.txt            依赖固定版本
  README.md                   本文件
  tools/
    sync_docs.py              把 meshtest 的三份 md 同步到 docs/（只同步 markdown）
  docs/                       ← 所有分页 markdown 都平铺在这里
    index.md                  首页（结论速查 + 记号约定）
    b2m.md                    反演方法：正演、网格、目标函数、λ 标定
    symmetry.md               对称面自动检测与镜像重合处理
    theory.md                 ← meshtest/THEORY.md（自动同步）
    experiments.md            ← meshtest/README.md（自动同步）
    experiment-results.md     ← meshtest/results.md（自动同步）
    datasets.md               数据集、几何、反演结果
    reproduce.md              环境、目录、全部复现命令
    javascripts/mathjax.js    MathJax 3 配置（Material + arithmatex 官方配方）
  site/                       构建产物（.gitignore 已忽略）
```

## 已知事项

- 公式在浏览器端由 MathJax 渲染，读者需要能访问 `unpkg.com`。要完全离线，
  把 `mkdocs.yml` 里 MathJax 的 URL 换成本地文件。
- `docs/` 下不放子文件夹（除 `javascripts/`）：`use_directory_urls: true` 时
  每页 URL 是 `/BtoM/<文件名>/`，同目录引用图片直接写文件名即可。
