"""MkDocs 构建钩子（在 mkdocs.yml 的 `hooks:` 里注册）。

只做一件事：在构建产物目录里放一个 `.nojekyll`。

为什么需要：GitHub Pages 用 **Deploy from a branch → main / docs** 发布时，
会先过一遍 Jekyll。Jekyll 会跳过所有下划线开头的文件/目录，并可能改动其它文件。
放一个 `.nojekyll` 就让它原样输出。MkDocs 本身不生成这个文件
（`mkdocs gh-deploy` 是靠 ghp-import 的 --no-jekyll 顺带加的，而我们现在不用 gh-deploy）。
"""
from pathlib import Path


def on_post_build(config, **kwargs):
    site_dir = Path(config["site_dir"])
    (site_dir / ".nojekyll").write_text("", encoding="utf-8")
    print(f"已写入 {site_dir / '.nojekyll'}")
