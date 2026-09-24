# -*- coding: utf-8 -*-
"""把仓库里的「事实来源」文档与图片同步到 MkDocs 站点目录 web/docs/。

事实来源（相对仓库根）          站点页面
---------------------------  ------------------------------------
meshtest/THEORY.md       ->  docs/theory/index.md
meshtest/README.md       ->  docs/experiments/index.md
meshtest/results.md      ->  docs/experiments/results.md
figures/*.png            ->  docs/assets/*.png
meshtest/figures/*.png   ->  docs/assets/*.png

设计原则：**仓库根的文档是唯一事实来源**，站点里这几页只是副本。
所有内容改动都发生在源文件里，改完重新同步即可，不需要手改两份。

用法（解释器 E:\\Python\\Miniforge\\envs\\HTML）：
    python web/tools/sync_docs.py              # 同步（只有内容真的变了才写盘）
    python web/tools/sync_docs.py --check      # 只检查：有漂移则返回 1（可挂 CI / 提交前钩子）
    python web/tools/sync_docs.py --no-assets  # 只同步 markdown，不动图片
"""
from __future__ import annotations

import argparse
import filecmp
import shutil
import sys
from pathlib import Path

WEB = Path(__file__).resolve().parent.parent      # web/
REPO = WEB.parent                                 # 仓库根
DOCS = WEB / "docs"
REPO_URL = "https://github.com/134zhou/BtoM/blob/main/"

# 源文件（相对仓库根） -> 站点页面（相对 web/）
MD_MAP = {
    "meshtest/THEORY.md": "docs/theory/index.md",
    "meshtest/README.md": "docs/experiments/index.md",
    "meshtest/results.md": "docs/experiments/results.md",
}

# 只把站点真正引用到的图拷进来（figures/ 与 meshtest/figures/ 都是 .gitignore 的产物目录）
ASSET_MAP = {
    "figures/symmetry_plane_detection.png": "docs/assets/symmetry_plane_detection.png",
    "figures/mirror_overlap.png": "docs/assets/mirror_overlap.png",
    "figures/B_field_section.png": "docs/assets/B_field_section.png",
    "figures/M_section.png": "docs/assets/M_section.png",
    "figures/M_Uspecial.png": "docs/assets/M_Uspecial.png",
    "figures/M_regularized_vs_noreg.png": "docs/assets/M_regularized_vs_noreg.png",
    "meshtest/figures/E1_discretization.png": "docs/assets/E1_discretization.png",
    "meshtest/figures/E2_lambda_sweep.png": "docs/assets/E2_lambda_sweep.png",
    "meshtest/figures/E5_lambda_rule.png": "docs/assets/E5_lambda_rule.png",
}


def banner(src_rel: str) -> str:
    url = REPO_URL + src_rel
    return (
        '!!! abstract "本页是自动同步的副本"\n'
        f'    内容来自仓库中的 [`{src_rel}`]({url})，由 `web/tools/sync_docs.py` 生成，'
        "**请不要直接改这一页**。\n"
        "    要改内容请改源文件，然后在 `web/` 目录执行：\n"
        "\n"
        "    ```bash\n"
        "    python tools/sync_docs.py\n"
        "    ```\n"
    )


def render(src_rel: str, text: str) -> str:
    """在标题之后插入同步说明；正文一字不改。"""
    text = text.replace("\r\n", "\n").rstrip("\n") + "\n"
    lines = text.split("\n")
    if lines and lines[0].startswith("# "):
        head, rest = lines[0], "\n".join(lines[1:]).lstrip("\n")
        out = f"{head}\n\n{banner(src_rel)}\n{rest}"
    else:
        out = f"{banner(src_rel)}\n{text}"
    return out.rstrip("\n") + "\n"


def read_text(p: Path) -> str:
    with open(p, encoding="utf-8", newline="") as f:
        return f.read()


def write_text(p: Path, s: str) -> None:
    p.parent.mkdir(parents=True, exist_ok=True)
    with open(p, "w", encoding="utf-8", newline="\n") as f:
        f.write(s)


def sync_md(check: bool) -> int:
    drift = 0
    for src_rel, dst_rel in MD_MAP.items():
        src, dst = REPO / src_rel, WEB / dst_rel
        if not src.is_file():
            print(f"  !! 源文件不存在：{src_rel}")
            drift += 1
            continue
        want = render(src_rel, read_text(src))
        have = read_text(dst) if dst.is_file() else None
        if have == want:
            print(f"  == 已最新  {dst_rel}  <-  {src_rel}")
        elif check:
            print(f"  !! 需要同步 {dst_rel}  <-  {src_rel}")
            drift += 1
        else:
            write_text(dst, want)
            print(f"  -> 已同步  {dst_rel}  <-  {src_rel}")
    return drift


def sync_assets(check: bool) -> int:
    drift = 0
    for src_rel, dst_rel in ASSET_MAP.items():
        src, dst = REPO / src_rel, WEB / dst_rel
        if not src.is_file():
            print(f"  !! 缺图（先跑出图脚本）：{src_rel}")
            drift += 1
            continue
        if dst.is_file() and filecmp.cmp(src, dst, shallow=False):
            print(f"  == 已最新  {dst_rel}")
        elif check:
            print(f"  !! 需要拷贝 {dst_rel}  <-  {src_rel}")
            drift += 1
        else:
            dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(src, dst)
            print(f"  -> 已拷贝  {dst_rel}  <-  {src_rel}")
    return drift


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="同步仓库文档/图片到 MkDocs 站点")
    ap.add_argument("--check", action="store_true", help="只检查，不写盘；有漂移返回 1")
    ap.add_argument("--no-assets", action="store_true", help="只同步 markdown")
    args = ap.parse_args(argv)

    print(f"仓库根：{REPO}")
    print(f"站点  ：{WEB}")
    print("\n[markdown]")
    drift = sync_md(args.check)
    if not args.no_assets:
        print("\n[图片]")
        drift += sync_assets(args.check)

    if drift and args.check:
        print(f"\n有 {drift} 项与源文件不一致，请执行：python web/tools/sync_docs.py")
        return 1
    print("\n完成。")
    return 0


if __name__ == "__main__":
    sys.exit(main())
