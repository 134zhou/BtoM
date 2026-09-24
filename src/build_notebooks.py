# -*- coding: utf-8 -*-
"""把带 # %% 标记的 .py 源文件转换为 .ipynb（仅用 stdlib，无需 nbformat）。"""
import json
import pathlib


def build(py_path, nb_path):
    py = pathlib.Path(py_path)
    nb_file = pathlib.Path(nb_path)
    # 防呆：目标 notebook 比源 .py 新时不重建，避免误覆盖手工维护过的 notebook
    if nb_file.exists() and nb_file.stat().st_mtime >= py.stat().st_mtime:
        print(f"skip {nb_path} (up to date)")
        return
    lines = py.read_text(encoding="utf-8").splitlines()
    cells = []
    cur_type = None
    cur_lines = []

    def flush():
        nonlocal cur_type, cur_lines
        if cur_type is not None:
            if cur_type == "markdown":
                md = []
                for ln in cur_lines:
                    if ln.startswith("# "):
                        md.append(ln[2:])
                    elif ln.startswith("#"):
                        md.append(ln[1:])
                    else:
                        md.append(ln)
                src = "\n".join(md).rstrip("\n")
                cells.append({"cell_type": "markdown", "metadata": {}, "source": src + "\n"})
            else:
                src = "\n".join(cur_lines).rstrip("\n")
                cells.append({"cell_type": "code", "execution_count": None,
                              "metadata": {}, "outputs": [], "source": src + "\n"})
            cur_lines = []

    for ln in lines:
        if ln.startswith("# %% [markdown]"):
            flush()
            cur_type = "markdown"
            cur_lines = []
        elif ln.startswith("# %%"):
            flush()
            cur_type = "code"
            cur_lines = []
        else:
            if cur_type is not None:
                cur_lines.append(ln)
    flush()

    nb = {
        "cells": cells,
        "metadata": {
            "kernelspec": {"display_name": "discretize", "language": "python", "name": "python3"},
            "language_info": {"name": "python", "version": "3.14.5"},
        },
        "nbformat": 4,
        "nbformat_minor": 5,
    }
    pathlib.Path(nb_path).write_text(json.dumps(nb, ensure_ascii=False, indent=1), encoding="utf-8")
    print(f"built {nb_path} ({len(cells)} cells)")


if __name__ == "__main__":
    NB = pathlib.Path(__file__).resolve().parent.parent / "notebooks"
    build(NB / "plot_2d_sections.py", NB / "plot_2d_sections.ipynb")
    build(NB / "B2M_noreg.py", NB / "B2M_noreg.ipynb")
    build(NB / "B2M_newprob.py", NB / "B2M_newprob.ipynb")
