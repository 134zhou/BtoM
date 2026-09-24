// MathJax 3 配置（MkDocs Material + pymdownx.arithmatex generic 模式的官方配方）。
//
// pymdownx.arithmatex 的 generic 模式会把 Markdown 里的 $...$ / $$...$$ 转成
// \(...\) / \[...\] 并包在 class="arithmatex" 的元素里，所以这里：
//   - inlineMath / displayMath 必须与之一致；
//   - options.processHtmlClass 只处理 arithmatex 生成的节点；
//   - tags: "ams" 让 \tag{...} 生效（本项目的理论推导里大量使用公式编号）。
window.MathJax = {
  tex: {
    inlineMath: [["\\(", "\\)"]],
    displayMath: [["\\[", "\\]"]],
    processEscapes: true,
    processEnvironments: true,
    tags: "ams"
  },
  options: {
    ignoreHtmlClass: ".*|",
    processHtmlClass: "arithmatex"
  }
};

// Material 的即时导航会替换页面内容，必须重新排版。
document$.subscribe(() => {
  MathJax.startup.output.clearCache();
  MathJax.typesetClear();
  MathJax.texReset();
  MathJax.typesetPromise();
});
