window.MathJax = {
  tex: {
    inlineMath: [['$', '$'], ['\\(', '\\)']],
    displayMath: [['$$', '$$'], ['\\[', '\\]']],
    processEscapes: true,
    processEnvironments: true
  },
  options: {
    // Process ALL HTML content including notebook cells
    ignoreHtmlClass: "no-mathjax",
    processHtmlClass: ".*"
  }
};

document$.subscribe(() => {
  MathJax.typesetPromise()
})
