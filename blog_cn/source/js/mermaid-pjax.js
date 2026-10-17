(function () {
  var MERMAID_URL = 'https://unpkg.com/mermaid/dist/mermaid.esm.min.mjs';
  var mermaidPromise = null;

  function loadMermaid() {
    if (!mermaidPromise) {
      mermaidPromise = import(MERMAID_URL).then(function (mod) {
        var mermaid = mod.default;
        mermaid.initialize({
          startOnLoad: true,
          flowchart: { curve: 'linear' }
        });
        return mermaid;
      });
    }
    return mermaidPromise;
  }

  function renderMermaid() {
    if (!document.querySelector('pre.mermaid:not([data-processed])')) {
      return;
    }
    loadMermaid()
      .then(function (mermaid) {
        var nodes = Array.prototype.slice.call(
          document.querySelectorAll('pre.mermaid:not([data-processed])')
        );
        if (nodes.length) {
          return mermaid.run({ nodes: nodes });
        }
      })
      .catch(function (error) {
        console.error('[mermaid] render failed', error);
      });
  }

  if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', renderMermaid);
  } else {
    renderMermaid();
  }

  document.addEventListener('pjax:complete', renderMermaid);
})();
