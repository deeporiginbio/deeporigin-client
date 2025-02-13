"""utility functions for working in Jupyter notebooks"""

from beartype import beartype
from IPython.display import HTML, display


@beartype
def render_mermaid(diagram_code: str) -> None:
    """
    Renders a Mermaid diagram in a Jupyter Notebook cell.

    Parameters:
      diagram_code (str): The Mermaid diagram definition, e.g.,
        'graph TD; A-->B;'
    """
    # Create the HTML for the diagram.
    diagram_html = f'<div class="mermaid">{diagram_code}</div>'
    display(HTML(diagram_html))

    # Check if mermaid is defined; if not, load it.
    # This snippet checks if window.mermaid exists, and if not, loads the script.
    display(
        HTML("""
    <script>
      if (typeof mermaid === 'undefined') {
        var script = document.createElement('script');
        script.src = "https://cdn.jsdelivr.net/npm/mermaid/dist/mermaid.min.js";
        script.onload = function() {
          mermaid.initialize({startOnLoad:true});
          mermaid.init(undefined, document.getElementsByClassName("mermaid"));
        };
        document.head.appendChild(script);
      } else {
          mermaid.init(undefined, document.getElementsByClassName("mermaid"));
      }
    </script>
    """)
    )
