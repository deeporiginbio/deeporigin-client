"""utility functions for working in Jupyter notebooks"""

from beartype import beartype
from IPython.display import HTML, display


@beartype
def show_progress_bar(
    *,
    completed: int,
    total: int,
    failed: int = 0,
    title: str = "Progress Report",
) -> None:
    """
    Displays a Bootstrap progress bar in a Jupyter Notebook with a title and external text.
    """

    progress_html = render_progress_bar(
        completed=completed,
        total=total,
        failed=failed,
        title=title,
    )

    display(HTML(progress_html))


@beartype
def render_progress_bar(
    *,
    completed: int,
    total: int,
    failed: int = 0,
    title: str = "Progress Report",
) -> str:
    """
    Displays a Bootstrap progress bar in a Jupyter Notebook with a title and external text.

    Parameters:
      completed (int): Total tasks attempted.
      total (int): Total tasks planned.
      failed (int): Number of failed tasks.
      title (str): Title to display above the progress bar.
    """
    if total <= 0:
        raise ValueError("Total must be a positive integer.")

    # Calculate passed (successful) and pending tasks.
    passed = max(completed - failed, 0)
    pending = max(total - completed, 0)

    # Calculate percentage for each segment relative to total.
    passed_pct = (passed / total) * 100
    failed_pct = (failed / total) * 100
    pending_pct = (pending / total) * 100

    # HTML for the title and text labels
    text_html = f"""
    <div style="margin-bottom: 5px;">
      <span>Completed: {passed}</span>
      <span style="margin-left: 15px;">Failed: {failed}</span>
      <span style="margin-left: 15px;">Remaining: {pending}</span>
    </div>
    """

    progress_html = f"""
    
    
    <h2>{title}</h2>
    {text_html}
    <div class="progress" style="height: 20px;">
      <div class="progress-bar bg-success" role="progressbar" style="width: {passed_pct:.1f}%"
           aria-valuenow="{passed}" aria-valuemin="0" aria-valuemax="{total}"></div>
      <div class="progress-bar bg-danger" role="progressbar" style="width: {failed_pct:.1f}%"
           aria-valuenow="{failed}" aria-valuemin="0" aria-valuemax="{total}"></div>
      <div class="progress-bar bg-secondary" role="progressbar" style="width: {pending_pct:.1f}%"
           aria-valuenow="{pending}" aria-valuemin="0" aria-valuemax="{total}"></div>
    </div>
    """

    return progress_html


def mermaid_to_html(diagram_code: str) -> str:
    """
    Converts a Mermaid diagram code to HTML.
    """

    html_code = f'<div class="mermaid">{diagram_code}</div>'

    html_code += """
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
    """

    return html_code


@beartype
def render_mermaid(diagram_code: str) -> None:
    """
    Renders a Mermaid diagram in a Jupyter Notebook cell.

    Parameters:
      diagram_code (str): The Mermaid diagram definition, e.g.,
        'graph TD; A-->B;'
    """

    # Check if mermaid is defined; if not, load it.
    # This snippet checks if window.mermaid exists, and if not, loads the script.
    display(HTML(mermaid_to_html(diagram_code)))
