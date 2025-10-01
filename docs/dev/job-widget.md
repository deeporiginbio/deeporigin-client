# Job widget (Shadow DOM)

The job status widget is rendered via a single Shadow DOM-based template (`src/templates/job_widget.html`). This isolates CSS/JS so it works consistently across JupyterLab, classic Notebook, VS Code notebooks, marimo, and web front-ends.

Key details:

- Uses a custom element `do-job-widget` which attaches an open shadow root.
- Loads Bootstrap 5 (CSS) inside the shadow via CDN, avoiding conflicts with legacy Bootstrap shipped by notebooks.
- Tabs are implemented with a small vanilla JavaScript controller inside the shadow; JavaScript is required.
- The Python code in `src/tools/job.py` renders this template in all environments.

CSP and offline notes:

- By default, assets load from `https://cdn.jsdelivr.net`. If your environment has a strict CSP or is offline, vendor these assets and update the `<link>` in the template accordingly.

Usage:

- In Python: `Job.from_id("<execution_id>").show()` or `.watch()`.


