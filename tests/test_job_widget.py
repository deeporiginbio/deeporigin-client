"""Tests for Shadow DOM Job widget rendering."""

from types import SimpleNamespace

from deeporigin.tools.job import Job


def test_render_uses_shadow_dom(monkeypatch):
    """Ensure the rendered HTML contains the Shadow DOM custom element."""

    def fake_sync(self):
        # Minimal state to render
        self._attributes = [SimpleNamespace(startedAt=None, completedAt=None)]
        self._status = ["Succeeded"]
        self._progress_reports = ["{}"]
        self._resource_ids = ["res-1"]
        self._inputs = [{}]
        self._outputs = [{}]
        self._metadata = [{}]
        self._tool = [{"key": "Docking"}]

    monkeypatch.setattr(Job, "sync", fake_sync, raising=True)

    job = Job.from_id("exec-1")
    html = job._render_job_view()

    assert "<do-job-widget>" in html
