from typing import Optional

import pytest

from deeporigin.tools.job import Job


class _DummyTask:
    def __init__(self) -> None:
        self.cancelled = False

    def cancel(self) -> None:
        self.cancelled = True


@pytest.fixture()
def dummy_task(monkeypatch: pytest.MonkeyPatch) -> _DummyTask:
    task = _DummyTask()

    def _create_task(_coro):
        # Do not run the coroutine in tests; just return a dummy task
        return task

    monkeypatch.setattr("deeporigin.tools.job.asyncio.create_task", _create_task)
    return task


@pytest.fixture()
def capture_update_display(monkeypatch: pytest.MonkeyPatch):
    calls: list[dict] = []

    def _fake_update_display(html_obj, *, display_id: Optional[str] = None):
        calls.append(
            {
                "html": html_obj.data if hasattr(html_obj, "data") else str(html_obj),
                "display_id": display_id,
            }
        )

    monkeypatch.setattr("deeporigin.tools.job.update_display", _fake_update_display)
    return calls


def test_stop_watching_triggers_final_non_auto_render(
    monkeypatch: pytest.MonkeyPatch, dummy_task: _DummyTask, capture_update_display
):
    # Avoid network/template dependencies before instantiation
    monkeypatch.setattr(Job, "sync", lambda self: None)
    monkeypatch.setattr(
        Job,
        "_render_job_view",
        lambda self, *, will_auto_update=False: f"auto={will_auto_update}",
    )

    # Also patch display() to avoid notebook I/O
    monkeypatch.setattr("deeporigin.tools.job.display", lambda *_args, **_kwargs: None)

    # Arrange: Job with minimal required state
    job = Job.from_id("job-1")

    # Ensure there is at least one non-terminal status so watch() proceeds
    job._status = ["Running"]

    # Act: start watching and then stop
    job.watch()

    # display_id should be set (UUID string)
    assert job._display_id is not None and isinstance(job._display_id, str)

    job.stop_watching()

    # Assert: update_display was called with the final non-auto-updating render
    assert capture_update_display, "Expected at least one update_display call"
    last = capture_update_display[-1]
    # After stop_watching, _display_id is cleared; compare before clear by allowing None
    assert last["html"] == "auto=False"


def test_compose_error_overlay_html_contains_message(monkeypatch: pytest.MonkeyPatch):
    # Avoid side effects in __post_init__
    monkeypatch.setattr(Job, "sync", lambda self: None)
    job = Job.from_id("job-1")

    html = job._compose_error_overlay_html(message="boom")
    assert "boom" in html
    assert "Will retry automatically" in html
