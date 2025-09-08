# Job watching

## stop_watching()

Calling `job.stop_watching()` cancels the background update task and performs a final render with `will_auto_update=False` using the stored display id, so the UI immediately reflects a stopped state without the auto-update spinner.
