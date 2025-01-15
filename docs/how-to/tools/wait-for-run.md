
When you start a tool run, you will receive a job ID. This job ID can be used to wait for the tool run to finish. This is useful if you want to block further code execution until the tool run is finished.

To wait for the tool run to finish, use:

```python
from deeporigin.tools.utils import wait_for_job
wait_for_job("job-id")
```
