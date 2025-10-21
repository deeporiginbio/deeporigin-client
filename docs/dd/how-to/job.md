# Job control

This document describes how to monitor, inspect, cancel and control jobs started on Deep Origin. 

## Job dataframes


### View all jobs 

To view all Jobs on Deep Origin, use:

```{.python notest}
from deeporigin.tools.job import get_dataframe

df = get_dataframe()
```
A dataframe with the following columns will be returned:

| Column Name | Description | 
| ---- | ---- |
| id | Job ID, e.g.: `9e9eb45e-93a6-4432-963e-669e982fde62` |
| created_at | Time stamp when Job request was received |
| execution_id | (internal) execution ID. | 
| completed_at | Time stamp of when the Job was completed |
| started_at | Time stamp of when the Job was started |
| status | One of `Succeeded` `Cancelled` `Failed` `Running` `Queued` `Created` `FailedQuotation` `Quoted` or `InsufficientFunds`|
| tool_key | Key of tool corresponding to this Job |
| tool_version | Version of tool |
| user_name | Name (or ID) of user that started this job |
| run_duration_minutes | Runtime in minutes, rounded to nearest minute |


### View all jobs by Status

Only jobs matching certain status(es) can be retrieved. For example,

```{.python notest}
from deeporigin.tools.job import get_dataframe

df = get_dataframe(only_with_status=["Running"])
```

only retrieves currently running jobs.

Multiple statuses can be retrieved using a single function call:

```{.python notest}
from deeporigin.tools.job import get_dataframe

df = get_dataframe(only_with_status=["Running", "Succeeded"])
```



### View Job metadata, inputs and outputs

By default, the job dataframe does not include information about metadata, inputs, and outputs. These can be included in the dataframe using:

```{.python notest}
from deeporigin.tools.job import get_dataframe

df = get_dataframe(
    include_metadata=True,
    include_inputs=True,
    include_outputs=True,
)
```

### Resolve user names

By default, user IDs are not resolved to names. To use user names in the dataframe, use:

```{.python notest}
from deeporigin.tools.job import get_dataframe

df = get_dataframe(
    resolve_user_names=True,
)
```

## The `Job` class

The `Job` class allows you to view jobs on Deep Origin, track and visualize their progress, and cancel them. 

### Constructing `Job` objects 


#### From Tools

Typically, tools will return `Job` objects when you run them. For example, starting a docking or ABFE run will return a job object:

```{.python notest}
# here, sim is a Complex
job = sim.docking.run(pocket=pockets[0])
```

#### From execution IDs

A `Job` can also be constructed from a single execution ID:

```{.python notest}
from deeporigin.tools.job import Job

# construct a job from a single ID
job = Job.from_id("job-id")
```

### Inspecting a Job

A job can be viewed by simply inspecting it:

```{.python notest}
job
```

A widget such as the following will be shown.

<iframe 
    src="../../images/job-docking.html" 
    width="100%" 
    height="400" 
    style="border:none;"
    title="Job widget"
></iframe>


??? danger "Jupyter notebook required"
    Inspecting jobs using the Job widget requires a Jupyter notebook (or lab). Other interactive environments, such as marimo, are not fully supported yet. 


### Progress monitoring 

To monitor the progress of a Job as it proceeds, use:

```{.python notest}
job.watch()
```

This causes the widget to auto-update on a fixed timer till the job reaches a terminal state. 

To to stop watching, use:


```{.python notest}
job.stop_watching()
```

### Job Status Display

When viewing a job, the status is displayed as a badge in the footer of the job widget:

- **Status Colors**: Each status has a distinct color:
  - `Running`: Blue (primary)
  - `Succeeded`: Green (success)  
  - `Failed`/`FailedQuotation`/`InsufficientFunds`: Red (danger)
  - `Cancelled`: Dark gray
  - `Created`: Gray (secondary)
  - `Queued`: Light blue (info)
  - `Quoted`: Yellow (warning)
  - Other statuses: Light gray

- **Auto-update Behavior**: The widget automatically stops updating when the job reaches a terminal state (`Succeeded`, `Failed`, `Cancelled`, `FailedQuotation`, `Quoted`, or `InsufficientFunds`)


### Cancelling jobs

To cancel a job, use:

```{python notest}
job.cancel()
```