

Variables and secrets specified in the Deep Orgin platform, either
at the organization level or the account level, can be installed onto
ComputeBenches. 

```bash
deeporigin variables install
```
Doing so prints a message that looks like:


```bash
No variables were modified

No variables were added
No variables were deleted
1 variables were unmodified:
  EnvironmentVariable: foo
```




!!! warning "Only on ComputeBenches"
    This functionality is meant to only function on Deep Origin ComputeBenches. 