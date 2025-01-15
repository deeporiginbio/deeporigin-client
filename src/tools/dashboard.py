import json
import os

import pandas as pd
import streamlit as st
from box import Box
from deeporigin.tools.utils import query_run_status
from deeporigin.utils.core import _ensure_do_folder

# Define the cache directory
cache_dir = _ensure_do_folder() / "jobs"


def read_jobs() -> list:
    # Initialize an empty list to store job data
    jobs = []

    # Ensure the directory exists
    if os.path.exists(cache_dir):
        # Iterate over all JSON files in the directory
        for file_name in os.listdir(cache_dir):
            if file_name.endswith(".json"):
                file_path = os.path.join(cache_dir, file_name)
                try:
                    # Load JSON content and append to the jobs list
                    with open(file_path, "r") as f:
                        job_data = json.load(f)
                        jobs.append(Box(job_data))
                except Exception as e:
                    print(f"Failed to load {file_path}: {e}")
    else:
        print(f"Directory {cache_dir} does not exist.")

    print(f"Loaded {len(jobs)} jobs.")

    return jobs


def update_all_jobs(jobs: list):
    for job in jobs:
        if job.attributes.status == "Running" or "Queued" or "Created":
            query_run_status(job.id)


def get_data():
    jobs = read_jobs()
    update_all_jobs(jobs)
    jobs = read_jobs()
    return pd.DataFrame(
        {
            "Job Id": [job.id for job in jobs],
            "Execution ID": [job.attributes.executionId for job in jobs],
            "Status": [job.attributes.status for job in jobs],
            "Tool": [job.attributes.toolId for job in jobs],
        }
    )


st.title("Deep Origin Tools Dashboard")


config = {
    "_index": st.column_config.TextColumn("Job ID"),
}


@st.fragment(func=None, run_every=5)
def show_dataframe():
    df = get_data()
    st.dataframe(df, column_config=config, hide_index=True)


show_dataframe()

# time.sleep(5)

# st.rerun()
