"""
Parallel execution utilities for running functions in parallel batches with retries.
"""

import concurrent.futures
import time
from typing import Any, Callable, Dict, List, Optional, TypeVar

from beartype import beartype

T = TypeVar("T")


@beartype
def run_func_in_parallel(
    *,
    func: Callable[..., T],
    batch_size: int = 10,
    max_retries: int = 3,
    sleep_between_batches: float = 0.1,
    args: List[Dict[str, Any]],
) -> Dict[str, Any]:
    """
    Run a function in parallel batches with retries for failures.

    This is a general wrapper that can parallelize any function call. It processes
    items in batches, retries failed calls, and provides comprehensive results.

    Args:
        func: The function to call in parallel. Must accept **kwargs.
        batch_size: Number of items to process in parallel per batch.
        max_retries: Maximum number of times to retry failed calls.
        sleep_between_batches: Delay between batches, in seconds.
        args: List of argument dictionaries. Each dict will be passed as **kwargs to func.

    Returns:
        Dict containing:
            - results: List of results from function calls (or None for failures)
            - total_failures: Total number of failures encountered (including retries)
            - permanent_failures: Indices of items that failed even after retries
            - elapsed_time: Time in seconds for entire process
            - durations: List of time spent on each call (or None if permanently failed)


    """
    if not args:
        return {
            "results": [],
            "total_failures": 0,
            "permanent_failures": [],
            "elapsed_time": 0,
            "durations": [],
        }

    total = len(args)
    results = [None] * total
    retries_left = [max_retries] * total
    total_failures = 0
    durations = [None] * total  # Track time spent per call

    def call_func_timed(idx: int) -> Optional[tuple]:
        nonlocal total_failures
        try:
            start = time.time()
            result = func(**args[idx])
            end = time.time()
            return (idx, result, end - start)
        except Exception as e:
            print(f"[Error] Call {idx} failed: {e}")
            total_failures += 1
            return (idx, None, None)

    start_time = time.time()

    def process_batch(batch_indices):
        """Process a batch of indices and update results."""
        with concurrent.futures.ThreadPoolExecutor() as executor:
            batch_results = list(executor.map(call_func_timed, batch_indices))

        for idx, result, duration in batch_results:
            if result is not None:
                results[idx] = result
                durations[idx] = duration
            else:
                retries_left[idx] -= 1

    # Process in batches until all successful or max retries exhausted
    while any(
        (result is None and retries > 0)
        for result, retries in zip(results, retries_left, strict=False)
    ):
        to_process = [
            i
            for i, (result, retries) in enumerate(
                zip(results, retries_left, strict=False)
            )
            if result is None and retries > 0
        ]

        # Process in batches
        for i in range(0, len(to_process), batch_size):
            batch = to_process[i : i + batch_size]
            process_batch(batch)
            time.sleep(sleep_between_batches)

    elapsed_time = time.time() - start_time
    permanent_failures = [i for i, result in enumerate(results) if result is None]

    return {
        "results": results,
        "total_failures": total_failures,
        "permanent_failures": permanent_failures,
        "elapsed_time": elapsed_time,
        "durations": durations,
    }
