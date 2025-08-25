"""
Parallel execution utilities for running functions in parallel batches with retries.

This module provides two approaches to parallel execution:

1. run_func_in_parallel: Uses concurrent.futures.ThreadPoolExecutor for true parallelism
2. run_func_in_parallel_async: Uses asyncio for concurrent execution with cooperative multitasking

Comparison:
- ThreadPoolExecutor: True parallelism, good for CPU-bound tasks, simpler error handling
- asyncio: Better for I/O-bound tasks, lower memory overhead, more complex but more efficient
  for network/file operations
"""

from collections.abc import Callable
import concurrent.futures
import time
from typing import Any, Optional, TypeVar

from beartype import beartype

T = TypeVar("T")


@beartype
def run_func_in_parallel(
    *,
    func: Callable,
    batch_size: int = 10,
    max_retries: int = 3,
    sleep_between_batches: float = 0.1,
    args: list[dict[str, Any]],
) -> dict[str, Any]:
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
        """Call a function and track the time it takes."""

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
        for result, retries in zip(results, retries_left, strict=True)
    ):
        to_process = [
            i
            for i, (result, retries) in enumerate(
                zip(results, retries_left, strict=True)
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


def run_func_in_parallel_async(
    *,
    func: Callable,
    batch_size: int = 10,
    max_retries: int = 3,
    sleep_between_batches: float = 0.1,
    args: list[dict[str, Any]],
) -> dict[str, Any]:
    """
    Run a function in parallel batches with retries for failures using asyncio.

    This is an async version of the parallel execution utility. It processes
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
    import asyncio

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
    durations = [None] * total

    async def call_func_timed(idx: int) -> Optional[tuple]:
        """Call a function and track the time it takes. Async version."""
        nonlocal total_failures
        try:
            start = time.time()
            # Run the function in a thread pool to avoid blocking
            loop = asyncio.get_event_loop()
            # Create a lambda to unpack kwargs properly
            result = await loop.run_in_executor(None, lambda: func(**args[idx]))
            end = time.time()
            return (idx, result, end - start)
        except Exception as e:
            print(f"[Error] Call {idx} failed: {e}")
            total_failures += 1
            return (idx, None, None)

    async def process_batch(batch_indices):
        """Process a batch of indices and update results."""
        tasks = [call_func_timed(idx) for idx in batch_indices]
        batch_results = await asyncio.gather(*tasks)

        for idx, result, duration in batch_results:
            if result is not None:
                results[idx] = result
                durations[idx] = duration
            else:
                retries_left[idx] -= 1

    async def process_batches():
        """Process all batches until completion."""
        while any(
            (result is None and retries > 0)
            for result, retries in zip(results, retries_left, strict=True)
        ):
            to_process = [
                i
                for i, (result, retries) in enumerate(
                    zip(results, retries_left, strict=True)
                )
                if result is None and retries > 0
            ]

            for i in range(0, len(to_process), batch_size):
                batch = to_process[i : i + batch_size]
                await process_batch(batch)
                await asyncio.sleep(sleep_between_batches)

    async def main():
        """Main async processing loop."""
        start_time = time.time()
        await process_batches()
        return time.time() - start_time

    # Run the async main function
    elapsed_time = asyncio.run(main())
    permanent_failures = [i for i, result in enumerate(results) if result is None]

    return {
        "results": results,
        "total_failures": total_failures,
        "permanent_failures": permanent_failures,
        "elapsed_time": elapsed_time,
        "durations": durations,
    }
