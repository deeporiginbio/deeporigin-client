"""Recording utilities for platform client interactions.

Provides a lightweight SQLite-backed recorder that can persist successful
request/response pairs for later replay in tests via a mock client.
"""

from __future__ import annotations

from dataclasses import asdict, is_dataclass
import hashlib
import json
from pathlib import Path
import sqlite3
import threading
import time
from typing import Any

DEFAULT_DB_PATH = (
    Path(__file__).parent.parent.parent / "tests" / "fixtures" / "requests.sqlite"
)


def _ensure_parent_dir(path: Path) -> None:
    """Ensure parent directory exists for the provided file path."""

    path.parent.mkdir(parents=True, exist_ok=True)


def _json_canonical_dumps(value: Any) -> str:
    """Serialize a Python object to a canonical JSON string.

    - Dataclasses converted via asdict
    - Sets converted to sorted lists
    - Path/UUID/datetime converted to string via str()
    - Non-JSON-serializable objects attempted via __dict__ then str()
    """

    def normalize(obj: Any) -> Any:
        if is_dataclass(obj):
            return normalize(asdict(obj))
        if isinstance(obj, dict):
            return {
                str(k): normalize(v)
                for k, v in sorted(obj.items(), key=lambda i: str(i[0]))
            }
        if isinstance(obj, (list, tuple)):
            return [normalize(i) for i in obj]
        if isinstance(obj, set):
            return [normalize(i) for i in sorted(obj, key=lambda x: str(x))]
        if isinstance(obj, (Path,)):
            return str(obj)
        try:
            json.dumps(obj)
            return obj
        except TypeError:
            if hasattr(obj, "to_dict"):
                return normalize(obj.to_dict())
            if hasattr(obj, "__dict__"):
                return normalize(obj.__dict__)
            return str(obj)

    normalized = normalize(value)
    return json.dumps(normalized, sort_keys=True, separators=(",", ":"))


def compute_request_hash(method: str, kwargs: dict) -> str:
    """Compute a stable SHA-256 hash for a request based on method and kwargs."""

    payload = f"{method}\n{_json_canonical_dumps(kwargs)}".encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


class RequestRecorder:
    """SQLite-backed recorder for request/response interactions.

    Thread-safe; maintains a single connection. Designed for tests and lightweight usage.
    """

    def __init__(self, db_path: Path | None = None):
        self.db_path = Path(db_path) if db_path else DEFAULT_DB_PATH
        _ensure_parent_dir(self.db_path)
        # Use a single connection; enable WAL for better concurrency
        self._conn = sqlite3.connect(str(self.db_path), check_same_thread=False)
        self._conn.execute("PRAGMA journal_mode=WAL;")
        self._conn.execute("PRAGMA synchronous=NORMAL;")
        self._lock = threading.Lock()
        self._init_schema()

    def _init_schema(self) -> None:
        with self._conn:  # implicit transaction
            self._conn.execute(
                """
                CREATE TABLE IF NOT EXISTS interactions (
                    id INTEGER PRIMARY KEY,
                    timestamp TEXT NOT NULL,
                    method TEXT NOT NULL,
                    request_json TEXT NOT NULL,
                    request_hash TEXT NOT NULL,
                    sequence_num INTEGER NOT NULL,
                    response_json TEXT NOT NULL,
                    duration_ms INTEGER NOT NULL
                );
                """
            )
            self._conn.execute(
                """
                CREATE UNIQUE INDEX IF NOT EXISTS idx_interactions_unique
                ON interactions(request_hash, sequence_num);
                """
            )
            self._conn.execute(
                """
                CREATE INDEX IF NOT EXISTS idx_interactions_method
                ON interactions(method);
                """
            )

    def record(
        self,
        *,
        method: str,
        kwargs: dict,
        response: Any,
        duration_ms: int,
        db_path: Path | None = None,
    ) -> None:
        """Record a successful interaction.

        Args:
            method: Fully qualified method name.
            kwargs: Request arguments.
            response: Successful response object (JSON-serializable or encodable).
            duration_ms: Duration of the call in milliseconds.
            db_path: Optional override DB path for this call.
        """

        if db_path and Path(db_path) != self.db_path:
            # If an override is provided, write to that DB lazily via a separate recorder
            override = RequestRecorder(Path(db_path))
            override.record(
                method=method, kwargs=kwargs, response=response, duration_ms=duration_ms
            )
            return

        request_json = _json_canonical_dumps(kwargs)
        request_hash = compute_request_hash(method, kwargs)
        response_json = _json_canonical_dumps(response)
        now_iso = time.strftime("%Y-%m-%dT%H:%M:%S", time.gmtime())

        with self._lock:
            cur = self._conn.cursor()
            cur.execute(
                "SELECT COALESCE(MAX(sequence_num), -1) + 1 FROM interactions WHERE request_hash=?",
                (request_hash,),
            )
            next_seq = cur.fetchone()[0]
            cur.execute(
                """
                INSERT INTO interactions(timestamp, method, request_json, request_hash, sequence_num, response_json, duration_ms)
                VALUES(?, ?, ?, ?, ?, ?, ?)
                """,
                (
                    now_iso,
                    method,
                    request_json,
                    request_hash,
                    next_seq,
                    response_json,
                    int(duration_ms),
                ),
            )
            self._conn.commit()

    def fetch_next(
        self,
        *,
        method: str,
        kwargs: dict,
        state: dict[str, int] | None = None,
    ) -> Any:
        """Fetch the next recorded response for the given request.

        Maintains an optional mutable state mapping request_hash->next_seq to achieve
        deterministic sequencing during tests.
        """

        request_hash = compute_request_hash(method, kwargs)
        with self._lock:
            cur = self._conn.cursor()
            # Determine next sequence to return
            next_seq = 0
            if state is not None:
                next_seq = state.get(request_hash, 0)
            cur.execute(
                "SELECT response_json FROM interactions WHERE request_hash=? AND sequence_num=?",
                (request_hash, next_seq),
            )
            row = cur.fetchone()
            if row is None:
                # If a state was provided and we missed, try last available for better error messages
                cur.execute(
                    "SELECT MAX(sequence_num) FROM interactions WHERE request_hash=?",
                    (request_hash,),
                )
                max_row = cur.fetchone()
                max_seq = max_row[0] if max_row else None
                available = []
                if max_seq is not None:
                    available = list(range(0, max_seq + 1))
                raise KeyError(
                    f"No recorded response for method={method} seq={next_seq}. Available sequences: {available}"
                )

            if state is not None:
                state[request_hash] = next_seq + 1

            return json.loads(row[0])
