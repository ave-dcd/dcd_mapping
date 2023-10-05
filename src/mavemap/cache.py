"""Provide utilities for caching requests."""
import os
from pathlib import Path

LOCAL_STORE_PATH = Path(
    os.environ.get("MAVEDB_STORAGE_DIR", Path.home() / ".local" / "share" / "mavemap")
)
if not LOCAL_STORE_PATH.exists():
    LOCAL_STORE_PATH.mkdir(exist_ok=True, parents=True)

# TODO add more robust cache management methods, esp. for API calls
