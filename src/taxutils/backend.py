"""Compatibility checks and error normalization for the required Rust backend."""

from __future__ import annotations

from typing import Any


BACKEND_API_VERSION = 2


class TaxutilsBackendError(RuntimeError):
    """Raised when the native backend cannot complete an operation."""


try:
    from . import _rust
except ImportError as exc:
    raise ImportError(
        "taxutils requires its native Rust extension; reinstall a supported "
        "wheel or build from source with Rust available"
    ) from exc

try:
    _API_VERSION = int(_rust.api_version())
except Exception as exc:
    raise ImportError("The taxutils Rust extension did not report its API version") from exc

if _API_VERSION != BACKEND_API_VERSION:
    raise ImportError(
        f"Rust backend API {_API_VERSION} is incompatible with required API "
        f"{BACKEND_API_VERSION}"
    )


def call_rust(name: str, *args: Any, **kwargs: Any):
    """Call one native operation and normalize its backend-specific exception."""
    try:
        return getattr(_rust, name)(*args, **kwargs)
    except _rust.TaxutilsBackendError as exc:
        raise TaxutilsBackendError(str(exc)) from exc


def backend_info() -> dict[str, object]:
    """Describe the required native backend."""
    return {
        "selected": "rust",
        "rust_version": str(_rust.crate_version()),
        "api_version": _API_VERSION,
    }
