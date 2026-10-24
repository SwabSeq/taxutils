"""Compatibility checks and error normalization for the required Rust backend."""

from __future__ import annotations

from typing import Any


class TaxutilsBackendError(RuntimeError):
    """Raised when the native backend cannot complete an operation."""


try:
    from . import _rust
except ImportError as exc:
    raise ImportError(
        "taxutils requires its native Rust extension; reinstall a supported "
        "wheel or build from source with Rust available"
    ) from exc


def _check_backend():
    try:
        api_version = int(_rust.api_version())
    except Exception as exc:
        raise ImportError("The taxutils Rust extension did not report its API version") from exc
    if api_version != 8:
        raise ImportError(
            f"Rust backend API {api_version} is incompatible with required API 8"
        )


_check_backend()


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
        "api_version": int(_rust.api_version()),
    }
