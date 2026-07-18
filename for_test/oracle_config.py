"""Configuration and identity helpers for the frozen Fortran oracle."""

import hashlib
import os
from pathlib import Path
from typing import Dict, Optional, Union


class OracleConfigurationError(RuntimeError):
    """Raised when the frozen Fortran oracle is not explicitly usable."""


PathLike = Union[str, Path]


def _configured_root(explicit_root: Optional[PathLike]) -> Path:
    configured = explicit_root or os.environ.get("OPENOCEANKRAKENC_FORTRAN_ROOT")
    if not configured:
        raise OracleConfigurationError(
            "Fortran oracle root is not configured; set "
            "OPENOCEANKRAKENC_FORTRAN_ROOT or pass an explicit root"
        )
    root = Path(configured).expanduser().resolve()
    if not root.is_dir():
        raise OracleConfigurationError(f"Fortran oracle root is not a directory: {root}")
    return root


def resolve_oracle_binary(name: str, explicit_root: Optional[PathLike] = None) -> Path:
    """Resolve *name* only below an explicitly configured oracle directory."""
    if Path(name).name != name:
        raise OracleConfigurationError(f"Oracle binary name must be a basename: {name}")
    root = _configured_root(explicit_root)
    binary = (root / name).resolve()
    if binary.parent != root or not binary.is_file():
        raise OracleConfigurationError(f"Missing Fortran oracle binary {name} under {root}")
    return binary


def oracle_identity(path: PathLike) -> Dict[str, str]:
    """Return the canonical path and SHA-256 of an oracle executable."""
    resolved = Path(path).resolve()
    if not resolved.is_file():
        raise OracleConfigurationError(f"Missing Fortran oracle binary: {resolved}")
    digest = hashlib.sha256()
    with resolved.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return {"path": str(resolved), "sha256": digest.hexdigest()}
