"""Locate named test cases in either flat or categorized fixture layouts."""

from pathlib import Path


def test_directory() -> Path:
    for ancestor in Path(__file__).resolve().parents:
        candidate = ancestor / "test"
        if candidate.is_dir() and any(candidate.rglob("*.env")):
            return candidate
    raise FileNotFoundError("unable to locate workspace test directory")


def case_root(name: str) -> Path:
    matches = list(test_directory().rglob(f"{name}.env"))
    if len(matches) != 1:
        raise FileNotFoundError(
            f"expected exactly one ENV named {name}.env, found {len(matches)}"
        )
    return matches[0].with_suffix("")
