#!/usr/bin/env python3
"""Verify the public snapshot with the Python standard library."""
import hashlib, json
from pathlib import Path

def main():
    root = Path(__file__).resolve().parent
    manifest = json.loads((root / "MANIFEST-PUBLIC.json").read_text())
    errors = []
    for name, record in manifest.items():
        rel = Path(name)
        if rel.is_absolute() or ".." in rel.parts:
            errors.append(f"Unsafe path: {name}")
            continue
        p = root / rel
        if not p.is_file():
            errors.append(f"Missing: {name}")
            continue
        with p.open("rb") as stream:
            actual = hashlib.file_digest(stream, "sha256").hexdigest()
        if p.stat().st_size != record["bytes"] or actual != record["sha256"]:
            errors.append(f"Changed: {name}")
    ignored_parts = {'.git', '.venv', '__pycache__'}
    actual_names = {p.relative_to(root).as_posix() for p in root.rglob('*')
                    if p.is_file() and not ignored_parts.intersection(p.relative_to(root).parts)
                    and p.suffix != '.pyc' and p.name != 'MANIFEST-PUBLIC.json'}
    unexpected = actual_names - set(manifest)
    errors.extend('Unexpected: ' + name for name in sorted(unexpected))
    if errors:
        raise SystemExit("\n".join(errors))
    print(f"PASS: {len(manifest)} files match the public snapshot.")

if __name__ == "__main__":
    main()
