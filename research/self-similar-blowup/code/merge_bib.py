"""Merge the per-section BibTeX files into bib/all.bib, de-duplicating keys.

Earlier files in ORDER take precedence (lit.bib is the verified master list).
Run from the project root: python3 code/merge_bib.py
"""
import re, pathlib

ORDER = ["lead", "lit", "mems", "power", "therm", "geom", "acou"]
root = pathlib.Path(__file__).resolve().parent.parent
seen, out, dropped = set(), [], []
for name in ORDER:
    path = root / "bib" / f"{name}.bib"
    if not path.exists():
        continue
    text = path.read_text()
    # split at entry starts (keep comments attached to the following entry)
    for m in re.finditer(r"@(\w+)\s*\{\s*([^,\s]+)\s*,.*?\n\}\s*$", text, re.S | re.M):
        key = m.group(2)
        if key in seen:
            dropped.append((name, key))
            continue
        seen.add(key)
        out.append(m.group(0).strip())
(root / "bib" / "all.bib").write_text("\n\n".join(out) + "\n")
print(f"{len(out)} entries written; {len(dropped)} duplicates dropped")
for name, key in dropped:
    print(f"  dropped {key} from {name}.bib")
