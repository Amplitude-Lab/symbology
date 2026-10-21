"""Pure table tests: no Mathematica or amplitude input required."""
import importlib.util
from pathlib import Path

spec = importlib.util.spec_from_file_location("ope_library", Path(__file__).resolve().parents[1] / "tools/ope_library.py")
lib = importlib.util.module_from_spec(spec)
spec.loader.exec_module(lib)
for kind, loop, bound, pair in [("MHV", 1, 4, 4), ("NMHV", 0, 4, 2)]:
    keys = list(lib.coordinates(kind, loop, bound, pair))
    rows = [list(key) + [0, 1] for key in keys]
    lib.validate_rows(rows, kind, loop, bound, pair)
    for invalid in [rows[:-1], rows + rows[:1], [rows[0][:-1] + [0]] + rows[1:]]:
        try:
            lib.validate_rows(invalid, kind, loop, bound, pair)
        except ValueError:
            pass
        else:
            raise AssertionError("Missing, duplicate or zero-denominator row accepted")
assert max(key[4] for key in lib.coordinates("MHV", 10, 2, 2)) == 18
assert max(key[4] for key in lib.coordinates("NMHV", 10, 2, 2)) == 20
print("PASS complete coordinate coverage, duplicates, denominators, L10 logs")
