#!/usr/bin/env python3
"""Print the sheet layout of a downloaded workbook (untrusted data; values only).
Usage: python3 -I peek_xlsx.py FILE.xlsx [max_rows] [max_cols]"""
import sys
import openpyxl
path = sys.argv[1]
max_rows = int(sys.argv[2]) if len(sys.argv) > 2 else 12
max_cols = int(sys.argv[3]) if len(sys.argv) > 3 else 40
wb = openpyxl.load_workbook(path, data_only=True)
for ws in wb.worksheets:
    print(f"=== sheet {ws.title!r} dims={ws.dimensions} max_row={ws.max_row} max_col={ws.max_column}")
    for i, row in enumerate(ws.iter_rows(values_only=True)):
        if i >= max_rows:
            break
        cells = ["" if c is None else (f"{c:.4g}" if isinstance(c, float) else str(c)[:28]) for c in row[:max_cols]]
        print(i + 1, "|", " | ".join(cells))
