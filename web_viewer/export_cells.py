#!/usr/bin/env python3
import argparse
import csv
import json
from pathlib import Path

import numpy as np
from scipy.ndimage import find_objects


def load_expr_cell_ids(expr_csv):
    ids = set()
    with open(expr_csv, newline="") as handle:
        reader = csv.reader(handle)
        next(reader)
        for row in reader:
            if row:
                ids.add(int(row[0]))
    return ids


def load_cell_ids(path):
    ids = set()
    with open(path, newline="") as handle:
        reader = csv.reader(handle)
        header = next(reader)
        if header and header[0] in {"cell", "cell_id", "id"}:
            for row in reader:
                if row:
                    ids.add(int(row[0]))
        else:
            if header:
                ids.add(int(header[0]))
            for row in reader:
                if row:
                    ids.add(int(row[0]))
    return ids


def summarize_mask(mask_path, expr_csv, sample, output_json, chunk_rows=512, cell_ids=None):
    mask = np.load(mask_path, mmap_mode="r")
    if mask.ndim != 2:
        raise ValueError(f"{mask_path} must be a 2D label mask, got shape {mask.shape}")

    max_id = int(mask.max())
    height, width = mask.shape
    expr_ids = set(range(1, max_id + 1)) if cell_ids is None else set(cell_ids)

    count = np.zeros(max_id + 1, dtype=np.int64)
    sum_x = np.zeros(max_id + 1, dtype=np.float64)
    sum_y = np.zeros(max_id + 1, dtype=np.float64)
    x_coords = np.arange(width, dtype=np.float64)

    for y0 in range(0, height, chunk_rows):
        chunk = np.asarray(mask[y0:y0 + chunk_rows])
        labels = chunk.reshape(-1)
        nz = labels > 0
        if not np.any(nz):
            continue

        labels = labels[nz]
        local_index = np.nonzero(nz)[0]
        ys = (local_index // width + y0).astype(np.float64)
        xs = np.take(x_coords, local_index % width)

        count += np.bincount(labels, minlength=max_id + 1)
        sum_x += np.bincount(labels, weights=xs, minlength=max_id + 1)
        sum_y += np.bincount(labels, weights=ys, minlength=max_id + 1)

        print(f"{sample}: scanned rows {min(y0 + chunk_rows, height)}/{height}", flush=True)

    print(f"{sample}: computing bounding boxes", flush=True)
    object_slices = find_objects(mask)

    cells = []
    for cell_id in sorted(expr_ids):
        if cell_id >= len(count) or count[cell_id] == 0:
            continue
        slices = object_slices[cell_id - 1] if cell_id - 1 < len(object_slices) else None
        if slices is None:
            continue
        y_slice, x_slice = slices
        area = int(count[cell_id])
        cells.append({
            "id": int(cell_id),
            "x": round(float(sum_x[cell_id] / area), 3),
            "y": round(float(sum_y[cell_id] / area), 3),
            "area": area,
            "bbox": [
                int(x_slice.start),
                int(y_slice.start),
                int(x_slice.stop - 1),
                int(y_slice.stop - 1),
            ],
        })

    payload = {
        "sample": sample,
        "width": int(width),
        "height": int(height),
        "cellCount": len(cells),
        "cells": cells,
    }

    output_json = Path(output_json)
    output_json.parent.mkdir(parents=True, exist_ok=True)
    output_json.write_text(json.dumps(payload, separators=(",", ":")), encoding="utf-8")
    print(f"Wrote {output_json} with {len(cells)} cells")


def main():
    parser = argparse.ArgumentParser(
        description="Export compact per-cell spatial metadata from a labelled cells.npy mask."
    )
    parser.add_argument("--mask", required=True)
    parser.add_argument("--expr")
    parser.add_argument("--cell-ids")
    parser.add_argument("--all-labels", action="store_true")
    parser.add_argument("--sample", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--chunk-rows", type=int, default=512)
    args = parser.parse_args()

    if args.cell_ids:
        cell_ids = load_cell_ids(args.cell_ids)
    elif args.expr:
        cell_ids = load_expr_cell_ids(args.expr)
    elif args.all_labels:
        cell_ids = None
    else:
        raise ValueError("Provide --expr, --cell-ids, or --all-labels")

    summarize_mask(args.mask, args.expr, args.sample, args.out, args.chunk_rows, cell_ids)


if __name__ == "__main__":
    main()
