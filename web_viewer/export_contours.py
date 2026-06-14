#!/usr/bin/env python3
import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path

import cv2
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


def contour_mode(name):
    if name == "none":
        return cv2.CHAIN_APPROX_NONE
    if name == "simple":
        return cv2.CHAIN_APPROX_SIMPLE
    raise ValueError(f"unsupported contour mode: {name}")


def tile_range(start, stop, tile_size):
    first = max(0, start // tile_size)
    last = max(0, stop // tile_size)
    return range(first, last + 1)


def cell_to_tiles(cell, tile_size):
    x0, y0, x1, y1 = cell["bbox"]
    for ty in tile_range(y0, y1, tile_size):
        for tx in tile_range(x0, x1, tile_size):
            yield tx, ty


def encode_contours(mask, cell_id, bbox, mode):
    x0, y0, x1, y1 = bbox
    crop = np.asarray(mask[y0:y1 + 1, x0:x1 + 1])
    binary = (crop == cell_id).astype(np.uint8)
    contours = cv2.findContours(binary, cv2.RETR_EXTERNAL, mode)[-2]

    encoded = []
    for contour in contours:
        points = contour.reshape(-1, 2)
        if len(points) < 2:
            continue
        encoded.append(points.astype(int).tolist())

    return encoded


def export_contours(mask_path, expr_csv, sample, out_dir, tile_size, chain, cell_ids=None, url_prefix="/data"):
    mask = np.load(mask_path, mmap_mode="r")
    if mask.ndim != 2:
        raise ValueError(f"{mask_path} must be a 2D label mask, got shape {mask.shape}")

    height, width = mask.shape
    expr_ids = set(range(1, int(mask.max()) + 1)) if cell_ids is None else set(cell_ids)
    objects = find_objects(mask)
    mode = contour_mode(chain)

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    tiles = defaultdict(list)
    exported = 0
    skipped_missing = 0
    skipped_empty = 0
    total_points = 0

    for index, cell_id in enumerate(sorted(expr_ids), 1):
        object_index = cell_id - 1
        if object_index < 0 or object_index >= len(objects) or objects[object_index] is None:
            skipped_missing += 1
            continue

        y_slice, x_slice = objects[object_index]
        bbox = [
            int(x_slice.start),
            int(y_slice.start),
            int(x_slice.stop - 1),
            int(y_slice.stop - 1),
        ]
        contours = encode_contours(mask, cell_id, bbox, mode)
        if not contours:
            skipped_empty += 1
            continue

        cell = {
            "id": int(cell_id),
            "bbox": bbox,
            "contours": contours,
        }
        for tile_key in cell_to_tiles(cell, tile_size):
            tiles[tile_key].append(cell)

        exported += 1
        total_points += sum(len(contour) for contour in contours)

        if index % 5000 == 0:
            print(
                f"{sample}: processed {index}/{len(expr_ids)} cells, "
                f"exported {exported}",
                flush=True,
            )

    tile_entries = []
    for tx, ty in sorted(tiles):
        payload = {
            "sample": sample,
            "tile": [tx, ty],
            "tileSize": int(tile_size),
            "cells": tiles[(tx, ty)],
        }
        filename = f"tile_{tx}_{ty}.json"
        (out_dir / filename).write_text(
            json.dumps(payload, separators=(",", ":")),
            encoding="utf-8",
        )

        tile_entries.append({
            "x": int(tx),
            "y": int(ty),
            "url": f"{url_prefix.rstrip('/')}/contours/{sample}/{filename}",
            "count": len(tiles[(tx, ty)]),
            "bounds": [
                int(tx * tile_size),
                int(ty * tile_size),
                int(min(width - 1, (tx + 1) * tile_size - 1)),
                int(min(height - 1, (ty + 1) * tile_size - 1)),
            ],
        })

    manifest = {
        "sample": sample,
        "width": int(width),
        "height": int(height),
        "tileSize": int(tile_size),
        "chain": chain,
        "cellCount": int(exported),
        "tileCount": len(tile_entries),
        "totalContourPoints": int(total_points),
        "skippedMissingMask": int(skipped_missing),
        "skippedEmptyContour": int(skipped_empty),
        "tiles": tile_entries,
    }
    (out_dir / "manifest.json").write_text(
        json.dumps(manifest, separators=(",", ":")),
        encoding="utf-8",
    )

    print(
        f"{sample}: wrote {len(tile_entries)} tiles, {exported} cells, "
        f"{total_points} contour points to {out_dir}",
        flush=True,
    )
    if skipped_missing or skipped_empty:
        print(
            f"{sample}: skipped {skipped_missing} cells missing from mask and "
            f"{skipped_empty} empty contours",
            flush=True,
        )


def main():
    parser = argparse.ArgumentParser(
        description="Export tiled cell contour JSON from a labelled cells.npy mask."
    )
    parser.add_argument("--mask", required=True)
    parser.add_argument("--expr")
    parser.add_argument("--cell-ids")
    parser.add_argument("--all-labels", action="store_true")
    parser.add_argument("--sample", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--url-prefix", default="/data")
    parser.add_argument("--tile-size", type=int, default=2048)
    parser.add_argument("--chain", choices=["simple", "none"], default="simple")
    args = parser.parse_args()

    if args.cell_ids:
        cell_ids = load_cell_ids(args.cell_ids)
    elif args.expr:
        cell_ids = load_expr_cell_ids(args.expr)
    elif args.all_labels:
        cell_ids = None
    else:
        raise ValueError("Provide --expr, --cell-ids, or --all-labels")

    export_contours(
        args.mask,
        args.expr,
        args.sample,
        args.out_dir,
        args.tile_size,
        args.chain,
        cell_ids,
        args.url_prefix,
    )


if __name__ == "__main__":
    main()
