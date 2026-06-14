#!/usr/bin/env python3
import argparse
import csv
import json
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path


def run_r_export(rda_path, object_name):
    with tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False) as handle:
        tsv_path = Path(handle.name)

    r_code = f"""
    suppressPackageStartupMessages(library(Seurat))
    load({json.dumps(str(Path(rda_path).resolve()))})
    obj <- get({json.dumps(object_name)})
    meta <- obj@meta.data
    cell_name <- rownames(meta)
    cell_id <- as.integer(sub(".*cell_", "", cell_name))
    out <- data.frame(rep=as.character(meta$orig.ident), cell_id=cell_id)
    write.table(out, {json.dumps(str(tsv_path))}, sep="\\t", quote=FALSE, row.names=FALSE)
    """

    subprocess.run(["Rscript", "-e", r_code], check=True)
    return tsv_path


def load_rep_cells(tsv_path):
    rep_cells = defaultdict(list)
    with open(tsv_path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            rep_cells[row["rep"]].append(int(row["cell_id"]))
    return rep_cells


def load_spatial_cells(path):
    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    return payload, {cell["id"]: cell for cell in payload["cells"]}


def parse_sample_rules(values):
    rules = []
    for value in values:
        if "=" not in value:
            raise ValueError(f"sample rule must be SAMPLE=prefix, got {value!r}")
        sample, prefix = value.split("=", 1)
        sample = sample.strip()
        prefix = prefix.strip().lower()
        if not sample or not prefix:
            raise ValueError(f"sample rule must be SAMPLE=prefix, got {value!r}")
        rules.append((sample, prefix))
    return rules


def sample_for_rep(rep, sample_rules):
    rep_lower = rep.lower()
    for sample, prefix in sample_rules:
        if rep_lower.startswith(prefix):
            return sample
    return None


def bbox_union(cells):
    if not cells:
        return None
    return [
        min(cell["bbox"][0] for cell in cells),
        min(cell["bbox"][1] for cell in cells),
        max(cell["bbox"][2] for cell in cells),
        max(cell["bbox"][3] for cell in cells),
    ]


def bbox_center(bbox):
    return ((bbox[0] + bbox[2]) / 2, (bbox[1] + bbox[3]) / 2)


def distance_to_bbox(point, bbox):
    x, y = point
    dx = max(bbox[0] - x, 0, x - bbox[2])
    dy = max(bbox[1] - y, 0, y - bbox[3])
    if dx == 0 and dy == 0:
        return 0
    return (dx * dx + dy * dy) ** 0.5


def tile_keys_for_bbox(bbox, tile_size, width, height):
    x0, y0, x1, y1 = bbox
    tx0 = max(0, x0 // tile_size)
    ty0 = max(0, y0 // tile_size)
    tx1 = min((width - 1) // tile_size, x1 // tile_size)
    ty1 = min((height - 1) // tile_size, y1 // tile_size)
    return [
        f"{tx},{ty}"
        for ty in range(ty0, ty1 + 1)
        for tx in range(tx0, tx1 + 1)
    ]


def natural_rep_label(rep):
    return rep.replace("_", " ").upper()


def assign_sample_replicates(sample, reps, cell_by_id, width, height, tile_size):
    cell_owners = defaultdict(list)
    for rep, ids in reps.items():
        for cell_id in ids:
            cell_owners[cell_id].append(rep)

    unique_bbox = {}
    for rep, ids in reps.items():
        cells = [
            cell_by_id[cell_id]
            for cell_id in ids
            if cell_id in cell_by_id and len(cell_owners[cell_id]) == 1
        ]
        unique_bbox[rep] = bbox_union(cells)

    assigned = {rep: [] for rep in reps}
    duplicate_assignments = []
    missing_ids = defaultdict(list)

    for rep, ids in reps.items():
        for cell_id in ids:
            if cell_id not in cell_by_id:
                missing_ids[rep].append(cell_id)
                continue
            if len(cell_owners[cell_id]) == 1:
                assigned[rep].append(cell_id)

    for cell_id, owners in sorted(cell_owners.items()):
        if len(owners) <= 1 or cell_id not in cell_by_id:
            continue

        point = bbox_center(cell_by_id[cell_id]["bbox"])
        ranked = []
        for rep in owners:
            bbox = unique_bbox.get(rep)
            if bbox is None:
                ranked.append((float("inf"), rep))
            else:
                ranked.append((distance_to_bbox(point, bbox), rep))
        ranked.sort()
        assigned_rep = ranked[0][1]
        assigned[assigned_rep].append(cell_id)
        duplicate_assignments.append({
            "cellId": int(cell_id),
            "reps": sorted(owners),
            "assignedRep": assigned_rep,
        })

    replicate_payloads = []
    for rep in sorted(reps):
        cell_ids = sorted(set(assigned[rep]))
        cells = [cell_by_id[cell_id] for cell_id in cell_ids if cell_id in cell_by_id]
        bbox = bbox_union(cells)
        if bbox is None:
            continue

        replicate_payloads.append({
            "id": rep,
            "label": natural_rep_label(rep),
            "sourceCellCount": len(reps[rep]),
            "assignedCellCount": len(cell_ids),
            "missingMaskCellCount": len(missing_ids[rep]),
            "bbox": bbox,
            "tileKeys": tile_keys_for_bbox(bbox, tile_size, width, height),
            "cellIds": cell_ids,
        })

    return replicate_payloads, duplicate_assignments


def export_replicates(rda, object_name, sample_inputs, contours_root, out_path, sample_rules):
    tsv_path = run_r_export(rda, object_name)
    try:
        rep_cells = load_rep_cells(tsv_path)
    finally:
        tsv_path.unlink(missing_ok=True)

    samples = {}
    duplicate_assignments = {}
    for sample, cells_path in sample_inputs:
        spatial, cell_by_id = load_spatial_cells(cells_path)
        manifest_path = Path(contours_root) / sample / "manifest.json"
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        reps = {
            rep: ids
            for rep, ids in rep_cells.items()
            if sample_for_rep(rep, sample_rules) == sample
        }
        replicate_payloads, duplicates = assign_sample_replicates(
            sample,
            reps,
            cell_by_id,
            int(spatial["width"]),
            int(spatial["height"]),
            int(manifest["tileSize"]),
        )
        samples[sample] = {
            "sample": sample,
            "width": int(spatial["width"]),
            "height": int(spatial["height"]),
            "tileSize": int(manifest["tileSize"]),
            "replicateCount": len(replicate_payloads),
            "replicates": replicate_payloads,
        }
        duplicate_assignments[sample] = duplicates

    payload = {
        "formatVersion": 1,
        "source": str(Path(rda).name),
        "samples": samples,
        "duplicateAssignments": duplicate_assignments,
    }

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(payload, separators=(",", ":")), encoding="utf-8")

    for sample, sample_payload in samples.items():
        print(
            f"{sample}: wrote {sample_payload['replicateCount']} replicates "
            f"to {out_path}",
            flush=True,
        )
        if duplicate_assignments[sample]:
            print(
                f"{sample}: assigned {len(duplicate_assignments[sample])} duplicate cell ids",
                flush=True,
            )


def main():
    parser = argparse.ArgumentParser(
        description="Export replicate metadata from a Seurat Rda object for the web viewer."
    )
    parser.add_argument("--rda", required=True)
    parser.add_argument("--object-name", default="st")
    parser.add_argument(
        "--sample",
        action="append",
        help="Sample cell metadata mapping SAMPLE=path/to/cells.json",
    )
    parser.add_argument("--cells-S1", default=None)
    parser.add_argument("--cells-S2", default=None)
    parser.add_argument("--contours-root", default="web_viewer/data/contours")
    parser.add_argument("--out", default="web_viewer/data/replicates.json")
    parser.add_argument(
        "--sample-rule",
        action="append",
        help="Map output sample id to orig.ident prefix, for example S1=s1_",
    )
    args = parser.parse_args()

    sample_inputs = []
    for value in args.sample or []:
        if "=" not in value:
            raise ValueError(f"--sample must be SAMPLE=path, got {value!r}")
        sample, path = value.split("=", 1)
        sample_inputs.append((sample.strip(), path.strip()))
    if not sample_inputs:
        if args.cells_S1:
            sample_inputs.append(("S1", args.cells_S1))
        else:
            sample_inputs.append(("S1", "web_viewer/data/S1_cells.json"))
        if args.cells_S2:
            sample_inputs.append(("S2", args.cells_S2))
        else:
            sample_inputs.append(("S2", "web_viewer/data/S2_cells.json"))

    sample_rules = parse_sample_rules(args.sample_rule) if args.sample_rule else [
        ("S1", "s1_"),
        ("S2", "s2_"),
    ]

    export_replicates(
        args.rda,
        args.object_name,
        sample_inputs,
        args.contours_root,
        args.out,
        sample_rules,
    )


if __name__ == "__main__":
    main()
