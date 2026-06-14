#!/usr/bin/env python3
import argparse
import csv
import json
import subprocess
import tempfile
from collections import Counter, defaultdict
from pathlib import Path


def sample_for_rep(rep):
    rep_lower = rep.lower()
    if rep_lower.startswith("s1_"):
        return "S1"
    if rep_lower.startswith("s2_"):
        return "S2"
    return None


def cluster_sort_key(cluster_id):
    try:
        return (0, int(cluster_id))
    except ValueError:
        return (1, cluster_id)


def run_r_export(rda_path, object_name, cluster_column):
    with tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False) as handle:
        tsv_path = Path(handle.name)

    r_code = f"""
    suppressPackageStartupMessages(library(Seurat))
    load({json.dumps(str(Path(rda_path).resolve()))})
    obj <- get({json.dumps(object_name)})
    meta <- obj@meta.data
    cluster_column <- {json.dumps(cluster_column)}
    if (!(cluster_column %in% colnames(meta))) {{
      stop(sprintf("Cluster column '%s' was not found in obj@meta.data", cluster_column))
    }}
    if (!("orig.ident" %in% colnames(meta))) {{
      stop("Column 'orig.ident' was not found in obj@meta.data")
    }}
    cell_name <- rownames(meta)
    cell_id <- suppressWarnings(as.integer(sub(".*cell_", "", cell_name)))
    if (any(is.na(cell_id))) {{
      stop("Could not parse integer cell ids from metadata row names")
    }}
    out <- data.frame(
      rep = as.character(meta$orig.ident),
      cell_id = cell_id,
      cluster_id = as.character(meta[[cluster_column]]),
      check.names = FALSE
    )
    write.table(out, {json.dumps(str(tsv_path))}, sep="\\t", quote=FALSE, row.names=FALSE)
    """

    subprocess.run(["Rscript", "-e", r_code], check=True)
    return tsv_path


def load_cluster_rows(tsv_path):
    rows = []
    with open(tsv_path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            cluster_id = str(row["cluster_id"])
            if not cluster_id:
                raise ValueError("Empty cluster id in R export")
            rows.append(
                {
                    "rep": row["rep"],
                    "sample": sample_for_rep(row["rep"]),
                    "cell_id": int(row["cell_id"]),
                    "cluster_id": cluster_id,
                }
            )
    return rows


def load_replicates(path):
    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    return payload["samples"]


def choose_cluster(row_by_rep_cell, rows_by_sample_cell, rep_id, sample, cell_id):
    direct = row_by_rep_cell.get((rep_id, cell_id))
    if direct:
        return direct

    candidates = rows_by_sample_cell.get((sample, cell_id), [])
    cluster_ids = sorted({row["cluster_id"] for row in candidates}, key=cluster_sort_key)
    if len(cluster_ids) == 1:
        return cluster_ids[0]
    return None


def build_payload(rows, replicates, rda, cluster_column):
    row_by_rep_cell = {
        (row["rep"], row["cell_id"]): row["cluster_id"]
        for row in rows
        if row["sample"]
    }
    rows_by_sample_cell = defaultdict(list)
    cluster_counts = Counter()
    for row in rows:
        if not row["sample"]:
            continue
        rows_by_sample_cell[(row["sample"], row["cell_id"])].append(row)
        cluster_counts[row["cluster_id"]] += 1

    cluster_ids = sorted(cluster_counts, key=cluster_sort_key)
    assigned_counts = Counter()
    samples = {}
    missing = defaultdict(list)

    for sample, sample_payload in replicates.items():
        cells = []
        for rep in sample_payload.get("replicates", []):
            rep_id = rep["id"]
            for cell_id in rep.get("cellIds", []):
                cluster_id = choose_cluster(
                    row_by_rep_cell,
                    rows_by_sample_cell,
                    rep_id,
                    sample,
                    int(cell_id),
                )
                if cluster_id is None:
                    missing[sample].append({"rep": rep_id, "cellId": int(cell_id)})
                    continue
                cells.append([int(cell_id), cluster_id])
                assigned_counts[cluster_id] += 1

        samples[sample] = {
            "sample": sample,
            "assignedCellCount": len(cells),
            "cells": cells,
        }

    return {
        "formatVersion": 1,
        "source": str(Path(rda).name),
        "clusterColumn": cluster_column,
        "clusters": [
            {
                "id": cluster_id,
                "label": cluster_id,
                "order": index,
                "cellCount": int(cluster_counts[cluster_id]),
                "assignedCellCount": int(assigned_counts[cluster_id]),
            }
            for index, cluster_id in enumerate(cluster_ids)
        ],
        "samples": samples,
        "missingAssignments": missing,
    }


def export_clusters(args):
    tsv_path = run_r_export(args.rda, args.object_name, args.cluster_column)
    try:
        rows = load_cluster_rows(tsv_path)
    finally:
        tsv_path.unlink(missing_ok=True)

    replicates = load_replicates(args.replicates_json)
    payload = build_payload(rows, replicates, args.rda, args.cluster_column)

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(
        json.dumps(payload, ensure_ascii=False, separators=(",", ":")),
        encoding="utf-8",
    )

    print(
        f"Wrote {len(payload['clusters'])} clusters to {out_path}",
        flush=True,
    )
    for sample, sample_payload in payload["samples"].items():
        print(
            f"{sample}: wrote {sample_payload['assignedCellCount']:,} cell cluster assignments",
            flush=True,
        )
    missing_total = sum(len(rows) for rows in payload["missingAssignments"].values())
    if missing_total:
        print(f"Warning: {missing_total:,} assigned spatial cells had no cluster", flush=True)


def main():
    parser = argparse.ArgumentParser(
        description="Export per-cell Seurat cluster assignments for the web viewer."
    )
    parser.add_argument("--rda", default="seurat_object.02-dims64.res0.6.Rda")
    parser.add_argument("--object-name", default="st")
    parser.add_argument("--cluster-column", default="seurat_clusters")
    parser.add_argument("--replicates-json", default="web_viewer/data/replicates.json")
    parser.add_argument("--out", default="web_viewer/data/clusters.json")
    args = parser.parse_args()

    export_clusters(args)


if __name__ == "__main__":
    main()
