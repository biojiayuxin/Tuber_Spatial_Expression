#!/usr/bin/env python3
import argparse
import csv
import json
import sqlite3
import subprocess
import tempfile
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path


TISSUES_TABLE = "tissues"
ASSIGNMENTS_TABLE = "tissue_cell_assignments"


def sample_for_rep(rep):
    rep_lower = rep.lower()
    if rep_lower.startswith("s1_"):
        return "S1"
    if rep_lower.startswith("s2_"):
        return "S2"
    return None


def run_r_export(rda_path, object_name, tissue_column):
    with tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False) as handle:
        tsv_path = Path(handle.name)

    r_code = f"""
    suppressPackageStartupMessages(library(Seurat))
    load({json.dumps(str(Path(rda_path).resolve()))})
    obj <- get({json.dumps(object_name)})
    meta <- obj@meta.data
    tissue_column <- {json.dumps(tissue_column)}
    if (!(tissue_column %in% colnames(meta))) {{
      stop(sprintf("Tissue column '%s' was not found in obj@meta.data", tissue_column))
    }}
    if (!("orig.ident" %in% colnames(meta))) {{
      stop("Column 'orig.ident' was not found in obj@meta.data")
    }}
    cell_name <- rownames(meta)
    cell_id <- suppressWarnings(as.integer(sub(".*cell_", "", cell_name)))
    if (any(is.na(cell_id))) {{
      stop("Could not parse integer cell ids from metadata row names")
    }}
    tissue_id <- as.character(meta[[tissue_column]])
    if (any(is.na(tissue_id)) || any(!nzchar(tissue_id))) {{
      stop(sprintf("Tissue column '%s' contains missing values", tissue_column))
    }}
    out <- data.frame(
      rep = as.character(meta$orig.ident),
      cell_id = cell_id,
      tissue_id = tissue_id,
      check.names = FALSE
    )
    write.table(out, {json.dumps(str(tsv_path))}, sep="\\t", quote=FALSE, row.names=FALSE)
    """

    subprocess.run(["Rscript", "-e", r_code], check=True)
    return tsv_path


def load_tissue_rows(tsv_path):
    rows = []
    with open(tsv_path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            tissue_id = str(row["tissue_id"])
            if not tissue_id:
                raise ValueError("Empty tissue id in R export")
            rep = row["rep"]
            rows.append(
                {
                    "rep": rep,
                    "sample": sample_for_rep(rep),
                    "cell_id": int(row["cell_id"]),
                    "tissue_id": tissue_id,
                }
            )
    return rows


def load_replicates(path):
    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    return payload["samples"]


def choose_tissue(row_by_rep_cell, rows_by_sample_cell, rep_id, sample, cell_id):
    direct = row_by_rep_cell.get((rep_id, cell_id))
    if direct:
        return direct

    candidates = rows_by_sample_cell.get((sample, cell_id), [])
    tissue_ids = sorted({row["tissue_id"] for row in candidates})
    if len(tissue_ids) == 1:
        return tissue_ids[0]
    return None


def build_import_rows(rows, replicates):
    row_by_rep_cell = {
        (row["rep"], row["cell_id"]): row["tissue_id"]
        for row in rows
        if row["sample"]
    }
    rows_by_sample_cell = defaultdict(list)
    tissue_counts = Counter()
    tissue_order = []
    seen_tissues = set()

    for row in rows:
        tissue_id = row["tissue_id"]
        tissue_counts[tissue_id] += 1
        if tissue_id not in seen_tissues:
            seen_tissues.add(tissue_id)
            tissue_order.append(tissue_id)
        if row["sample"]:
            rows_by_sample_cell[(row["sample"], row["cell_id"])].append(row)

    assignment_rows = []
    assigned_counts = Counter()
    missing = defaultdict(list)

    for sample, sample_payload in replicates.items():
        seen_sample_cells = set()
        for rep in sample_payload.get("replicates", []):
            rep_id = rep["id"]
            for raw_cell_id in rep.get("cellIds", []):
                cell_id = int(raw_cell_id)
                if cell_id in seen_sample_cells:
                    continue
                tissue_id = choose_tissue(
                    row_by_rep_cell,
                    rows_by_sample_cell,
                    rep_id,
                    sample,
                    cell_id,
                )
                if tissue_id is None:
                    missing[sample].append({"rep": rep_id, "cellId": cell_id})
                    continue
                seen_sample_cells.add(cell_id)
                assignment_rows.append((sample, cell_id, tissue_id))
                assigned_counts[tissue_id] += 1

    tissue_rows = [
        (
            tissue_id,
            tissue_id,
            index,
            int(tissue_counts[tissue_id]),
            int(assigned_counts[tissue_id]),
        )
        for index, tissue_id in enumerate(tissue_order)
    ]
    return tissue_rows, assignment_rows, missing


def connect_db(path):
    conn = sqlite3.connect(path)
    conn.execute("PRAGMA temp_store=MEMORY")
    conn.execute("PRAGMA cache_size=-200000")
    return conn


def create_schema(conn):
    conn.executescript(
        f"""
        DROP TABLE IF EXISTS {ASSIGNMENTS_TABLE};
        DROP TABLE IF EXISTS {TISSUES_TABLE};

        CREATE TABLE {TISSUES_TABLE} (
          tissue_id TEXT PRIMARY KEY,
          tissue_label TEXT NOT NULL,
          tissue_order INTEGER NOT NULL,
          cell_count INTEGER NOT NULL,
          assigned_cell_count INTEGER NOT NULL
        );

        CREATE TABLE {ASSIGNMENTS_TABLE} (
          sample TEXT NOT NULL,
          cell_id INTEGER NOT NULL,
          tissue_id TEXT NOT NULL,
          PRIMARY KEY (sample, cell_id),
          FOREIGN KEY (tissue_id) REFERENCES {TISSUES_TABLE}(tissue_id)
        );

        CREATE INDEX idx_tissue_cell_assignments_tissue
          ON {ASSIGNMENTS_TABLE}(tissue_id, sample, cell_id);
        """
    )


def import_to_db(conn, tissue_rows, assignment_rows, args):
    create_schema(conn)
    conn.executemany(
        f"""
        INSERT INTO {TISSUES_TABLE}(
          tissue_id,
          tissue_label,
          tissue_order,
          cell_count,
          assigned_cell_count
        )
        VALUES (?, ?, ?, ?, ?)
        """,
        tissue_rows,
    )
    conn.executemany(
        f"""
        INSERT INTO {ASSIGNMENTS_TABLE}(sample, cell_id, tissue_id)
        VALUES (?, ?, ?)
        """,
        assignment_rows,
    )
    conn.executemany(
        "INSERT OR REPLACE INTO metadata(key, value) VALUES (?, ?)",
        [
            ("tissue_column", args.tissue_column),
            ("tissue_source", str(Path(args.rda).name)),
            ("tissue_imported_at", datetime.now(timezone.utc).isoformat()),
        ],
    )
    conn.execute("ANALYZE")
    conn.commit()


def import_tissues(args):
    if not Path(args.db).exists():
        raise FileNotFoundError(f"{args.db} not found; build expression.sqlite first")

    tsv_path = run_r_export(args.rda, args.object_name, args.tissue_column)
    try:
        rows = load_tissue_rows(tsv_path)
    finally:
        tsv_path.unlink(missing_ok=True)

    replicates = load_replicates(args.replicates_json)
    tissue_rows, assignment_rows, missing = build_import_rows(rows, replicates)

    with connect_db(args.db) as conn:
        import_to_db(conn, tissue_rows, assignment_rows, args)

    print(
        f"Imported {len(tissue_rows)} tissues and {len(assignment_rows):,} cell assignments "
        f"into {args.db}",
        flush=True,
    )
    for tissue_id, _, _, cell_count, assigned_count in tissue_rows:
        print(
            f"{tissue_id}: {assigned_count:,} assigned spatial cells "
            f"({cell_count:,} Seurat cells)",
            flush=True,
        )
    missing_total = sum(len(rows) for rows in missing.values())
    if missing_total:
        print(f"Warning: {missing_total:,} assigned spatial cells had no tissue", flush=True)


def main():
    parser = argparse.ArgumentParser(
        description="Import per-cell tissue assignments from a Seurat Rda object into SQLite."
    )
    parser.add_argument("--rda", default="seurat_object.celltype.Rda")
    parser.add_argument("--object-name", default="st")
    parser.add_argument("--tissue-column", default="celltype")
    parser.add_argument("--replicates-json", default="web_viewer/data/replicates.json")
    parser.add_argument("--db", default="web_viewer/data/expression.sqlite")
    args = parser.parse_args()

    import_tissues(args)


if __name__ == "__main__":
    main()
