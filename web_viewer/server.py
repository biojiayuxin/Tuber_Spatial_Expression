#!/usr/bin/env python3
import argparse
import csv
import json
import math
import sqlite3
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import parse_qs, urlparse

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
WEB_ROOT = Path(__file__).resolve().parent
STATIC_ROOT = WEB_ROOT / "static"
CACHE_ROOT = WEB_ROOT / "cache"
DATASETS_CONFIG_PATH = WEB_ROOT / "datasets.json"
CATEGORY_COLORS_PATH = ROOT / "colors.txt"
CACHE_VERSION = 2
DOTPLOT_CLUSTER_COLUMN = "seurat_clusters"
TISSUE_COLUMN = "celltype"
TISSUES_TABLE = "tissues"
TISSUE_ASSIGNMENTS_TABLE = "tissue_cell_assignments"


def load_dataset_catalog():
    payload = json.loads(DATASETS_CONFIG_PATH.read_text(encoding="utf-8"))
    datasets = []
    for dataset in payload.get("datasets", []):
        item = dict(dataset)
        item["data_root"] = (WEB_ROOT / item["dataRoot"]).resolve()
        item["db_path"] = item["data_root"] / "expression.sqlite"
        item["cache_root"] = CACHE_ROOT / item["id"]
        samples = []
        for sample in item.get("samples", []):
            sample_item = dict(sample)
            if sample_item.get("csv"):
                sample_item["csv_path"] = (ROOT / sample_item["csv"]).resolve()
            samples.append(sample_item)
        item["samples"] = samples
        item["sample_ids"] = [sample["id"] for sample in samples]
        datasets.append(item)

    dataset_by_id = {dataset["id"]: dataset for dataset in datasets}
    default_id = payload.get("defaultDataset") or (datasets[0]["id"] if datasets else "")
    return {
        "defaultDataset": default_id,
        "datasets": datasets,
        "datasetById": dataset_by_id,
    }


def get_dataset(dataset_id=None):
    catalog = load_dataset_catalog()
    selected_id = dataset_id or catalog["defaultDataset"]
    dataset = catalog["datasetById"].get(selected_id)
    if dataset is None:
        raise KeyError(f"dataset {selected_id!r} not found")
    return dataset


def public_dataset_payload():
    catalog = load_dataset_catalog()
    return {
        "defaultDataset": catalog["defaultDataset"],
        "datasets": [
            {
                "id": dataset["id"],
                "label": dataset.get("label", dataset["id"]),
                "dataPath": dataset.get("dataPath", f"/dataset-data/{dataset['id']}"),
                "defaultSample": dataset.get("defaultSample") or dataset["sample_ids"][0],
                "defaultGene": dataset.get("defaultGene", ""),
                "samples": [
                    {
                        "id": sample["id"],
                        "label": sample.get("label", sample["id"]),
                        "columns": int(sample.get("columns") or 0),
                        "contoursPath": sample.get("contoursPath", ""),
                    }
                    for sample in dataset["samples"]
                ],
            }
            for dataset in catalog["datasets"]
        ],
    }


def sample_config(dataset, sample_id):
    for sample in dataset["samples"]:
        if sample["id"] == sample_id:
            return sample
    raise KeyError(f"sample {sample_id!r} not found in dataset {dataset['id']!r}")


def connect_db(dataset):
    return sqlite3.connect(f"file:{dataset['db_path']}?mode=ro", uri=True)


def json_response(handler, status, payload):
    body = json.dumps(payload, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
    handler.send_response(status)
    handler.send_header("Content-Type", "application/json; charset=utf-8")
    handler.send_header("Content-Length", str(len(body)))
    handler.send_header("Cache-Control", "no-store")
    handler.end_headers()
    handler.wfile.write(body)


def normalize_hex_color(value):
    color = value.strip()
    if len(color) == 7 and color[0] == "#":
        digits = color[1:]
    elif len(color) == 6:
        digits = color
        color = f"#{color}"
    else:
        return None

    if all(char in "0123456789abcdefABCDEF" for char in digits):
        return color.upper()
    return None


def load_category_colors():
    payload = {
        "formatVersion": 1,
        "clusters": {},
        "tissues": {},
    }
    if not CATEGORY_COLORS_PATH.exists():
        return payload

    section = None
    for raw_line in CATEGORY_COLORS_PATH.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if not line:
            continue

        if line.startswith("#"):
            heading = line.lstrip("#").strip().lower()
            if "cluster" in heading:
                section = "clusters"
            elif "tissue" in heading:
                section = "tissues"
            continue

        if section not in payload:
            continue

        try:
            label, color_value = line.rsplit(None, 1)
        except ValueError:
            continue

        color = normalize_hex_color(color_value)
        if color is not None:
            payload[section][label.strip()] = color

    return payload


def read_gene_names(csv_path):
    with open(csv_path, newline="") as handle:
        reader = csv.reader(handle)
        header = next(reader)
    return header[1:]


def get_gene_names(dataset):
    if dataset["db_path"].exists():
        try:
            with connect_db(dataset) as conn:
                rows = conn.execute("SELECT gene FROM genes ORDER BY gene_id").fetchall()
            if rows:
                return [row[0] for row in rows]
        except sqlite3.Error as exc:
            print(f"SQLite gene list unavailable, falling back to JSON/CSV: {exc}", flush=True)

    genes_path = dataset["data_root"] / "genes.json"
    if genes_path.exists():
        return json.loads(genes_path.read_text(encoding="utf-8"))

    csv_samples = [sample for sample in dataset["samples"] if sample.get("csv_path")]
    if not csv_samples:
        raise RuntimeError(
            f"gene list unavailable for dataset {dataset['id']}; build expression.sqlite first"
        )
    genes = read_gene_names(csv_samples[0]["csv_path"])
    dataset["data_root"].mkdir(parents=True, exist_ok=True)
    genes_path.write_text(json.dumps(genes, separators=(",", ":")), encoding="utf-8")
    return genes


def cache_path(dataset, sample, gene):
    safe_gene = gene.replace("/", "_").replace("\\", "_")
    return dataset["cache_root"] / f"{sample}_{safe_gene}.json"


def load_expression_from_db(dataset, sample, gene):
    with connect_db(dataset) as conn:
        gene_row = conn.execute(
            "SELECT gene_id FROM genes WHERE gene = ?",
            (gene,),
        ).fetchone()
        if gene_row is None:
            raise ValueError(f"{gene} not found")

        gene_id = int(gene_row[0])
        stats = conn.execute(
            """
            SELECT vmin, vmax, nonzero
            FROM sample_genes
            WHERE sample = ? AND gene_id = ?
            """,
            (sample, gene_id),
        ).fetchone()
        if stats is None:
            raise ValueError(f"{gene} not found in {sample}")

        rows = conn.execute(
            """
            SELECT cell_id, value
            FROM expression_values
            WHERE sample = ? AND gene_id = ?
            ORDER BY cell_id
            """,
            (sample, gene_id),
        ).fetchall()

    return {
        "formatVersion": CACHE_VERSION,
        "source": "sqlite",
        "sample": sample,
        "gene": gene,
        "vmin": float(stats[0]),
        "vmax": float(stats[1]),
        "nonzero": int(stats[2]),
        "values": [[int(cell_id), float(value)] for cell_id, value in rows],
    }


def load_expression_from_csv(dataset, sample, gene):
    path = cache_path(dataset, sample, gene)
    if path.exists():
        payload = json.loads(path.read_text(encoding="utf-8"))
        if (
            payload.get("formatVersion") == CACHE_VERSION
            and "vmin" in payload
            and "vmax" in payload
            and "values" in payload
        ):
            return payload

    csv_path = sample_config(dataset, sample).get("csv_path")
    if csv_path is None:
        raise RuntimeError(
            f"CSV fallback unavailable for sample {sample} in dataset {dataset['id']}"
        )
    df = pd.read_csv(csv_path, usecols=["cell", gene])
    values_all = df[gene].astype(float)
    df_nonzero = df[values_all > 0]

    values = [
        [int(cell_id), float(value)]
        for cell_id, value in zip(df_nonzero["cell"].to_numpy(), df_nonzero[gene].to_numpy())
    ]

    payload = {
        "formatVersion": CACHE_VERSION,
        "sample": sample,
        "gene": gene,
        "vmin": float(values_all.min()) if len(values_all) else 0.0,
        "vmax": float(values_all.max()) if len(values_all) else 0.0,
        "nonzero": len(values),
        "values": values,
    }

    dataset["cache_root"].mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, separators=(",", ":")), encoding="utf-8")
    return payload


def load_expression(dataset, sample, gene):
    if dataset["db_path"].exists():
        try:
            return load_expression_from_db(dataset, sample, gene)
        except sqlite3.Error as exc:
            print(f"SQLite expression lookup failed, falling back to CSV/cache: {exc}", flush=True)

    return load_expression_from_csv(dataset, sample, gene)


def db_table_exists(conn, name):
    return conn.execute(
        "SELECT 1 FROM sqlite_master WHERE type = 'table' AND name = ?",
        (name,),
    ).fetchone() is not None


def get_metadata_value(conn, key, default=None):
    try:
        row = conn.execute(
            "SELECT value FROM metadata WHERE key = ?",
            (key,),
        ).fetchone()
    except sqlite3.Error:
        return default
    return row[0] if row is not None else default


def load_dotplot_from_db(dataset, gene):
    if not dataset["db_path"].exists():
        raise RuntimeError(
            f"dotplot database not found for dataset {dataset['id']}; build expression.sqlite first"
        )

    with connect_db(dataset) as conn:
        if not (
            db_table_exists(conn, "dotplot_clusters")
            and db_table_exists(conn, "dotplot_gene_cluster_stats")
        ):
            raise RuntimeError(
                "dotplot tables not found; run web_viewer/import_dotplot_stats.py"
            )

        gene_row = conn.execute(
            "SELECT gene_id FROM genes WHERE gene = ?",
            (gene,),
        ).fetchone()
        if gene_row is None:
            raise ValueError(f"{gene} not found")

        gene_id = int(gene_row[0])
        cluster_column = get_metadata_value(
            conn,
            "dotplot_cluster_column",
            DOTPLOT_CLUSTER_COLUMN,
        )
        rows = conn.execute(
            """
            SELECT
              c.cluster_id,
              c.cluster_order,
              c.cell_count,
              s.avg_expr,
              s.pct_expr,
              s.expressing_count
            FROM dotplot_clusters c
            JOIN dotplot_gene_cluster_stats s
              ON s.cluster_id = c.cluster_id
             AND s.gene_id = ?
            ORDER BY c.cluster_order
            """,
            (gene_id,),
        ).fetchall()

    if not rows:
        raise ValueError(f"dotplot stats for {gene} not found")

    avg_values = [float(row[3]) for row in rows]
    logged_avg = [math.log1p(value) for value in avg_values]
    if len(logged_avg) > 1:
        mean_logged = sum(logged_avg) / len(logged_avg)
        variance = sum((value - mean_logged) ** 2 for value in logged_avg) / (len(logged_avg) - 1)
        sd_logged = math.sqrt(variance)
    else:
        sd_logged = 0.0
    scaled_avg = [
        max(-2.5, min(2.5, (value - mean_logged) / sd_logged)) if sd_logged > 0 else 0.0
        for value in logged_avg
    ]

    return {
        "gene": gene,
        "clusterColumn": cluster_column,
        "clusters": [
            {
                "id": str(cluster_id),
                "label": str(cluster_id),
                "order": int(cluster_order),
                "cellCount": int(cell_count),
                "avgExpr": float(avg_expr),
                "avgExprScaled": float(scaled_avg[index]),
                "pctExpr": float(pct_expr),
                "expressingCount": int(expressing_count),
            }
            for index, (
                cluster_id,
                cluster_order,
                cell_count,
                avg_expr,
                pct_expr,
                expressing_count,
            ) in enumerate(rows)
        ],
    }


def load_tissues_from_db(dataset):
    if not dataset["db_path"].exists():
        raise RuntimeError(
            f"tissue database not found for dataset {dataset['id']}; build expression.sqlite first"
        )

    with connect_db(dataset) as conn:
        if not (
            db_table_exists(conn, TISSUES_TABLE)
            and db_table_exists(conn, TISSUE_ASSIGNMENTS_TABLE)
        ):
            raise RuntimeError("tissue tables not found; run web_viewer/import_tissues.py")

        tissue_column = get_metadata_value(conn, "tissue_column", TISSUE_COLUMN)
        source = get_metadata_value(conn, "tissue_source", "")
        tissue_rows = conn.execute(
            f"""
            SELECT
              tissue_id,
              tissue_label,
              tissue_order,
              cell_count,
              assigned_cell_count
            FROM {TISSUES_TABLE}
            ORDER BY tissue_order, tissue_label
            """
        ).fetchall()
        assignment_rows = conn.execute(
            f"""
            SELECT sample, cell_id, tissue_id
            FROM {TISSUE_ASSIGNMENTS_TABLE}
            ORDER BY sample, cell_id
            """
        ).fetchall()

    samples = {}
    for sample in dataset["sample_ids"]:
        samples[sample] = {
            "sample": sample,
            "assignedCellCount": 0,
            "cells": [],
        }

    for sample, cell_id, tissue_id in assignment_rows:
        sample_payload = samples.setdefault(
            sample,
            {
                "sample": sample,
                "assignedCellCount": 0,
                "cells": [],
            },
        )
        sample_payload["cells"].append([int(cell_id), str(tissue_id)])
        sample_payload["assignedCellCount"] += 1

    return {
        "formatVersion": 1,
        "source": source,
        "tissueColumn": tissue_column,
        "tissues": [
            {
                "id": str(tissue_id),
                "label": str(tissue_label),
                "order": int(tissue_order),
                "cellCount": int(cell_count),
                "assignedCellCount": int(assigned_cell_count),
            }
            for tissue_id, tissue_label, tissue_order, cell_count, assigned_cell_count
            in tissue_rows
        ],
        "samples": samples,
    }


def get_json_document(dataset, name, fallback_path):
    if dataset["db_path"].exists():
        try:
            with connect_db(dataset) as conn:
                row = conn.execute(
                    "SELECT payload FROM json_documents WHERE name = ?",
                    (name,),
                ).fetchone()
            if row is not None:
                return json.loads(row[0])
        except sqlite3.Error as exc:
            print(f"SQLite JSON document unavailable, falling back to file: {exc}", flush=True)

    path = Path(fallback_path)
    if path.exists():
        return json.loads(path.read_text(encoding="utf-8"))
    return None


class Handler(SimpleHTTPRequestHandler):
    def translate_path(self, path):
        parsed = urlparse(path)
        clean_path = parsed.path.lstrip("/")

        if clean_path.startswith("data/"):
            return str(WEB_ROOT / clean_path)

        if clean_path.startswith("dataset-data/"):
            parts = clean_path.split("/", 2)
            if len(parts) >= 3:
                dataset_id = parts[1]
                try:
                    dataset = get_dataset(dataset_id)
                except KeyError:
                    return str(WEB_ROOT / "__missing_dataset__")
                target = (dataset["data_root"] / parts[2]).resolve()
                if target == dataset["data_root"] or dataset["data_root"] in target.parents:
                    return str(target)
                return str(WEB_ROOT / "__invalid_dataset_path__")
            return str(WEB_ROOT / "__missing_dataset_file__")

        if clean_path == "":
            return str(STATIC_ROOT / "index.html")

        return str(STATIC_ROOT / clean_path)

    def do_GET(self):
        parsed = urlparse(self.path)
        params = parse_qs(parsed.query)

        if parsed.path == "/api/datasets":
            json_response(self, 200, public_dataset_payload())
            return

        try:
            dataset = get_dataset((params.get("dataset") or [""])[0].strip() or None)
        except KeyError as exc:
            json_response(self, 404, {"error": str(exc)})
            return

        if parsed.path == "/api/genes":
            try:
                genes = get_gene_names(dataset)
            except RuntimeError as exc:
                json_response(self, 503, {"error": str(exc)})
                return
            json_response(self, 200, {"dataset": dataset["id"], "genes": genes})
            return

        if parsed.path == "/api/colors":
            json_response(self, 200, load_category_colors())
            return

        if parsed.path == "/api/replicates":
            replicates = get_json_document(dataset, "replicates", dataset["data_root"] / "replicates.json")
            if replicates is None:
                json_response(self, 404, {"error": "replicates not found"})
                return
            json_response(self, 200, replicates)
            return

        if parsed.path == "/api/tissues":
            try:
                payload = load_tissues_from_db(dataset)
            except RuntimeError as exc:
                json_response(self, 503, {"error": str(exc)})
                return
            except sqlite3.Error as exc:
                json_response(self, 500, {"error": f"tissue SQLite error: {exc}"})
                return

            json_response(self, 200, payload)
            return

        if parsed.path == "/api/gene":
            gene = (params.get("gene") or [""])[0].strip()
            if not gene:
                json_response(self, 400, {"error": "missing gene"})
                return

            try:
                genes = get_gene_names(dataset)
            except RuntimeError as exc:
                json_response(self, 503, {"error": str(exc)})
                return
            if gene not in genes:
                json_response(self, 404, {"error": f"{gene} not found"})
                return

            try:
                samples = {
                    sample: load_expression(dataset, sample, gene)
                    for sample in dataset["sample_ids"]
                }
                payload = {
                    "dataset": dataset["id"],
                    "gene": gene,
                    "range": {
                        "vmin": min(sample_payload["vmin"] for sample_payload in samples.values()),
                        "vmax": max(sample_payload["vmax"] for sample_payload in samples.values()),
                    },
                    "samples": samples,
                }
            except ValueError as exc:
                json_response(self, 404, {"error": str(exc)})
                return
            except RuntimeError as exc:
                json_response(self, 503, {"error": str(exc)})
                return

            json_response(self, 200, payload)
            return

        if parsed.path == "/api/dotplot":
            gene = (params.get("gene") or [""])[0].strip()
            if not gene:
                json_response(self, 400, {"error": "missing gene"})
                return

            try:
                payload = load_dotplot_from_db(dataset, gene)
            except ValueError as exc:
                json_response(self, 404, {"error": str(exc)})
                return
            except RuntimeError as exc:
                json_response(self, 503, {"error": str(exc)})
                return
            except sqlite3.Error as exc:
                json_response(self, 500, {"error": f"dotplot SQLite error: {exc}"})
                return

            json_response(self, 200, payload)
            return

        super().do_GET()


def main():
    parser = argparse.ArgumentParser(description="Run the spatial expression web viewer.")
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", type=int, default=8000)
    args = parser.parse_args()

    CACHE_ROOT.mkdir(parents=True, exist_ok=True)
    for dataset in load_dataset_catalog()["datasets"]:
        dataset["data_root"].mkdir(parents=True, exist_ok=True)
        dataset["cache_root"].mkdir(parents=True, exist_ok=True)

    server = ThreadingHTTPServer((args.host, args.port), Handler)
    print(f"Serving http://{args.host}:{args.port}")
    server.serve_forever()


if __name__ == "__main__":
    main()
