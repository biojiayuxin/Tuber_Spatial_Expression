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
DATA_ROOT = WEB_ROOT / "data"
CACHE_ROOT = WEB_ROOT / "cache"
DB_PATH = DATA_ROOT / "expression.sqlite"
CACHE_VERSION = 2
DOTPLOT_CLUSTER_COLUMN = "seurat_clusters"

SAMPLES = {
    "S1": ROOT / "S1_all_genes.csv",
    "S2": ROOT / "S2_all_genes.csv",
}


def connect_db():
    return sqlite3.connect(f"file:{DB_PATH}?mode=ro", uri=True)


def json_response(handler, status, payload):
    body = json.dumps(payload, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
    handler.send_response(status)
    handler.send_header("Content-Type", "application/json; charset=utf-8")
    handler.send_header("Content-Length", str(len(body)))
    handler.send_header("Cache-Control", "no-store")
    handler.end_headers()
    handler.wfile.write(body)


def read_gene_names(csv_path):
    with open(csv_path, newline="") as handle:
        reader = csv.reader(handle)
        header = next(reader)
    return header[1:]


def get_gene_names():
    if DB_PATH.exists():
        try:
            with connect_db() as conn:
                rows = conn.execute("SELECT gene FROM genes ORDER BY gene_id").fetchall()
            if rows:
                return [row[0] for row in rows]
        except sqlite3.Error as exc:
            print(f"SQLite gene list unavailable, falling back to JSON/CSV: {exc}", flush=True)

    genes_path = DATA_ROOT / "genes.json"
    if genes_path.exists():
        return json.loads(genes_path.read_text(encoding="utf-8"))

    genes = read_gene_names(SAMPLES["S1"])
    DATA_ROOT.mkdir(parents=True, exist_ok=True)
    genes_path.write_text(json.dumps(genes, separators=(",", ":")), encoding="utf-8")
    return genes


def cache_path(sample, gene):
    safe_gene = gene.replace("/", "_").replace("\\", "_")
    return CACHE_ROOT / f"{sample}_{safe_gene}.json"


def load_expression_from_db(sample, gene):
    with connect_db() as conn:
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


def load_expression_from_csv(sample, gene):
    path = cache_path(sample, gene)
    if path.exists():
        payload = json.loads(path.read_text(encoding="utf-8"))
        if (
            payload.get("formatVersion") == CACHE_VERSION
            and "vmin" in payload
            and "vmax" in payload
            and "values" in payload
        ):
            return payload

    csv_path = SAMPLES[sample]
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

    CACHE_ROOT.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, separators=(",", ":")), encoding="utf-8")
    return payload


def load_expression(sample, gene):
    if DB_PATH.exists():
        try:
            return load_expression_from_db(sample, gene)
        except sqlite3.Error as exc:
            print(f"SQLite expression lookup failed, falling back to CSV/cache: {exc}", flush=True)

    return load_expression_from_csv(sample, gene)


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


def load_dotplot_from_db(gene):
    if not DB_PATH.exists():
        raise RuntimeError("dotplot database not found; build web_viewer/data/expression.sqlite first")

    with connect_db() as conn:
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


def get_json_document(name, fallback_path):
    if DB_PATH.exists():
        try:
            with connect_db() as conn:
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

        if clean_path == "":
            return str(STATIC_ROOT / "index.html")

        return str(STATIC_ROOT / clean_path)

    def do_GET(self):
        parsed = urlparse(self.path)

        if parsed.path == "/api/genes":
            genes = get_gene_names()
            json_response(self, 200, {"genes": genes})
            return

        if parsed.path == "/api/replicates":
            replicates = get_json_document("replicates", DATA_ROOT / "replicates.json")
            if replicates is None:
                json_response(self, 404, {"error": "replicates not found"})
                return
            json_response(self, 200, replicates)
            return

        if parsed.path == "/api/gene":
            params = parse_qs(parsed.query)
            gene = (params.get("gene") or [""])[0].strip()
            if not gene:
                json_response(self, 400, {"error": "missing gene"})
                return

            genes = get_gene_names()
            if gene not in genes:
                json_response(self, 404, {"error": f"{gene} not found"})
                return

            try:
                samples = {
                    sample: load_expression(sample, gene)
                    for sample in SAMPLES
                }
                payload = {
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

            json_response(self, 200, payload)
            return

        if parsed.path == "/api/dotplot":
            params = parse_qs(parsed.query)
            gene = (params.get("gene") or [""])[0].strip()
            if not gene:
                json_response(self, 400, {"error": "missing gene"})
                return

            try:
                payload = load_dotplot_from_db(gene)
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

    DATA_ROOT.mkdir(parents=True, exist_ok=True)
    CACHE_ROOT.mkdir(parents=True, exist_ok=True)

    server = ThreadingHTTPServer((args.host, args.port), Handler)
    print(f"Serving http://{args.host}:{args.port}")
    server.serve_forever()


if __name__ == "__main__":
    main()
