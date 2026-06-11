#!/usr/bin/env python3
import argparse
import csv
import json
import sqlite3
import time
from pathlib import Path

import numpy as np
import pandas as pd


SCHEMA_VERSION = 1


def connect_db(path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(path)
    conn.execute("PRAGMA journal_mode=OFF")
    conn.execute("PRAGMA synchronous=OFF")
    conn.execute("PRAGMA temp_store=MEMORY")
    conn.execute("PRAGMA cache_size=-200000")
    return conn


def create_schema(conn):
    conn.executescript(
        """
        DROP TABLE IF EXISTS metadata;
        DROP TABLE IF EXISTS samples;
        DROP TABLE IF EXISTS genes;
        DROP TABLE IF EXISTS sample_genes;
        DROP TABLE IF EXISTS expression_values;
        DROP TABLE IF EXISTS json_documents;

        CREATE TABLE metadata (
          key TEXT PRIMARY KEY,
          value TEXT NOT NULL
        );

        CREATE TABLE samples (
          sample TEXT PRIMARY KEY,
          source_path TEXT NOT NULL,
          cell_count INTEGER NOT NULL
        );

        CREATE TABLE genes (
          gene_id INTEGER PRIMARY KEY,
          gene TEXT NOT NULL UNIQUE
        );

        CREATE TABLE sample_genes (
          sample TEXT NOT NULL,
          gene_id INTEGER NOT NULL,
          vmin REAL NOT NULL,
          vmax REAL NOT NULL,
          nonzero INTEGER NOT NULL,
          PRIMARY KEY (sample, gene_id)
        );

        CREATE TABLE expression_values (
          sample TEXT NOT NULL,
          gene_id INTEGER NOT NULL,
          cell_id INTEGER NOT NULL,
          value REAL NOT NULL
        );

        CREATE TABLE json_documents (
          name TEXT PRIMARY KEY,
          payload TEXT NOT NULL
        );
        """
    )
    conn.execute(
        "INSERT INTO metadata(key, value) VALUES (?, ?)",
        ("schema_version", str(SCHEMA_VERSION)),
    )


def create_indexes(conn):
    conn.executescript(
        """
        CREATE INDEX IF NOT EXISTS idx_expression_gene
          ON expression_values(sample, gene_id, cell_id);
        CREATE INDEX IF NOT EXISTS idx_genes_gene ON genes(gene);
        """
    )


def read_header(csv_path):
    with open(csv_path, newline="") as handle:
        reader = csv.reader(handle)
        header = next(reader)
    if not header or header[0] != "cell":
        raise ValueError(f"{csv_path} must start with a 'cell' column")
    return header[1:]


def insert_genes(conn, genes):
    rows = [(index + 1, gene) for index, gene in enumerate(genes)]
    conn.executemany(
        "INSERT INTO genes(gene_id, gene) VALUES (?, ?)",
        rows,
    )
    return {gene: index + 1 for index, gene in enumerate(genes)}


def import_expression_csv(conn, sample, csv_path, genes, gene_ids, chunk_rows):
    csv_path = Path(csv_path)
    header_genes = read_header(csv_path)
    if header_genes != genes:
        raise ValueError(
            f"{csv_path} gene header does not match the first imported sample"
        )

    gene_count = len(genes)
    gene_id_array = np.array([gene_ids[gene] for gene in genes], dtype=np.int64)
    vmin = np.full(gene_count, np.inf, dtype=np.float64)
    vmax = np.full(gene_count, -np.inf, dtype=np.float64)
    nonzero = np.zeros(gene_count, dtype=np.int64)
    cell_count = 0
    inserted = 0
    start = time.time()

    print(f"{sample}: importing {csv_path}", flush=True)
    for chunk_index, chunk in enumerate(pd.read_csv(csv_path, chunksize=chunk_rows), 1):
        cells = chunk["cell"].to_numpy(dtype=np.int64, copy=False)
        values = chunk.iloc[:, 1:].to_numpy(dtype=np.float64, copy=False)

        if values.size:
            vmin = np.minimum(vmin, values.min(axis=0))
            vmax = np.maximum(vmax, values.max(axis=0))

        row_index, col_index = np.nonzero(values > 0)
        if len(row_index):
            nonzero += np.bincount(col_index, minlength=gene_count)
            rows = zip(
                [sample] * len(row_index),
                gene_id_array[col_index].astype(int).tolist(),
                cells[row_index].astype(int).tolist(),
                values[row_index, col_index].astype(float).tolist(),
            )
            conn.executemany(
                """
                INSERT INTO expression_values(sample, gene_id, cell_id, value)
                VALUES (?, ?, ?, ?)
                """,
                rows,
            )
            inserted += len(row_index)

        cell_count += len(chunk)
        if chunk_index % 10 == 0:
            elapsed = time.time() - start
            print(
                f"{sample}: {cell_count:,} rows, {inserted:,} nonzero values "
                f"({elapsed:.1f}s)",
                flush=True,
            )

    vmin[~np.isfinite(vmin)] = 0
    vmax[~np.isfinite(vmax)] = 0
    stats_rows = [
        (sample, int(gene_id_array[i]), float(vmin[i]), float(vmax[i]), int(nonzero[i]))
        for i in range(gene_count)
    ]
    conn.executemany(
        """
        INSERT INTO sample_genes(sample, gene_id, vmin, vmax, nonzero)
        VALUES (?, ?, ?, ?, ?)
        """,
        stats_rows,
    )
    conn.execute(
        "INSERT INTO samples(sample, source_path, cell_count) VALUES (?, ?, ?)",
        (sample, str(csv_path), int(cell_count)),
    )
    print(
        f"{sample}: imported {cell_count:,} rows and {inserted:,} nonzero values",
        flush=True,
    )


def import_json_document(conn, name, path):
    path = Path(path)
    if not path.exists():
        print(f"Skipping missing JSON document: {path}", flush=True)
        return
    payload = json.loads(path.read_text(encoding="utf-8"))
    conn.execute(
        "INSERT INTO json_documents(name, payload) VALUES (?, ?)",
        (name, json.dumps(payload, ensure_ascii=False, separators=(",", ":"))),
    )
    print(f"Imported JSON document {name} from {path}", flush=True)


def build_database(args):
    out_path = Path(args.out)
    if out_path.exists():
        if not args.replace:
            raise FileExistsError(f"{out_path} already exists; use --replace to rebuild it")
        out_path.unlink()

    conn = connect_db(out_path)
    try:
        create_schema(conn)
        s1_genes = read_header(args.s1_csv)
        s2_genes = read_header(args.s2_csv)
        if s1_genes != s2_genes:
            raise ValueError("S1 and S2 CSV files must have the same gene columns")

        gene_ids = insert_genes(conn, s1_genes)
        conn.commit()

        import_expression_csv(conn, "S1", args.s1_csv, s1_genes, gene_ids, args.chunk_rows)
        conn.commit()
        import_expression_csv(conn, "S2", args.s2_csv, s1_genes, gene_ids, args.chunk_rows)
        conn.commit()

        import_json_document(conn, "replicates", args.replicates_json)
        conn.commit()

        print("Creating SQLite indexes...", flush=True)
        create_indexes(conn)
        conn.execute("ANALYZE")
        conn.commit()
    finally:
        conn.close()

    print(f"Wrote {out_path}", flush=True)


def main():
    parser = argparse.ArgumentParser(
        description="Build the SQLite data source used by web_viewer/server.py."
    )
    parser.add_argument("--s1-csv", default="S1_all_genes.csv")
    parser.add_argument("--s2-csv", default="S2_all_genes.csv")
    parser.add_argument("--replicates-json", default="web_viewer/data/replicates.json")
    parser.add_argument("--out", default="web_viewer/data/expression.sqlite")
    parser.add_argument("--chunk-rows", type=int, default=256)
    parser.add_argument("--replace", action="store_true")
    args = parser.parse_args()

    build_database(args)


if __name__ == "__main__":
    main()
