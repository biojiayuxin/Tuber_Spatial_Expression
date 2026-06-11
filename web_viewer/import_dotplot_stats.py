#!/usr/bin/env python3
import argparse
import csv
import json
import math
import sqlite3
import subprocess
import tempfile
from datetime import datetime, timezone
from pathlib import Path


DOTPLOT_CLUSTERS = "dotplot_clusters"
DOTPLOT_STATS = "dotplot_gene_cluster_stats"
IMPORT_CLUSTERS = "dotplot_clusters_import"
IMPORT_STATS = "dotplot_gene_cluster_stats_import"


def r_literal(value):
    if value is None:
        return "NULL"
    return json.dumps(str(value))


def connect_db(path):
    conn = sqlite3.connect(path)
    conn.execute("PRAGMA temp_store=MEMORY")
    conn.execute("PRAGMA cache_size=-200000")
    return conn


def table_exists(conn, name):
    row = conn.execute(
        "SELECT 1 FROM sqlite_master WHERE type = 'table' AND name = ?",
        (name,),
    ).fetchone()
    return row is not None


def load_gene_ids(conn):
    rows = conn.execute("SELECT gene_id, gene FROM genes").fetchall()
    if not rows:
        raise ValueError("No genes found in SQLite table 'genes'")
    return {gene: int(gene_id) for gene_id, gene in rows}


def create_import_schema(conn):
    conn.executescript(
        f"""
        DROP TABLE IF EXISTS {IMPORT_CLUSTERS};
        DROP TABLE IF EXISTS {IMPORT_STATS};

        CREATE TABLE {IMPORT_CLUSTERS} (
          cluster_id TEXT PRIMARY KEY,
          cluster_order INTEGER NOT NULL,
          cell_count INTEGER NOT NULL
        );

        CREATE TABLE {IMPORT_STATS} (
          gene_id INTEGER NOT NULL,
          cluster_id TEXT NOT NULL,
          avg_expr REAL NOT NULL,
          pct_expr REAL NOT NULL,
          expressing_count INTEGER NOT NULL,
          cell_count INTEGER NOT NULL,
          PRIMARY KEY (gene_id, cluster_id)
        );
        """
    )
    conn.commit()


def build_r_script(args, clusters_path, stats_path, metadata_path):
    return f"""
    suppressPackageStartupMessages({{
      library(Seurat)
      library(Matrix)
    }})

    rda_path <- {r_literal(Path(args.rda).resolve())}
    object_name <- {r_literal(args.object_name)}
    cluster_column <- {r_literal(args.cluster_column)}
    assay_name <- {r_literal(args.assay)}
    layer_name <- {r_literal(args.layer)}
    clusters_path <- {r_literal(clusters_path)}
    stats_path <- {r_literal(stats_path)}
    metadata_path <- {r_literal(metadata_path)}

    load(rda_path)
    if (!exists(object_name)) {{
      stop(sprintf("Object '%s' was not found in %s", object_name, rda_path))
    }}
    obj <- get(object_name)
    meta <- obj@meta.data
    if (!(cluster_column %in% colnames(meta))) {{
      stop(sprintf("Cluster column '%s' was not found in obj@meta.data", cluster_column))
    }}

    if (is.null(assay_name) || !nzchar(assay_name)) {{
      assay_name <- Seurat::DefaultAssay(obj)
    }}
    message(sprintf("Using assay '%s' and layer/slot '%s'", assay_name, layer_name))

    expr <- tryCatch(
      Seurat::GetAssayData(obj, assay = assay_name, layer = layer_name),
      error = function(layer_error) {{
        message("Layer lookup failed, retrying slot lookup: ", conditionMessage(layer_error))
        Seurat::GetAssayData(obj, assay = assay_name, slot = layer_name)
      }}
    )

    if (!inherits(expr, "sparseMatrix")) {{
      expr <- as(expr, "dgCMatrix")
    }}

    metadata_df <- data.frame(
      key = c("assay", "layer", "average_expression"),
      value = c(assay_name, layer_name, "Seurat DotPlot avg.exp: mean(expm1(data))"),
      check.names = FALSE
    )
    write.table(
      metadata_df,
      metadata_path,
      sep = "\\t",
      quote = FALSE,
      row.names = FALSE,
      col.names = TRUE
    )

    cells <- colnames(expr)
    meta_index <- match(cells, rownames(meta))
    if (any(is.na(meta_index))) {{
      stop("Some expression matrix cells are missing from obj@meta.data row names")
    }}

    cluster_values <- as.character(meta[[cluster_column]][meta_index])
    if (any(is.na(cluster_values)) || any(!nzchar(cluster_values))) {{
      stop(sprintf("Cluster column '%s' contains missing values", cluster_column))
    }}

    cluster_numeric <- suppressWarnings(as.integer(cluster_values))
    if (any(is.na(cluster_numeric))) {{
      stop(sprintf("Cluster column '%s' must contain numeric cluster ids", cluster_column))
    }}

    cluster_ids <- as.character(sort(unique(cluster_numeric)))
    cluster_counts <- tabulate(match(cluster_values, cluster_ids), nbins = length(cluster_ids))
    clusters_df <- data.frame(
      cluster_id = cluster_ids,
      cluster_order = seq_along(cluster_ids) - 1L,
      cell_count = as.integer(cluster_counts),
      check.names = FALSE
    )
    write.table(
      clusters_df,
      clusters_path,
      sep = "\\t",
      quote = FALSE,
      row.names = FALSE,
      col.names = TRUE
    )

    if (file.exists(stats_path)) {{
      file.remove(stats_path)
    }}

    gene_names <- rownames(expr)
    for (i in seq_along(cluster_ids)) {{
      cluster_id <- cluster_ids[[i]]
      cluster_cells <- which(cluster_values == cluster_id)
      cell_count <- length(cluster_cells)
      submat <- expr[, cluster_cells, drop = FALSE]
      expressing_count <- Matrix::rowSums(submat > 0)
      submat_linear <- submat
      submat_linear@x <- expm1(submat_linear@x)
      avg_expr <- Matrix::rowMeans(submat_linear)
      stats_df <- data.frame(
        gene = gene_names,
        cluster_id = cluster_id,
        avg_expr = as.numeric(avg_expr),
        pct_expr = as.numeric(expressing_count) / cell_count * 100,
        expressing_count = as.integer(expressing_count),
        cell_count = as.integer(cell_count),
        check.names = FALSE
      )
      write.table(
        stats_df,
        stats_path,
        sep = "\\t",
        quote = FALSE,
        row.names = FALSE,
        col.names = i == 1L,
        append = i > 1L
      )
      message(sprintf(
        "Cluster %s: %d cells, wrote %d genes",
        cluster_id,
        cell_count,
        length(gene_names)
      ))
    }}
    """


def run_r_export(args, clusters_path, stats_path, metadata_path, script_path):
    script_path.write_text(
        build_r_script(args, clusters_path, stats_path, metadata_path),
        encoding="utf-8",
    )
    subprocess.run(["Rscript", str(script_path)], check=True)


def read_r_metadata(metadata_path):
    metadata = {}
    with open(metadata_path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            metadata[row["key"]] = row["value"]
    return metadata


def import_clusters(conn, clusters_path):
    rows = []
    with open(clusters_path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            rows.append(
                (
                    row["cluster_id"],
                    int(row["cluster_order"]),
                    int(row["cell_count"]),
                )
            )

    if not rows:
        raise ValueError("R export produced no clusters")

    conn.executemany(
        f"""
        INSERT INTO {IMPORT_CLUSTERS}(cluster_id, cluster_order, cell_count)
        VALUES (?, ?, ?)
        """,
        rows,
    )
    conn.commit()
    return rows


def insert_stat_batch(conn, batch):
    conn.executemany(
        f"""
        INSERT INTO {IMPORT_STATS}(
          gene_id,
          cluster_id,
          avg_expr,
          pct_expr,
          expressing_count,
          cell_count
        )
        VALUES (?, ?, ?, ?, ?, ?)
        """,
        batch,
    )


def import_stats(conn, stats_path, gene_ids, batch_size):
    total = 0
    skipped = 0
    batch = []

    with open(stats_path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            gene_id = gene_ids.get(row["gene"])
            if gene_id is None:
                skipped += 1
                continue

            avg_expr = float(row["avg_expr"])
            pct_expr = float(row["pct_expr"])
            if not math.isfinite(avg_expr) or not math.isfinite(pct_expr):
                raise ValueError(f"Non-finite dotplot statistic for gene {row['gene']}")

            batch.append(
                (
                    gene_id,
                    row["cluster_id"],
                    avg_expr,
                    pct_expr,
                    int(float(row["expressing_count"])),
                    int(float(row["cell_count"])),
                )
            )

            if len(batch) >= batch_size:
                insert_stat_batch(conn, batch)
                total += len(batch)
                conn.commit()
                print(f"Imported {total:,} dotplot rows", flush=True)
                batch.clear()

    if batch:
        insert_stat_batch(conn, batch)
        total += len(batch)
        conn.commit()

    return total, skipped


def validate_import(conn, gene_count, clusters, args):
    cluster_count = len(clusters)
    expected_rows = gene_count * cluster_count
    actual_rows = conn.execute(f"SELECT COUNT(*) FROM {IMPORT_STATS}").fetchone()[0]
    if actual_rows != expected_rows:
        raise ValueError(
            f"Dotplot stats row count mismatch: expected {expected_rows:,}, got {actual_rows:,}"
        )

    cluster_total = sum(row[2] for row in clusters)
    if args.expect_cells is not None and cluster_total != args.expect_cells:
        raise ValueError(
            f"Cluster cell count mismatch: expected {args.expect_cells:,}, got {cluster_total:,}"
        )

    if args.expect_clusters is not None and cluster_count != args.expect_clusters:
        raise ValueError(
            f"Cluster count mismatch: expected {args.expect_clusters}, got {cluster_count}"
        )

    sample_total = conn.execute("SELECT COALESCE(SUM(cell_count), 0) FROM samples").fetchone()[0]
    if sample_total and cluster_total != int(sample_total):
        raise ValueError(
            f"Cluster cell total {cluster_total:,} does not match SQLite samples total {int(sample_total):,}"
        )


def swap_import_tables(conn, args, clusters, r_metadata):
    imported_at = datetime.now(timezone.utc).isoformat(timespec="seconds")
    cluster_total = sum(row[2] for row in clusters)

    conn.commit()
    try:
        conn.execute("BEGIN")
        conn.execute(f"DROP TABLE IF EXISTS {DOTPLOT_STATS}")
        conn.execute(f"DROP TABLE IF EXISTS {DOTPLOT_CLUSTERS}")
        conn.execute(f"ALTER TABLE {IMPORT_CLUSTERS} RENAME TO {DOTPLOT_CLUSTERS}")
        conn.execute(f"ALTER TABLE {IMPORT_STATS} RENAME TO {DOTPLOT_STATS}")
        conn.executemany(
            "INSERT OR REPLACE INTO metadata(key, value) VALUES (?, ?)",
            [
                ("dotplot_cluster_column", args.cluster_column),
                ("dotplot_source_rda", str(Path(args.rda))),
                ("dotplot_assay", r_metadata.get("assay", args.assay or "")),
                ("dotplot_layer", r_metadata.get("layer", args.layer)),
                ("dotplot_average_expression", r_metadata.get("average_expression", "Seurat DotPlot avg.exp")),
                ("dotplot_cluster_count", str(len(clusters))),
                ("dotplot_cell_count", str(cluster_total)),
                ("dotplot_imported_at", imported_at),
            ],
        )
        conn.commit()
    except Exception:
        conn.rollback()
        raise


def import_dotplot(args):
    db_path = Path(args.db)
    if not db_path.exists():
        raise FileNotFoundError(f"{db_path} does not exist; build expression.sqlite first")

    with connect_db(db_path) as conn:
        final_exists = table_exists(conn, DOTPLOT_CLUSTERS) or table_exists(conn, DOTPLOT_STATS)
        if final_exists and not args.replace:
            raise FileExistsError(
                "Dotplot tables already exist; use --replace to rebuild them"
            )
        gene_ids = load_gene_ids(conn)

    with tempfile.TemporaryDirectory(dir=args.tmp_dir) as temp_dir:
        temp_dir = Path(temp_dir)
        clusters_path = temp_dir / "dotplot_clusters.tsv"
        stats_path = temp_dir / "dotplot_gene_cluster_stats.tsv"
        metadata_path = temp_dir / "dotplot_metadata.tsv"
        script_path = temp_dir / "export_dotplot_stats.R"

        print("Exporting cluster dotplot statistics from Seurat...", flush=True)
        run_r_export(args, clusters_path, stats_path, metadata_path, script_path)

        with connect_db(db_path) as conn:
            r_metadata = read_r_metadata(metadata_path)
            create_import_schema(conn)
            clusters = import_clusters(conn, clusters_path)
            total, skipped = import_stats(conn, stats_path, gene_ids, args.batch_size)
            print(f"Imported {total:,} dotplot rows; skipped {skipped:,} rows not in SQLite genes", flush=True)
            validate_import(conn, len(gene_ids), clusters, args)
            swap_import_tables(conn, args, clusters, r_metadata)
            conn.execute("ANALYZE")
            conn.commit()

    print(
        f"Wrote {DOTPLOT_CLUSTERS} and {DOTPLOT_STATS} to {db_path}",
        flush=True,
    )


def main():
    parser = argparse.ArgumentParser(
        description="Import single-gene cluster dotplot statistics into expression.sqlite."
    )
    parser.add_argument("--rda", default="seurat_object.02-dims64.res0.6.Rda")
    parser.add_argument("--db", default="web_viewer/data/expression.sqlite")
    parser.add_argument("--object-name", default="st")
    parser.add_argument("--cluster-column", default="seurat_clusters")
    parser.add_argument("--assay", default=None)
    parser.add_argument("--layer", default="data")
    parser.add_argument("--batch-size", type=int, default=10000)
    parser.add_argument("--tmp-dir", default=None)
    parser.add_argument("--expect-cells", type=int, default=None)
    parser.add_argument("--expect-clusters", type=int, default=None)
    parser.add_argument("--replace", action="store_true")
    args = parser.parse_args()

    import_dotplot(args)


if __name__ == "__main__":
    main()
