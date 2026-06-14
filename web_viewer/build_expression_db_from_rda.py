#!/usr/bin/env python3
import argparse
import csv
import json
import sqlite3
import subprocess
import tempfile
from pathlib import Path


SCHEMA_VERSION = 1


def r_literal(value):
    if value is None:
        return "NULL"
    return json.dumps(str(value))


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


def build_r_script(args, genes_path, stats_path, values_path, metadata_path):
    sample_rules = "list(" + ",".join(
        f"{json.dumps(sample)}={json.dumps(prefix)}"
        for sample, prefix in parse_sample_rules(args.sample_rule)
    ) + ")"
    return f"""
    suppressPackageStartupMessages(library(Matrix))

    rda_path <- {r_literal(Path(args.rda).resolve())}
    object_name <- {r_literal(args.object_name)}
    assay_name <- {r_literal(args.assay)}
    layer_name <- {r_literal(args.layer)}
    sample_rules <- {sample_rules}
    genes_path <- {r_literal(genes_path)}
    stats_path <- {r_literal(stats_path)}
    values_path <- {r_literal(values_path)}
    metadata_path <- {r_literal(metadata_path)}

    load(rda_path)
    if (!exists(object_name)) {{
      stop(sprintf("Object '%s' was not found in %s", object_name, rda_path))
    }}
    obj <- get(object_name)
    assays <- methods::slot(obj, "assays")
    if (!(assay_name %in% names(assays))) {{
      stop(sprintf("Assay '%s' was not found", assay_name))
    }}
    assay <- assays[[assay_name]]
    expr <- methods::slot(assay, layer_name)
    if (!inherits(expr, "sparseMatrix")) {{
      expr <- as(expr, "dgCMatrix")
    }}
    expr <- as(expr, "dgCMatrix")

    meta <- methods::slot(obj, "meta.data")
    cells <- colnames(expr)
    meta_index <- match(cells, rownames(meta))
    if (any(is.na(meta_index))) {{
      stop("Some expression matrix cells are missing from obj@meta.data row names")
    }}
    if (!("orig.ident" %in% colnames(meta))) {{
      stop("Column 'orig.ident' was not found in obj@meta.data")
    }}
    reps <- as.character(meta$orig.ident[meta_index])
    cell_ids <- suppressWarnings(as.integer(sub(".*cell_", "", cells)))
    if (any(is.na(cell_ids))) {{
      stop("Could not parse integer cell ids from expression matrix column names")
    }}

    sample_for_rep <- function(rep) {{
      for (sample in names(sample_rules)) {{
        if (startsWith(tolower(rep), tolower(sample_rules[[sample]]))) {{
          return(sample)
        }}
      }}
      return(NA_character_)
    }}
    samples <- vapply(reps, sample_for_rep, character(1))
    if (any(is.na(samples))) {{
      stop(sprintf("Unmapped orig.ident values: %s", paste(unique(reps[is.na(samples)]), collapse=", ")))
    }}

    genes_df <- data.frame(
      gene_id = seq_len(nrow(expr)),
      gene = rownames(expr),
      check.names = FALSE
    )
    write.table(genes_df, genes_path, sep="\\t", quote=FALSE, row.names=FALSE)

    if (file.exists(stats_path)) file.remove(stats_path)
    if (file.exists(values_path)) file.remove(values_path)

    sample_order <- names(sample_rules)
    for (sample in sample_order) {{
      sample_cols <- which(samples == sample)
      if (!length(sample_cols)) next
      submat <- expr[, sample_cols, drop=FALSE]
      cell_count <- length(sample_cols)
      sm <- summary(submat)
      nonzero <- tabulate(sm$i, nbins=nrow(expr))
      vmin <- numeric(nrow(expr))
      vmax <- numeric(nrow(expr))
      if (nrow(sm) > 0) {{
        min_values <- tapply(sm$x, sm$i, min)
        max_values <- tapply(sm$x, sm$i, max)
        min_index <- as.integer(names(min_values))
        max_index <- as.integer(names(max_values))
        vmin[min_index[nonzero[min_index] == cell_count]] <- as.numeric(min_values[nonzero[min_index] == cell_count])
        vmax[max_index] <- as.numeric(max_values)
      }}
      stats_df <- data.frame(
        sample = sample,
        gene_id = seq_len(nrow(expr)),
        vmin = as.numeric(vmin),
        vmax = as.numeric(vmax),
        nonzero = as.integer(nonzero),
        check.names = FALSE
      )
      write.table(
        stats_df,
        stats_path,
        sep="\\t",
        quote=FALSE,
        row.names=FALSE,
        col.names=!file.exists(stats_path),
        append=file.exists(stats_path)
      )

      values_df <- data.frame(
        sample = sample,
        gene_id = as.integer(sm$i),
        cell_id = as.integer(cell_ids[sample_cols[sm$j]]),
        value = as.numeric(sm$x),
        check.names = FALSE
      )
      write.table(
        values_df,
        values_path,
        sep="\\t",
        quote=FALSE,
        row.names=FALSE,
        col.names=!file.exists(values_path),
        append=file.exists(values_path)
      )
      message(sprintf("%s: exported %d cells and %d nonzero values", sample, cell_count, nrow(values_df)))
    }}

    sample_counts <- as.integer(tabulate(match(samples, sample_order), nbins=length(sample_order)))
    metadata_df <- data.frame(
      key = c("source_rda", "object_name", "assay", "layer"),
      value = c(rda_path, object_name, assay_name, layer_name),
      check.names = FALSE
    )
    samples_df <- data.frame(
      sample = sample_order,
      source_path = rda_path,
      cell_count = sample_counts,
      check.names = FALSE
    )
    write.table(
      metadata_df,
      metadata_path,
      sep="\\t",
      quote=FALSE,
      row.names=FALSE,
      col.names=TRUE
    )
    write.table(
      samples_df,
      paste0(metadata_path, ".samples"),
      sep="\\t",
      quote=FALSE,
      row.names=FALSE,
      col.names=TRUE
    )
    """


def parse_sample_rules(values):
    rules = []
    for value in values:
        if "=" not in value:
            raise ValueError(f"sample rule must be SAMPLE=prefix, got {value!r}")
        sample, prefix = value.split("=", 1)
        sample = sample.strip()
        prefix = prefix.strip()
        if not sample or not prefix:
            raise ValueError(f"sample rule must be SAMPLE=prefix, got {value!r}")
        rules.append((sample, prefix))
    return rules


def run_r_export(args, paths):
    script_path = paths["script"]
    script_path.write_text(
        build_r_script(
            args,
            paths["genes"],
            paths["stats"],
            paths["values"],
            paths["metadata"],
        ),
        encoding="utf-8",
    )
    subprocess.run(["Rscript", str(script_path)], check=True)


def insert_genes(conn, genes_path):
    with open(genes_path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = [(int(row["gene_id"]), row["gene"]) for row in reader]
    conn.executemany("INSERT INTO genes(gene_id, gene) VALUES (?, ?)", rows)
    return len(rows)


def insert_tsv(conn, table, columns, path, batch_size):
    total = 0
    batch = []
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            batch.append(tuple(convert_value(row[column]) for column in columns))
            if len(batch) >= batch_size:
                insert_batch(conn, table, columns, batch)
                total += len(batch)
                conn.commit()
                print(f"{table}: imported {total:,} rows", flush=True)
                batch.clear()
    if batch:
        insert_batch(conn, table, columns, batch)
        total += len(batch)
        conn.commit()
    return total


def convert_value(value):
    if value == "":
        return value
    try:
        if "." not in value and "e" not in value.lower():
            return int(value)
    except ValueError:
        return value
    try:
        return float(value)
    except ValueError:
        return value


def insert_batch(conn, table, columns, rows):
    placeholders = ",".join(["?"] * len(columns))
    column_sql = ",".join(columns)
    conn.executemany(
        f"INSERT INTO {table}({column_sql}) VALUES ({placeholders})",
        rows,
    )


def import_metadata(conn, metadata_path):
    with open(metadata_path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = [(row["key"], row["value"]) for row in reader]
    conn.executemany("INSERT OR REPLACE INTO metadata(key, value) VALUES (?, ?)", rows)

    with open(str(metadata_path) + ".samples", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        sample_rows = [
            (row["sample"], row["source_path"], int(row["cell_count"]))
            for row in reader
            if int(row["cell_count"]) > 0
        ]
    conn.executemany(
        "INSERT INTO samples(sample, source_path, cell_count) VALUES (?, ?, ?)",
        sample_rows,
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

    with tempfile.TemporaryDirectory(dir=args.tmp_dir) as temp_dir:
        temp_dir = Path(temp_dir)
        paths = {
            "genes": temp_dir / "genes.tsv",
            "stats": temp_dir / "sample_genes.tsv",
            "values": temp_dir / "expression_values.tsv",
            "metadata": temp_dir / "metadata.tsv",
            "script": temp_dir / "export_expression.R",
        }
        run_r_export(args, paths)

        conn = connect_db(out_path)
        try:
            create_schema(conn)
            gene_count = insert_genes(conn, paths["genes"])
            conn.commit()
            stats_total = insert_tsv(
                conn,
                "sample_genes",
                ["sample", "gene_id", "vmin", "vmax", "nonzero"],
                paths["stats"],
                args.batch_size,
            )
            values_total = insert_tsv(
                conn,
                "expression_values",
                ["sample", "gene_id", "cell_id", "value"],
                paths["values"],
                args.batch_size,
            )
            import_metadata(conn, paths["metadata"])
            import_json_document(conn, "replicates", args.replicates_json)
            conn.commit()
            print("Creating SQLite indexes...", flush=True)
            create_indexes(conn)
            conn.execute("ANALYZE")
            conn.commit()
        finally:
            conn.close()

    print(
        f"Wrote {out_path} with {gene_count:,} genes, {stats_total:,} sample-gene rows, "
        f"and {values_total:,} expression values",
        flush=True,
    )


def main():
    parser = argparse.ArgumentParser(
        description="Build expression.sqlite directly from a Seurat Rda object."
    )
    parser.add_argument("--rda", required=True)
    parser.add_argument("--object-name", default="st")
    parser.add_argument("--assay", default="SCT")
    parser.add_argument("--layer", default="data")
    parser.add_argument(
        "--sample-rule",
        action="append",
        required=True,
        help="Map output sample id to orig.ident prefix, for example S1=s1_",
    )
    parser.add_argument("--replicates-json", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--batch-size", type=int, default=100000)
    parser.add_argument("--tmp-dir", default=None)
    parser.add_argument("--replace", action="store_true")
    args = parser.parse_args()

    build_database(args)


if __name__ == "__main__":
    main()
