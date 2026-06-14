#!/usr/bin/env python3
import argparse
import json
import subprocess
from pathlib import Path


def r_literal(value):
    return json.dumps(str(value))


def export_cell_ids(args):
    r_code = f"""
    load({r_literal(Path(args.rda).resolve())})
    obj <- get({r_literal(args.object_name)})
    meta <- methods::slot(obj, "meta.data")
    if (!("orig.ident" %in% colnames(meta))) {{
      stop("Column 'orig.ident' was not found in obj@meta.data")
    }}
    ids <- suppressWarnings(as.integer(sub(".*cell_", "", rownames(meta))))
    if (any(is.na(ids))) {{
      stop("Could not parse integer cell ids from metadata row names")
    }}
    keep <- startsWith(tolower(as.character(meta$orig.ident)), tolower({r_literal(args.rep_prefix)}))
    out <- {r_literal(Path(args.out).resolve())}
    dir.create(dirname(out), recursive=TRUE, showWarnings=FALSE)
    writeLines(c("cell_id", as.character(sort(unique(ids[keep])))), out)
    cat(sprintf("Wrote %d cell ids to %s\\n", sum(keep), out))
    """
    subprocess.run(["Rscript", "-e", r_code], check=True)


def main():
    parser = argparse.ArgumentParser(
        description="Export numeric cell ids from Seurat Rda metadata for a rep prefix."
    )
    parser.add_argument("--rda", required=True)
    parser.add_argument("--object-name", default="st")
    parser.add_argument("--rep-prefix", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    export_cell_ids(args)


if __name__ == "__main__":
    main()
