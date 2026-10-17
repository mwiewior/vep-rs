#!/usr/bin/env python3
"""Convert a Perl VEP cache's co-located variants to vep-rs's JSON variation shards.

A Perl VEP cache with `var_type tabix` keeps each chromosome's known variants in
`<chr>/all_vars.gz`, one tab-separated row per variant in the column order its
`info.txt` names under `variation_cols`. This writes them as the shards vep-rs's
`load_all_variations` reads, `<output_dir>/variations/<chr>/<start>-<end>.json`
(1 Mb bins), each a JSON array of objects keyed by those column names, with the
cache's `.` (unset) omitted. Together with `storable_to_json.pl`, which converts the
transcripts, it gives a JSON cache with the same transcripts and co-located variants
as the Storable cache, for `--check_existing`, the `--af*` flags and the rest of
`--everything`.

Usage:
  all_vars_to_json.py <cache_dir> <output_dir> [chr1,chr2,...] [--gzip]

`cache_dir` is the Storable cache (for example `~/.vep/homo_sapiens/116_GRCh38`);
`output_dir` is the JSON cache directory `storable_to_json.pl` wrote, or a new one.
Without a chromosome list every chromosome directory holding an `all_vars.gz` is
converted. `--gzip` writes `.json.gz` shards, about a tenth of the size; vep-rs
releases up to 0.4.0 read only plain `.json` variation shards.
"""

from __future__ import annotations

import argparse
import gzip
import json
import sys
from collections import defaultdict
from pathlib import Path

BIN = 1_000_000


def variation_cols(cache_dir: Path) -> list[str]:
    for line in open(cache_dir / "info.txt"):
        key, _, value = line.rstrip("\n").partition("\t")
        if key == "variation_cols":
            return value.split(",")
    sys.exit(f"{cache_dir}/info.txt names no variation_cols; is this a tabix variation cache?")


def convert(cache_dir: Path, chrom: str, output_dir: Path, cols: list[str], gz: bool) -> int:
    bins: dict[int, list[dict]] = defaultdict(list)
    with gzip.open(cache_dir / chrom / "all_vars.gz", "rt") as fh:
        for line in fh:
            row = {c: v for c, v in zip(cols, line.rstrip("\n").split("\t")) if v not in (".", "") and c != "chr"}
            bins[(int(row["start"]) - 1) // BIN * BIN + 1].append(row)
    out = output_dir / "variations" / chrom
    out.mkdir(parents=True, exist_ok=True)
    n = 0
    for start, rows in sorted(bins.items()):
        name = f"{start}-{start + BIN - 1}.json"
        if gz:
            with gzip.open(out / f"{name}.gz", "wt", compresslevel=6) as fo:
                json.dump(rows, fo)
        else:
            with open(out / name, "w") as fo:
                json.dump(rows, fo)
        n += len(rows)
    return n


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("cache_dir", type=Path)
    p.add_argument("output_dir", type=Path)
    p.add_argument("chromosomes", nargs="?", help="comma-separated list (default: every chromosome with an all_vars.gz)")
    p.add_argument("--gzip", action="store_true", help="write .json.gz shards")
    a = p.parse_args(argv)
    cols = variation_cols(a.cache_dir)
    chroms = (
        a.chromosomes.split(",")
        if a.chromosomes
        else sorted(d.name for d in a.cache_dir.iterdir() if (d / "all_vars.gz").exists())
    )
    for chrom in chroms:
        n = convert(a.cache_dir, chrom, a.output_dir, cols, a.gzip)
        print(f"chr{chrom}: {n:,} variants", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
