#!/usr/bin/env python3
"""
merge_window_counts.py

Collect the per-library window coverage files produced by BEDTOOLS_COVERAGE into a
single bins x libraries count matrix, and build the companion per-library metrics
table (spike-in pairs, scale factors, control pairing) that the downstream
normalisation needs.

Library identity is resolved by matching each coverage filename against the `id`
column of the pipeline metadata, longest match wins. Unmatched or ambiguous files
are a hard error: a silently mis-assigned library is worse than a failed run.

Outputs
-------
window_counts.csv.gz   chrom,start,end,<library_id> ... raw counts per window. These are
                       bedtools coverage counts of BAM records, i.e. each mate of a
                       pair is counted (about 2 x fragments); the column
                       window_fragments in the metrics keeps its historical name.
library_metrics.csv    one row per library, with everything needed to build offsets
"""

import argparse
import csv
import glob
import gzip
import os
import re
import sys


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--coverage-dir", required=True,
                   help="directory holding the BEDTOOLS_COVERAGE *.bed outputs")
    p.add_argument("--spikein-log-dir", default=None,
                   help="directory holding the spike-in bowtie2 *.log files")
    p.add_argument("--metadata", required=True,
                   help="pipeline metadata CSV (needs at least an `id` column)")
    p.add_argument("--out-counts", default="window_counts.csv.gz")
    p.add_argument("--out-metrics", default="library_metrics.csv")
    p.add_argument("--normalisation-c", type=float, default=10000.0,
                   help="numerator of the spike-in scale factor (params.normalisation_c)")
    p.add_argument("--count-column", type=int, default=3,
                   help="0-based column of the read count in the coverage BED "
                        "(bedtools coverage on a 3-column window BED puts it at 3)")
    return p.parse_args()


def read_metadata(path):
    with open(path, newline="", encoding="utf-8-sig") as fh:
        rows = [{k: (v.strip() if isinstance(v, str) else v) for k, v in r.items()}
                for r in csv.DictReader(fh)]
    if not rows:
        sys.exit(f"ERROR: metadata {path} has no rows")
    if "id" not in rows[0]:
        sys.exit(f"ERROR: metadata {path} has no `id` column")
    ids = [r["id"] for r in rows]
    if len(set(ids)) != len(ids):
        dupes = sorted({i for i in ids if ids.count(i) > 1})
        sys.exit(f"ERROR: duplicate ids in metadata: {dupes}")
    return rows, ids


def assign_library(fname, ids):
    """Longest id that the basename starts with. Hard error on no/ambiguous match."""
    base = os.path.basename(fname)
    hits = [i for i in ids if base == i or base.startswith(i + "_") or base.startswith(i + ".")]
    if not hits:
        return None, f"no metadata id matches '{base}'"
    hits.sort(key=len, reverse=True)
    if len(hits) > 1 and len(hits[0]) == len(hits[1]):
        return None, f"'{base}' matches several ids equally well: {hits[:2]}"
    return hits[0], None


def parse_spikein_log(path):
    """Concordantly aligned pairs, exactly as Bam_Norm.nf counts them."""
    total = 0
    pat = re.compile(r"^\s*(\d+)\s+\([\d.]+%\)\s+aligned concordantly "
                     r"(?:exactly 1 time|>1 times)")
    with open(path) as fh:
        for line in fh:
            m = pat.match(line)
            if m:
                total += int(m.group(1))
    return total


def read_coverage(path, count_col):
    coords, counts = [], []
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.rstrip("\n").split("\t")
            coords.append((f[0], int(f[1]), int(f[2])))
            counts.append(int(float(f[count_col])))
    return coords, counts


def main():
    a = parse_args()
    meta_rows, ids = read_metadata(a.metadata)

    cov_files = sorted(glob.glob(os.path.join(a.coverage_dir, "*.bed")) +
                       glob.glob(os.path.join(a.coverage_dir, "*.bed.gz")))
    if not cov_files:
        sys.exit(f"ERROR: no *.bed found in {a.coverage_dir}")

    cov_by_lib, problems = {}, []
    for f in cov_files:
        lib, err = assign_library(f, ids)
        if err:
            problems.append(err)
            continue
        if lib in cov_by_lib:
            problems.append(f"two coverage files map to library '{lib}': "
                            f"{os.path.basename(cov_by_lib[lib])}, {os.path.basename(f)}")
            continue
        cov_by_lib[lib] = f
    if problems:
        sys.exit("ERROR: coverage files could not be assigned to libraries:\n  - "
                 + "\n  - ".join(problems))

    missing = [i for i in ids if i not in cov_by_lib]
    if missing:
        print(f"WARNING: {len(missing)} metadata libraries have no coverage file "
              f"and are dropped: {missing}", file=sys.stderr)
    libs = [i for i in ids if i in cov_by_lib]

    # --- counts matrix, coordinates taken from the first library and enforced on the rest
    ref_coords, mat = None, {}
    for lib in libs:
        coords, counts = read_coverage(cov_by_lib[lib], a.count_column)
        if ref_coords is None:
            ref_coords = coords
        elif coords != ref_coords:
            sys.exit(f"ERROR: window coordinates of '{lib}' differ from "
                     f"'{libs[0]}'. All libraries must use the same window BED.")
        mat[lib] = counts
    n_win = len(ref_coords)

    out_counts = a.out_counts
    opener = gzip.open if out_counts.endswith(".gz") else open
    with opener(out_counts, "wt", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["chrom", "start", "end"] + libs)
        for k, (c, s, e) in enumerate(ref_coords):
            w.writerow([c, s, e] + [mat[lib][k] for lib in libs])

    # --- spike-in counts
    spike = {}
    if a.spikein_log_dir and os.path.isdir(a.spikein_log_dir):
        for f in sorted(glob.glob(os.path.join(a.spikein_log_dir, "*.log"))):
            lib, err = assign_library(f, ids)
            if err:
                print(f"WARNING: spike-in log skipped ({err})", file=sys.stderr)
                continue
            spike[lib] = parse_spikein_log(f)

    # --- metrics table
    by_id = {r["id"]: r for r in meta_rows}
    controls = {r["id"]: r.get("control", "") or "" for r in meta_rows}
    is_control = set(v for v in controls.values() if v)

    extra = [c for c in ("group", "replicate", "control") if c in meta_rows[0]]
    header = (["library_id"] + extra +
              ["is_control", "n_windows", "window_fragments",
               "spikein_pairs", "spikein_scale_factor"])
    with open(a.out_metrics, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(header)
        for lib in libs:
            tot = sum(mat[lib])
            sp = spike.get(lib, "")
            sf = (a.normalisation_c / sp) if isinstance(sp, int) and sp > 0 else ""
            w.writerow([lib] + [by_id[lib].get(c, "") for c in extra] +
                       [int(lib in is_control), n_win, tot, sp, sf])

    n_spike = sum(1 for lib in libs if isinstance(spike.get(lib), int) and spike[lib] > 0)
    print(f"merged {len(libs)} libraries x {n_win} windows -> {out_counts}")
    print(f"spike-in pairs recovered for {n_spike}/{len(libs)} libraries")
    print(f"control pairing declared for "
          f"{sum(1 for l in libs if controls.get(l))}/{len(libs)} libraries")


if __name__ == "__main__":
    main()
