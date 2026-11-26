#!/usr/bin/env python3
"""
make_test_data.py - regenerate the synthetic CUT&RUN test dataset used by -profile test.

Everything is simulated with a fixed seed, so the files are reproducible byte for byte
(gzip mtime is pinned to 0).

Genomes
  worm.fa   I (120 kb), II (100 kb), X (80 kb), MtDNA (6 kb). A 2 kb block of I is copied
            into II so that reads there are multi-mapping and removed by the MAPQ filter.
  ecoli.fa  one 60 kb contig standing in for the E. coli carry-over spike-in.

Libraries (2 x 50 bp, ~150 bp fragments)
  WT_K27me3_1   H3K27me3 domains on I:20-60 kb, II:30-50 kb and X:10-70 kb
  WT_noAb_1     uniform background (control of WT_K27me3_1)
  MUT_K27me3_1  same domains, the one on I at half strength
  MUT_noAb_1    uniform background (control of MUT_K27me3_1)
Each library also carries ~4 % MtDNA and a library-specific E. coli fraction, so the
spike-in scale factors differ between libraries.

Usage: python make_test_data.py [outdir]   (default: the directory of this script)
"""
import gzip, os, random, sys

SEED = 20261126
READ_LEN = 50
OUT = sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(os.path.abspath(__file__))
rng = random.Random(SEED)
COMP = str.maketrans("ACGT", "TGCA")


def rand_seq(n):
    return "".join(rng.choice("ACGT") for _ in range(n))


def revcomp(s):
    return s.translate(COMP)[::-1]


def write_fasta(path, recs):
    with open(path, "w") as fh:
        for name, seq in recs:
            fh.write(f">{name}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i + 60] + "\n")


worm = {"I": rand_seq(120_000), "II": rand_seq(100_000), "X": rand_seq(80_000), "MtDNA": rand_seq(6_000)}
# segmental duplication I:90,000-92,000 -> II:80,000-82,000 (multi-mapping block)
dup = worm["I"][90_000:92_000]
worm["II"] = worm["II"][:80_000] + dup + worm["II"][82_000:]
ecoli = {"ecoli": rand_seq(60_000)}
write_fasta(os.path.join(OUT, "worm.fa"), worm.items())
write_fasta(os.path.join(OUT, "ecoli.fa"), ecoli.items())

nuclear = ["I", "II", "X"]
nuc_len = {c: len(worm[c]) for c in nuclear}
domains = [("I", 20_000, 60_000), ("II", 30_000, 50_000), ("X", 10_000, 70_000)]


def frag_len():
    return max(60, min(500, int(rng.gauss(150, 35))))


def sample_frag(genome, chrom):
    L = frag_len()
    s = rng.randrange(0, len(genome[chrom]) - L)
    return genome[chrom][s:s + L]


def sample_background():
    tot = sum(nuc_len.values())
    r = rng.randrange(tot)
    for c in nuclear:
        if r < nuc_len[c]:
            return c
        r -= nuc_len[c]


def sample_domain(weights):
    tot = sum(w * (e - s) for (c, s, e), w in zip(domains, weights))
    r = rng.uniform(0, tot)
    for (c, s, e), w in zip(domains, weights):
        span = w * (e - s)
        if r < span:
            L = frag_len()
            st = rng.randrange(s, e - L)
            return worm[c][st:st + L]
        r -= span
    raise RuntimeError


def library(n_pairs, signal_frac, weights, ecoli_frac, mito_frac=0.04):
    frags = []
    for _ in range(n_pairs):
        u = rng.random()
        if u < ecoli_frac:
            frags.append(sample_frag(ecoli, "ecoli"))
        elif u < ecoli_frac + mito_frac:
            frags.append(sample_frag(worm, "MtDNA"))
        elif u < ecoli_frac + mito_frac + signal_frac:
            frags.append(sample_domain(weights))
        else:
            frags.append(sample_frag(worm, sample_background()))
    return frags


def write_fastq(prefix, frags):
    qual = "I" * READ_LEN
    paths = [os.path.join(OUT, f"{prefix}_R{i}.fastq.gz") for i in (1, 2)]
    fhs = [gzip.GzipFile(p, "wb", mtime=0) for p in paths]
    for k, f in enumerate(frags):
        if rng.random() < 0.5:   # random strand
            f = revcomp(f)
        r1, r2 = f[:READ_LEN], revcomp(f)[:READ_LEN]
        name = f"@{prefix}:{k}"
        fhs[0].write(f"{name}/1\n{r1}\n+\n{qual}\n".encode())
        fhs[1].write(f"{name}/2\n{r2}\n+\n{qual}\n".encode())
    for fh in fhs:
        fh.close()
    return [os.path.basename(p) for p in paths]


libs = [
    # id,           group,  ctrl,          pairs, signal, domain weights,  E. coli fraction
    ("WT_K27me3_1",  "WT",  "WT_noAb_1",   6000,  0.60,   (1.0, 1.0, 1.0), 0.020),
    ("WT_noAb_1",    "WT",  "",            3000,  0.00,   (1.0, 1.0, 1.0), 0.060),
    ("MUT_K27me3_1", "MUT", "MUT_noAb_1",  6000,  0.45,   (0.5, 1.0, 1.0), 0.030),
    ("MUT_noAb_1",   "MUT", "",            3000,  0.00,   (1.0, 1.0, 1.0), 0.050),
]
rows = ["id,group,replicate,control,fastq_path_1,fastq_path_2"]
for lid, grp, ctrl, n, sig, w, ef in libs:
    r1, r2 = write_fastq(lid, library(n, sig, w, ef))
    rows.append(f"{lid},{grp},1,{ctrl},{r1},{r2}")
with open(os.path.join(OUT, "samplesheet_test.csv"), "w") as fh:
    fh.write("\n".join(rows) + "\n")
print("wrote test data to", OUT)
