# fxstubbe/CUTRUN-nf: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.1.0 - 2026-11-26

### `Added`

- Multi-resolution genomic windows: `--window_sizes` (default `1000,5000,10000,20000,50000,100000`).
  Windows are generated from the FASTA index (`MAKE_WINDOWS`), so no pre-made BED is needed.
- `MERGE_WINDOW_COUNTS` per window size: `04_windows/<bin>/window_counts_<bin>.csv.gz` and
  `library_metrics_<bin>.csv` (spike-in pairs and scale factor per library).
- Spike-in scaled bigWig (`03_Norm/bigwig/<id>.spikein.bigWig`), converted from the spike-in bedGraph.
- SEACR peak calling enabled, paired with the samplesheet `control` (threshold mode without one).
- Samplesheet validation before anything runs (unique, file-name-safe ids, valid controls, FASTQs
  exist). FASTQ paths may be relative to the samplesheet.
- `--fasta` / `--spike_fasta` with optional pre-built indices; Bowtie2 indices are built when missing.
- `-profile test` on a synthetic dataset (`assets/test_data`), plus `conda`, `docker`, `apptainer` profiles.
- MultiQC now also reads TrimGalore, Bowtie2 and samtools stats/flagstat/idxstats.

### `Fixed`

- Window coverage files were named `<id>_[:]_windows.bed`: `ext.prefix` referenced an undefined `splitsize`.
- Process label resources were never applied (`conf/resources.config` nested a `process {}` block
  inside `process {}`): every task ran on 1 CPU locally and requested 12 CPUs / 72 GB / 6 h on Baobab.
  Replaced by `conf/base.config` plus per-profile `resourceLimits`.
- TrimGalore reports were never published (pattern `*.report.txt` vs `*_trimming_report.txt`).
- `normalisation_mode` and `normalisation_to_spikein` had no defaults (bamCoverage received `--normalizeUsing null` without params.yml).
- Boolean switches given on the command line (`--split false`) were ignored with the strict syntax
  parser (Nextflow >= 25.10), which passes them as the string `"false"`.
- MultiQC config path resolved against the launch directory instead of `projectDir`.
- Reference FASTA, index, `.fai` and chromosome sizes are value channels: `CUSTOM_GETCHROMSIZES`
  ran once per library before.
- The mitochondrial `-e` filter moved into `ext.args` so it is parsed as an option on every platform.
- Fragment-length table sorted numerically.

### `Changed`

- Per-library window BEDs moved from `02_alignment/BED/` to `04_windows/<bin>/per_library/`.
- bigWigs moved from `03_Norm/bedgraph/` to `03_Norm/bigwig/`.
- The intermediate `*.filtered.bam` is no longer published by default (`--save_bam_filtered_unsorted`);
  the sorted `*.filter.sort.bam` is.
- Windows exclude the mitochondrial chromosome when `remove_mitochondrial_reads` is set.
- `scale_factors.csv` gains a third column, `spikein_pairs`.

### `Removed`

- `params.splitwindows` / `params.splitsize` (replaced by `window_sizes`), `conf/resources.config`,
  the unused `VALIDATE_METADATA` subworkflow and the empty `workflows/test_validation.nf`.
