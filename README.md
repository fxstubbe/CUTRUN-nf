# CUTRUN-nf

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A524.10.5-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with apptainer](https://img.shields.io/badge/run%20with-apptainer-1d355c.svg?labelColor=000000)](https://apptainer.org/)

CUT&RUN analysis with spike-in normalisation on the *E. coli* DNA carried over from the
in-house pA/G-MNase, plus genome-wide bins x libraries count matrices at several resolutions.
Built on the nf-core template (3.3.2) and nf-core/modules.

## Pipeline

1. Adapter and quality trimming (Trim Galore), read QC (FastQC)
2. Alignment to the target genome (Bowtie 2, end-to-end, `--very-sensitive --no-mixed --no-discordant --minins 10 --maxins 700 --dovetail`)
3. Sort, index, stats (samtools); filter MAPQ >= 30, proper pairs, primary alignments, mitochondrial reads removed
4. Alignment of the read pairs that did **not** align to the target genome to the spike-in genome, so no read is counted in both
5. Spike-in scale factor = `normalisation_c` / concordant spike-in pairs (`scale_factors.csv`)
6. Tracks: spike-in scaled fragment bedGraph (bedtools genomecov `-bg -pc -scale`), the same values as bigWig, and a depth-normalised bigWig (deepTools bamCoverage, `--normalizeUsing CPM`)
7. Genomic windows at every size in `--window_sizes` (default 1, 5, 10, 20, 50 and 100 kb): per-library counts and one merged matrix per size
8. Peak calling (SEACR 1.3) on the spike-in scaled bedGraphs, against the control declared in the samplesheet
9. Fragment length distribution, MultiQC report

## Requirements and installation

* Nextflow >= 24.10.5 (tested with 24.10 syntax and with the strict parser of Nextflow 26.10)
* Java 17 or later
* One of Apptainer/Singularity, Docker or conda for the tools (each process has its own container)

```bash
git clone https://github.com/fxstubbe/CUTRUN-nf.git
cd CUTRUN-nf
git checkout v1.1.0          # pin a release for reproducible runs (see CHANGELOG.md)
```

Or let Nextflow fetch a pinned version directly:

```bash
nextflow run fxstubbe/CUTRUN-nf -r v1.1.0 -profile test,apptainer
```

## Quick start

```bash
# end-to-end test on the bundled synthetic dataset (~1 MB, a few minutes)
nextflow run main.nf -profile test,docker            # or test,apptainer / test,conda

# real data
nextflow run main.nf -profile unige_baobab -params-file params.yml
```

`params.yml` holds the run-specific values:

```yaml
metadata:      '/path/to/samplesheet.csv'
outdir:        '/path/to/results'
igenomes_base: '/path/to/igenomes'      # expects the iGenomes layout for WBcel235 and EB1
genome:        'WBcel235'
spike_genome:  'EB1'
window_sizes:  '1000,5000,10000,20000,50000,100000'
```

### On Baobab (UNIGE)

`submit_cr_nextflow.sh` submits the Nextflow head job to SLURM (loads `Java/17` and `Nextflow`,
runs `-profile unige_baobab -params-file ./params.yml` and writes trace, report, timeline and DAG
to `../reports/`). Edit `params.yml` and the `#SBATCH` header, then:

```bash
mkdir -p logs
sbatch submit_cr_nextflow.sh
```

Each process is submitted as its own SLURM job to `shared-cpu`, with the CPU, memory and time
of its label in `conf/base.config`, capped at 12 CPUs / 72 GB / 12 h. Add `-resume` to the
`nextflow run` line to restart a run from the last completed step.

### Other references

Instead of iGenomes, any FASTA can be given with `--fasta` / `--spike_fasta`; without
`--bowtie2_index` / `--spike_bowtie2_index` the index is built by the pipeline.

## Samplesheet

```csv
id,group,replicate,control,fastq_path_1,fastq_path_2
N2_K27me3_1,N2,1,N2_noAb_1,/data/N2_K27me3_1_R1.fastq.gz,/data/N2_K27me3_1_R2.fastq.gz
N2_noAb_1,N2,1,,/data/N2_noAb_1_R1.fastq.gz,/data/N2_noAb_1_R2.fastq.gz
```

* `id` must be unique and contain only letters, digits, `.`, `_` and `-`.
* `control` is optional; when set it must be another `id`. Libraries named as a control are not peak-called themselves.
* FASTQ paths may be absolute, remote (`s3://`, `https://`) or relative to the samplesheet.
* All problems are reported together before anything runs.

## Main parameters

| Parameter | Default | |
|---|---|---|
| `window_sizes` | `1000,5000,10000,20000,50000,100000` | bp; `10kb`/`1mb` suffixes accepted |
| `split` | `true` | build the window matrices |
| `save_window_beds` | `true` | also publish the per-library window BEDs |
| `normalisation_to_spikein` | `true` | `false`: unscaled (`.raw`) bedGraph, no spike-in bigWig |
| `normalisation_mode` | `CPM` | bamCoverage `--normalizeUsing` for the depth-normalised bigWig |
| `normalisation_c` | `10000` | numerator of the spike-in scale factor |
| `minimum_alignment_q_score` | `30` | |
| `remove_mitochondrial_reads` | `true` | also drops the mitochondrion from the windows |
| `run_peak_calling` | `true` | |
| `use_control` | `true` | `false`: SEACR numeric threshold (`seacr_peak_threshold`) for every target |
| `seacr_norm` / `seacr_stringent` | `non` / `stringent` | |

## Outputs

```
<outdir>/
  01_QC/            trimgalore/  fastqc/  multiqc/  Fragments_length_fromBAM/
  02_alignment/     <genome>/{BAM,log,unmapped}   <spike_genome>/{BAM,log}
  03_Norm/          bedgraph/<id>.spikein.bedgraph
                    bigwig/<id>.spikein.bigWig   bigwig/<id>.CPM.bigWig
  04_windows/<bin>/ <genome>.<bin>.windows.bed
                    window_counts_<bin>.csv.gz     chrom,start,end + one column per library
                    library_metrics_<bin>.csv      library_id,group,replicate,control,is_control,
                                                   n_windows,window_fragments,spikein_pairs,spikein_scale_factor
                    per_library/<id>_<bin>_windows.bed   (bedtools coverage columns)
  05_PeakCalling/SEACR/<id>.seacr.peaks.stringent.bed
  scale_factors.csv sample_id,scale_factor,spikein_pairs
```

Notes on the numbers:

* **bedGraph and spike-in bigWig carry identical values**: fragment coverage (each read pair once,
  from the leftmost to the rightmost mate) multiplied by the scale factor, at 1 bp resolution.
  These are what SEACR uses.
* **The CPM bigWig is on a different scale.** bamCoverage counts both mates and every fragment
  overlapping a 50 bp bin, so it should be compared with other CPM tracks only.
* **Window counts are raw** bedtools coverage counts of BAM records: each mate of a pair counts
  once, so a window total is about twice its fragment count, and a read crossing a window edge
  counts in both windows. Summing 1 kb windows therefore gives a few percent more than the
  matching 10 kb window. The `window_fragments` metrics column keeps its historical name.
* Windows follow the FASTA order. The last window of each chromosome is truncated at its end.

## Test data

`-profile test` runs every step on `assets/test_data` (about 1 MB, a few minutes):
a toy "worm" genome (I 120 kb, II 100 kb, X 80 kb, MtDNA 6 kb; a 2 kb block of I duplicated
into II to exercise the MAPQ filter), a 60 kb toy *E. coli* spike-in, and four libraries
of 2 x 50 bp reads: `WT_K27me3_1` and `MUT_K27me3_1`, each with its no-antibody control. The
simulated H3K27me3 domains sit on I:20-60 kb, II:30-50 kb and X:10-70 kb, with the domain on I
at half strength in the mutant. `assets/test_data/make_test_data.py` regenerates the files
byte for byte (fixed seed).

## Troubleshooting

* **A switch on the command line seems ignored** (`--split false`). Fixed in 1.1.0: from Nextflow
  25.10 the strict parser passes command-line booleans as text, which the pipeline now converts.
  Values given in `params.yml` were never affected.
* **`Samplesheet ... failed validation`**: every problem is listed; fix them all and relaunch.
* **No spike-in pairs**: the run continues with a warning and a scale factor of 1 for that
  library, which is then *not* spike-in normalised. Check `02_alignment/<spike_genome>/log/`.
* **Window files from before 1.1.0** were named `<id>_[:]_windows.bed` and lived in
  `02_alignment/BED/`; they now live in `04_windows/<bin>/per_library/`.

## Profiles

| Profile | |
|---|---|
| `test` | synthetic dataset in `assets/test_data` (regenerate with `make_test_data.py`) |
| `unige_baobab` | SLURM `shared-cpu`, Apptainer; per-task resources from the process labels, capped at 12 CPUs / 72 GB / 12 h |
| `my_laptop` | local, conda, capped at 8 CPUs / 12 GB |
| `my_tower` | local, Docker |
| `conda`, `docker`, `apptainer` | software only, to combine with `test` or a custom config |

## Changes

See [`CHANGELOG.md`](CHANGELOG.md). Version 1.1.0 adds multi-resolution windows, SEACR peak
calling, the spike-in bigWig, samplesheet validation and the test profile, and fixes the
process resources, window file names and TrimGalore report publishing.

## Citation

If you use CUTRUN-nf, please cite this repository and the release you ran (`-r v1.1.0`).
A Zenodo DOI will be added here with the first archived release.

## Credits

Written by François-Xavier Stubbe. Uses code and infrastructure from the [nf-core](https://nf-co.re)
community under the [MIT license](https://github.com/nf-core/tools/blob/main/LICENSE):
Ewels P. *et al.* The nf-core framework for community-curated bioinformatics pipelines.
*Nat Biotechnol.* 2020. doi:[10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
Tool references are in [`CITATIONS.md`](CITATIONS.md).
