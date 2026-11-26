/*
 * Multi-resolution genomic windows
 *
 * For every size in params.window_sizes: tile the genome (bedtools makewindows on the
 * FASTA index), count reads per window per library (bedtools coverage -sorted), and merge
 * the libraries into one bins x libraries matrix plus a per-library metrics table.
 *
 *   <outdir>/04_windows/<bin>/<genome>.<bin>.windows.bed
 *   <outdir>/04_windows/<bin>/window_counts_<bin>.csv.gz
 *   <outdir>/04_windows/<bin>/library_metrics_<bin>.csv
 *   <outdir>/04_windows/<bin>/per_library/<id>_<bin>_windows.bed   (save_window_beds)
 */

include { MAKE_WINDOWS        } from '../../modules/local/make_windows/main'
include { BEDTOOLS_COVERAGE   } from '../../modules/nf-core/bedtools/coverage/main'
include { MERGE_WINDOW_COUNTS } from '../../modules/local/merge_windows/main'

// 1000 -> '1kb', 100000 -> '100kb', 1500 -> '1500bp'
def windowLabel(long size) {
    return size % 1000 == 0 ? "${size.intdiv(1000)}kb".toString() : "${size}bp".toString()
}

// '1000,5000,10kb' (string), 10000 (number) or [1000, 5000] (list) -> sorted unique sizes in bp
def parseWindowSizes(spec) {
    def tokens = spec instanceof List ? spec.collect { it.toString() } : spec.toString().tokenize(',')
    def mult = [bp: 1L, kb: 1000L, mb: 1000000L]
    def sizes = tokens.collect { it.trim().toLowerCase() }.findAll { it }.collect { tok ->
        def m = tok =~ /^(\d+)\s*(bp|kb|mb)?$/
        if (!m.matches()) {
            error "Invalid entry '${tok}' in --window_sizes (expected e.g. 1000,5000,10kb)"
        }
        long n = m.group(1).toLong() * mult[m.group(2) ?: 'bp']
        if (n <= 0) {
            error "Window size must be positive, got '${tok}'"
        }
        n
    }
    if (!sizes) {
        error "--window_sizes is empty"
    }
    return sizes.unique().sort()
}

workflow SPLIT_WINDOWS {

    take:
    ch_bam_bai      // channel: [ meta, bam, bai ]
    ch_fai          // value  : fasta.fai
    ch_spike_log    // channel: [ meta, spike-in bowtie2 log ]
    ch_metadata     // value  : samplesheet
    exclude_chroms  // string : chromosomes left out of the windows (e.g. 'MtDNA')

    main:

    ch_window_def = Channel.fromList(
        parseWindowSizes(params.window_sizes).collect { s -> [bin: windowLabel(s), size: s] }
    )

    MAKE_WINDOWS(ch_window_def.combine(ch_fai), exclude_chroms)

    // every library x every window size
    ch_coverage_in = ch_bam_bai
        .combine(MAKE_WINDOWS.out.bed)
        .map { meta, bam, _bai, wmeta, bed -> [ meta + [bin: wmeta.bin, bin_size: wmeta.size], bed, bam ] }

    BEDTOOLS_COVERAGE(ch_coverage_in, ch_fai)

    ch_by_bin = BEDTOOLS_COVERAGE.out.bed
        .map { meta, bed -> [ meta.bin, bed ] }
        .groupTuple()

    MERGE_WINDOW_COUNTS(
        ch_by_bin,
        ch_spike_log.map { _meta, logFile -> logFile }.collect(),
        ch_metadata
    )

    emit:
    windows  = MAKE_WINDOWS.out.bed              // [ wmeta, windows.bed ]
    coverage = BEDTOOLS_COVERAGE.out.bed         // [ meta (+bin), coverage.bed ]
    counts   = MERGE_WINDOW_COUNTS.out.counts    // [ bin, window_counts_<bin>.csv.gz ]
    metrics  = MERGE_WINDOW_COUNTS.out.metrics   // [ bin, library_metrics_<bin>.csv ]
}
