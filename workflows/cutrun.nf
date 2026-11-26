/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CUTRUN-nf main workflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    TrimGalore -> FastQC -> Bowtie2 (target) -> sort -> MAPQ/proper-pair/mito filter
    Bowtie2 (spike-in, on the pairs that did not align to the target)
    -> spike-in scale factor -> bedGraph + bigWigs -> multi-resolution windows
    -> SEACR peaks -> fragment length -> MultiQC
----------------------------------------------------------------------------------------
*/

include { TRIMGALORE                             } from '../modules/nf-core/trimgalore/main'
include { FASTQC                                 } from '../modules/nf-core/fastqc/main'
include { MULTIQC                                } from '../modules/nf-core/multiqc/main'
include { BOWTIE2_BUILD as BOWTIE2_BUILD_TARGET  } from '../modules/nf-core/bowtie2/build/main'
include { BOWTIE2_BUILD as BOWTIE2_BUILD_SPIKEIN } from '../modules/nf-core/bowtie2/build/main'
include { BOWTIE2_ALIGN as BOWTIE2_TARGET_ALIGN  } from '../modules/nf-core/bowtie2/align/main'
include { BOWTIE2_ALIGN as BOWTIE2_SPIKEIN_ALIGN } from '../modules/nf-core/bowtie2/align/main'
include { SAMTOOLS_VIEW                          } from '../modules/nf-core/samtools/view/main'
include { CUSTOM_GETCHROMSIZES                   } from '../modules/nf-core/custom/getchromsizes/main'
include { SEACR_CALLPEAK                         } from '../modules/nf-core/seacr/callpeak/main'
include { FRAGMENT_LEN                           } from '../modules/local/Fragment_len/main'

include { BAM_SORT_STATS_SAMTOOLS                         } from '../subworkflows/nf-core/bam_sort_stats_samtools/main.nf'
include { BAM_SORT_STATS_SAMTOOLS as SAMTOOLS_VIEW_SORT   } from '../subworkflows/nf-core/bam_sort_stats_samtools/main.nf'
include { BAM_NORM                                        } from '../subworkflows/local/Bam_Norm.nf'
include { SPLIT_WINDOWS                                   } from '../subworkflows/local/Window_bed.nf'
include { SAMPLE_CONTROL_PAIRING                          } from '../subworkflows/local/Emit_pairs_PeakCalling.nf'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    HELPERS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def genomeEntry(name) {
    return (params.genomes && name && params.genomes.containsKey(name)) ? params.genomes[name] : [:]
}

// Absolute or remote (s3://, https://, ...) paths are used as is; relative paths are
// resolved against the directory that holds the samplesheet.
def resolveFastq(sheet, String p) {
    def remote = p =~ /^[A-Za-z][A-Za-z0-9+.-]*:\/\//
    return (p.startsWith('/') || remote.find()) ? file(p) : sheet.parent.resolve(p)
}

/*
 * Read and validate the samplesheet: id,group,replicate,control,fastq_path_1,fastq_path_2
 *   - ids unique and file-name safe; every non-empty control is an existing id
 *   - FASTQ paths absolute, remote (s3://, https://) or relative to the samplesheet
 *   - every FASTQ exists
 * All problems are reported together.
 */
def readSamplesheet(sheetPath) {
    def sheet    = file(sheetPath, checkIfExists: true)
    def rows     = sheet.splitCsv(header: true, strip: true)
    def required = ['id', 'group', 'replicate', 'fastq_path_1', 'fastq_path_2']
    def problems = []

    if (!rows) {
        error "Samplesheet ${sheet} has no rows"
    }
    def missingCols = required.findAll { col -> !rows[0].containsKey(col) }
    if (missingCols) {
        error "Samplesheet ${sheet} is missing column(s): ${missingCols.join(', ')}"
    }

    def ids = rows.collect { row -> row.id }
    ids.countBy { it }.findAll { _id, n -> n > 1 }.each { id, n -> problems << "id '${id}' appears ${n} times" }
    ids.findAll { id -> !(id ==~ /[A-Za-z0-9._-]+/) }.each { id -> problems << "id '${id}' may only contain letters, digits, '.', '_' and '-'" }

    def controlIds = rows.collect { row -> row.control ?: '' }.findAll { c -> c }.toSet()
    controlIds.findAll { c -> !(c in ids) }.each { c -> problems << "control '${c}' is not an id in the samplesheet" }

    def samples = rows.collect { row ->
        def reads = ['fastq_path_1', 'fastq_path_2'].collect { col ->
            def p = row[col]
            if (!p) {
                problems << "${row.id}: ${col} is empty"
                return null
            }
            def f = resolveFastq(sheet, p)
            if (!f.exists()) {
                problems << "${row.id}: ${col} not found -> ${f}"
            }
            f
        }
        def meta = [
            id        : row.id,
            group     : row.group,
            replicate : row.replicate,
            control   : row.control ?: '',
            is_control: row.id in controlIds,
            single_end: false
        ]
        [ meta, reads ]
    }

    if (problems) {
        error "Samplesheet ${sheet} failed validation:\n  - " + problems.join('\n  - ')
    }
    log.info "Samplesheet OK: ${samples.size()} libraries, ${controlIds.size()} used as control"
    return samples
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow CutRun {

    main:

    // ---- Parameters that must be set
    if (!params.metadata) {
        error "Please provide a samplesheet with --metadata (or in params.yml)"
    }
    if (!params.outdir) {
        error "Please provide an output directory with --outdir (or in params.yml)"
    }

    // ---- References: explicit FASTA/index first, iGenomes entry otherwise
    def target       = genomeEntry(params.genome)
    def spike        = genomeEntry(params.spike_genome)
    def target_fasta = params.fasta ?: target.fasta
    def spike_fasta  = params.spike_fasta ?: spike.fasta
    def target_index = params.bowtie2_index ?: (params.fasta ? null : target.bowtie2)
    def spike_index  = params.spike_bowtie2_index ?: (params.spike_fasta ? null : spike.bowtie2)
    if (!target_fasta) {
        error "No FASTA for target genome '${params.genome}': set --fasta or a known --genome"
    }
    if (!spike_fasta) {
        error "No FASTA for spike-in genome '${params.spike_genome}': set --spike_fasta or a known --spike_genome"
    }
    def mito           = params.mito_name ?: target.mito_name
    def exclude_chroms = (params.remove_mitochondrial_reads.toString().toBoolean() && mito) ? mito.toString() : ''

    ch_fasta       = Channel.value([ [id: params.genome],       file(target_fasta, checkIfExists: true) ])
    ch_spike_fasta = Channel.value([ [id: params.spike_genome], file(spike_fasta,  checkIfExists: true) ])

    ch_index = Channel.empty()
    if (target_index) {
        ch_index = Channel.value([ [id: params.genome], file(target_index, checkIfExists: true) ])
    } else {
        BOWTIE2_BUILD_TARGET(ch_fasta)
        ch_index = BOWTIE2_BUILD_TARGET.out.index
    }
    ch_spike_index = Channel.empty()
    if (spike_index) {
        ch_spike_index = Channel.value([ [id: params.spike_genome], file(spike_index, checkIfExists: true) ])
    } else {
        BOWTIE2_BUILD_SPIKEIN(ch_spike_fasta)
        ch_spike_index = BOWTIE2_BUILD_SPIKEIN.out.index
    }

    // FASTA index and chromosome sizes, computed once
    CUSTOM_GETCHROMSIZES(ch_fasta)
    ch_fai        = CUSTOM_GETCHROMSIZES.out.fai.map   { _meta, fai   -> fai   }
    ch_sizes      = CUSTOM_GETCHROMSIZES.out.sizes.map { _meta, sizes -> sizes }
    ch_fasta_path = ch_fasta.map { _meta, fasta -> fasta }

    // ---- Samplesheet
    ch_metadata = Channel.value(file(params.metadata, checkIfExists: true))
    ch_reads    = Channel.fromList(readSamplesheet(params.metadata))

    // ####################################################################################
    // Trimming and read QC
    // ####################################################################################

    TRIMGALORE(ch_reads)
    FASTQC(TRIMGALORE.out.reads)

    // ####################################################################################
    // Alignment to the target genome, sort, filter
    // ####################################################################################

    BOWTIE2_TARGET_ALIGN(
        TRIMGALORE.out.reads,
        ch_index,
        ch_fasta,
        true,   // save_unaligned: these pairs feed the spike-in alignment
        false   // sort_bam: sorted below
    )

    BAM_SORT_STATS_SAMTOOLS(BOWTIE2_TARGET_ALIGN.out.bam, ch_fasta)

    SAMTOOLS_VIEW(
        BAM_SORT_STATS_SAMTOOLS.out.bam.join(BAM_SORT_STATS_SAMTOOLS.out.bai),
        ch_fasta,
        [],
        []
    )

    SAMTOOLS_VIEW_SORT(SAMTOOLS_VIEW.out.bam, ch_fasta)

    ch_bam_bai = SAMTOOLS_VIEW_SORT.out.bam.join(SAMTOOLS_VIEW_SORT.out.bai)   // [ meta, bam, bai ]

    // ####################################################################################
    // Spike-in alignment on the pairs that did not align to the target genome
    // ####################################################################################

    BOWTIE2_SPIKEIN_ALIGN(
        BOWTIE2_TARGET_ALIGN.out.fastq,
        ch_spike_index,
        ch_spike_fasta,
        false,  // save_unaligned
        true    // sort_bam
    )

    // ####################################################################################
    // Scale factors, bedGraph, bigWig
    // ####################################################################################

    BAM_NORM(
        ch_bam_bai,
        ch_fasta_path,
        ch_fai,
        BOWTIE2_SPIKEIN_ALIGN.out.log,
        ch_sizes
    )

    // ####################################################################################
    // Genomic windows at every size in params.window_sizes
    // ####################################################################################

    if (params.split.toString().toBoolean()) {
        SPLIT_WINDOWS(
            ch_bam_bai,
            ch_fai,
            BOWTIE2_SPIKEIN_ALIGN.out.log,
            ch_metadata,
            exclude_chroms
        )
    }

    // ####################################################################################
    // Fragment length distribution
    // ####################################################################################

    if (params.fragment_len_reporting.toString().toBoolean()) {
        FRAGMENT_LEN(ch_bam_bai)
    }

    // ####################################################################################
    // Peak calling: SEACR on the (spike-in scaled) bedGraphs, against the declared control
    // ####################################################################################

    if (params.run_peak_calling.toString().toBoolean()) {
        SAMPLE_CONTROL_PAIRING(BAM_NORM.out.bedgraph)
        SEACR_CALLPEAK(SAMPLE_CONTROL_PAIRING.out.paired_ch, params.seacr_peak_threshold)
    }

    // ####################################################################################
    // MultiQC
    // ####################################################################################

    ch_multiqc_files = Channel.empty()
        .mix(FASTQC.out.zip.map              { _meta, f -> f })
        .mix(TRIMGALORE.out.log.map          { _meta, f -> f })
        .mix(BOWTIE2_TARGET_ALIGN.out.log.map { _meta, f -> f })
        .mix(SAMTOOLS_VIEW_SORT.out.stats.map    { _meta, f -> f })
        .mix(SAMTOOLS_VIEW_SORT.out.flagstat.map { _meta, f -> f })
        .mix(SAMTOOLS_VIEW_SORT.out.idxstats.map { _meta, f -> f })

    MULTIQC(
        ch_multiqc_files.collect(),
        Channel.value(file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true)),
        [],
        [],
        [],
        []
    )
}
