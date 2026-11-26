/*
 * Spike-in scale factors, normalised bedGraph and bigWig tracks
 *
 * scale factor = params.normalisation_c / concordantly aligned spike-in pairs, where the
 * spike-in alignment is run on the read pairs that did NOT align to the target genome.
 */

include { DEEPTOOLS_BAMCOVERAGE                                 } from '../../modules/nf-core/deeptools/bamcoverage/main'
include { UCSC_BEDGRAPHTOBIGWIG                                 } from '../../modules/nf-core/ucsc/bedgraphtobigwig/main'
include { BEDTOOLS_GENOMECOV                                    } from '../../modules/nf-core/bedtools/genomecov/main'

// Concordant pairs (exactly 1 time + >1 times) from a paired-end bowtie2 log
def countSpikeinPairs(logFile) {
    long n = 0
    logFile.eachLine { line ->
        def m = line =~ /^\s*(\d+)\s+\([\d.]+%\)\s+aligned concordantly (exactly 1 time|>1 times)/
        if (m.find()) {
            n += m.group(1).toLong()
        }
    }
    return n
}

workflow BAM_NORM {

    take:
    ch_bam_bai      // channel: [ meta, bam, bai ]   filtered, coordinate-sorted
    ch_fasta        // value  : fasta
    ch_fai          // value  : fasta.fai
    ch_spike_log    // channel: [ meta, spike-in bowtie2 log ]
    ch_sizes        // value  : chromosome sizes

    main:

    ch_scale = ch_spike_log.map { meta, logFile ->
        def pairs = countSpikeinPairs(logFile)
        if (pairs == 0) {
            log.warn "[${meta.id}] no concordant spike-in pairs: scale factor set to 1, this library is NOT spike-in normalised"
        }
        def sf = pairs > 0 ? (params.normalisation_c as double) / pairs : 1.0d
        [ meta, pairs, sf ]
    }

    // <outdir>/scale_factors.csv  (first two columns unchanged from earlier versions)
    ch_scale
        .map { meta, pairs, sf -> "${meta.id},${sf},${pairs}" }
        .collectFile(name: 'scale_factors.csv', newLine: true, sort: true,
                     seed: 'sample_id,scale_factor,spikein_pairs', storeDir: params.outdir)

    ch_bam_scale = ch_bam_bai.join(ch_scale)    // [ meta, bam, bai, pairs, sf ]

    // bedGraph (SEACR input): spike-in scaled, or raw when normalisation_to_spikein = false
    ch_genomecov_in = ch_bam_scale.map { meta, bam, _bai, _pairs, sf ->
        params.normalisation_to_spikein.toString().toBoolean() ?
            [ meta + [norm: 'spikein'], bam, sf ] :
            [ meta + [norm: 'raw'],     bam, 1  ]
    }
    BEDTOOLS_GENOMECOV(ch_genomecov_in, ch_sizes, 'bedgraph', true)

    // Depth-normalised bigWig (params.normalisation_mode), always produced
    DEEPTOOLS_BAMCOVERAGE(ch_bam_bai, ch_fasta, ch_fai)

    // Spike-in scaled bigWig: the bedGraph itself converted, so the browser track carries
    // exactly the values SEACR sees (fragment coverage x scale factor, 1 bp resolution).
    // A bamCoverage --scaleFactor track would not: bamCoverage counts both mates and every
    // fragment overlapping a bin, which inflates it ~2 x (L + bin) / L for fragment length L.
    ch_bigwig_spikein = Channel.empty()
    if (params.normalisation_to_spikein.toString().toBoolean()) {
        UCSC_BEDGRAPHTOBIGWIG(BEDTOOLS_GENOMECOV.out.genomecov, ch_sizes)
        ch_bigwig_spikein = UCSC_BEDGRAPHTOBIGWIG.out.bigwig
    }

    emit:
    bedgraph       = BEDTOOLS_GENOMECOV.out.genomecov   // [ meta, bedgraph ]
    scale          = ch_scale                           // [ meta, spike-in pairs, scale factor ]
    bigwig         = DEEPTOOLS_BAMCOVERAGE.out.bigwig
    bigwig_spikein = ch_bigwig_spikein
}
