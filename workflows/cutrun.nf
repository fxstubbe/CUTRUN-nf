


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Imports the MODULES process
include { FASTQC } from '../modules/nf-core/fastqc/main'
include { MULTIQC } from '../modules/nf-core/multiqc/main'
include { TRIMGALORE } from '../modules/nf-core/trimgalore/main'
include { BOWTIE2_ALIGN as BOWTIE2_TARGET_ALIGN } from '../modules/nf-core/bowtie2/align/main'
include { BOWTIE2_ALIGN as BOWTIE2_SPIKEIN_ALIGN } from '../modules/nf-core/bowtie2/align/main'
include { SAMTOOLS_SORT } from '../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_VIEW } from '../modules/nf-core/samtools/view/main'
include { SAMTOOLS_FLAGSTAT } from '../modules/nf-core/samtools/flagstat/main'
include { SAMTOOLS_IDXSTATS } from '../modules/nf-core/samtools/idxstats/main'
include { PICARD_MARKDUPLICATES   } from '../modules/nf-core/picard/markduplicates/main'
include { PICARD_ADDORREPLACEREADGROUPS } from '../modules/nf-core/picard/addorreplacereadgroups/main'
include { CUSTOM_GETCHROMSIZES } from '../modules/nf-core/custom/getchromsizes/main'
include { DEEPTOOLS_BAMCOVERAGE } from '../modules/nf-core/deeptools/bamcoverage/main'
include { SEACR_CALLPEAK } from '../modules/nf-core/seacr/callpeak/main'

// Homemade module
include{FRAGMENT_LEN} from '../modules/local/Fragment_len/main'

//Import SUBWORFLOWS
include { BAM_SORT_STATS_SAMTOOLS   } from '../subworkflows/nf-core/bam_sort_stats_samtools/main.nf'
include { BAM_SORT_STATS_SAMTOOLS as SAMTOOLS_VIEW_SORT } from '../subworkflows/nf-core/bam_sort_stats_samtools/main.nf'
include { SAMPLE_CONTROL_PAIRING } from '../subworkflows/local/Emit_pairs_PeakCalling.nf'
include {VALIDATE_METADATA} from '../subworkflows/local/Validate_metadata.nf'   
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
//
//



// Define the workflow
workflow CutRun {

    // ####################################################################################
    // Validate metadata from CSV file
    // ####################################################################################
    
    // VALIDATE_METADATA(params.metadata)


    // ####################################################################################
    // INPUTS :  
    // Takes a .csv file [id, group, replicate, fastq_path_1, fastq_path_2]
    // Generates a Groovy list [id, [fastq_path_1, fastq_path_2]]
    // ####################################################################################

    // Prepare the input channel from the metadata CSV file
    input_Reads_Channel = Channel
                            .fromPath(params.metadata)
                            .splitCsv(header: true)
                            .map { row ->
                                def meta = [
                                    id: row.id,
                                    group: row.group,
                                    replicate: row.replicate,
                                    control: row.control
                                ]
                                def reads = [file(row.fastq_path_1), file(row.fastq_path_2)]
                                [meta, reads]
                            }
                                            

    // ####################################################################################
    // TRIM-GALORE : Adapter and Quality trimming
    // ####################################################################################

    TRIMGALORE(input_Reads_Channel) // for debug TRIMGALORE.out.reads.view()

    // ####################################################################################
    // FASTQC : QC reports 
    // Compile into a MULTI_QC report
    // ####################################################################################

    ch_multiqc_files = Channel.empty()

    FASTQC(TRIMGALORE.out.reads)

    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect { _sample_id, zip_file -> zip_file }) //Get path to ZIP file
    ch_multiqc_config = Channel.fromPath("./assets/multiqc_config.yml", checkIfExists: true) //Config file for fastqc (OPTIONAL)

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        [],
        [],
        [],
        []
    )

    // ####################################################################################
    // BOWTIE2 : Alignment 
    // ####################################################################################

    // 1) Run bowtie2 on target genome
    // ----------------------------------------------------

    // Prepare channels for BOWTIE2 process
    fasta_channel = input_Reads_Channel.map { meta, _reads -> [meta, file(params.genomes[params.genome].fasta)] }
    index_channel = input_Reads_Channel.map { meta, _reads -> [meta, file(params.genomes[params.genome].bowtie2)] } 


    // Run BOWTIE2_TARGET_ALIGN process
   BOWTIE2_TARGET_ALIGN (
       TRIMGALORE.out.reads, // Paths of trimmed reads generated by trim-galore
       index_channel, // Channel poiting to the bowtie2 index
       fasta_channel, // Channel poiting to the FASTA file 
       true,  // save_unaligned; can be set to save_align when figure out conditional mapping
       false   // sort_basort BAM, False because already included in BAM_SORT_STATS_SAMTOOLS subworkflow
   )
    // Sort output bam
    BAM_SORT_STATS_SAMTOOLS( BOWTIE2_TARGET_ALIGN.out.bam, fasta_channel)

    // Filter based on q-score
    samtool_ch = BAM_SORT_STATS_SAMTOOLS.out.bam
                                        .join(BAM_SORT_STATS_SAMTOOLS.out.bai)  // join on meta
                                        .map { meta, bam, bai ->
                                            [meta, bam, bai ]   // bai replaces your []
                                        }
    SAMTOOLS_VIEW(samtool_ch, fasta_channel, [], [])

    // Sort output bam
    SAMTOOLS_VIEW_SORT(SAMTOOLS_VIEW.out.bam, fasta_channel)


    // ####################################################################################
    //(OPTIONAL)  BOWTIE2 : Alignment on spike-in genome for normalization
    // ####################################################################################

    // // Prepare channels for BOWTIE2 process
    spike_fasta_channel = input_Reads_Channel.map { meta, _reads -> [meta, file(params.genomes[params.spike_genome].fasta)] }
    spike_index_channel = input_Reads_Channel.map { meta, _reads -> [meta, file(params.genomes[params.spike_genome].bowtie2)] } 

    // Run BOWTIE2_TARGET_ALIGN process
    BOWTIE2_SPIKEIN_ALIGN (
        BOWTIE2_TARGET_ALIGN.out.fastq, // Paths of trimmed reads generated by trim-galore
        spike_index_channel, // Channel poiting to the bowtie2 index
        spike_fasta_channel, // Channel poiting to the FASTA file 
        true,  // save_unaligned
        true   // sort_bam
    )
   


    // ####################################################################################
    // Deeptools: Read Normalization
    // ####################################################################################

    // Prepare channels for DEEPTOOLS_BAMCOVERAGE process
    deeptools_ch = SAMTOOLS_VIEW_SORT.out.bam
                                    .join(SAMTOOLS_VIEW_SORT.out.bai)  // join on meta
                                    .map { meta, bam, bai ->
                                        tuple (meta, bam, bai)    // bai replaces your []
                                    }
     // Prepare fasta channel for DEEPTOOLS_BAMCOVERAGE                                
    fasta_only_channel_dp = fasta_channel.map { meta, fasta -> fasta }

    // Prepare fasta fai channel for DEEPTOOLS_BAMCOVERAGE
    fasta_fai_channel = fasta_only_channel_dp.map { fasta -> file("${fasta}.fai") }

    // Run DEEPTOOLS_BAMCOVERAGE process
    // First Call to get bedgraph for SEACR
    DEEPTOOLS_BAMCOVERAGE(
        deeptools_ch,
        fasta_only_channel_dp,
        fasta_fai_channel
    )
    DEEPTOOLS_BAMCOVERAGE.out.bedgraph.view()

    // Prepare Chromosome size file for UCSC BEDGRAPHTOBIGWIG
    // CUSTOM_GETCHROMSIZES(fasta_channel)

    // Convert bedgraph to bigwig using UCSC BEDGRAPHTOBIGWIG
    // UCSC_BEDGRAPHTOBIGWIG(

    //     DEEPTOOLS_BAMCOVERAGE.out.bedgraph, // Input bedgraph channel
    //     CUSTOM_GETCHROMSIZES.out.sizes.map { meta, sizes -> sizes }

    // )

    // ####################################################################################
    // Extract fragment length for QC
    // ####################################################################################

    // Get fragment length
    FRAGMENT_LEN(deeptools_ch)

    // ####################################################################################
    // PEAK CALLING
    // ####################################################################################

    // Emit pairs for peak calling using SEACR
    SAMPLE_CONTROL_PAIRING(DEEPTOOLS_BAMCOVERAGE.out.bedgraph)
    SAMPLE_CONTROL_PAIRING.out.paired_ch.view() //Sanity checl

    // Run SEACR peak calling
    SEACR_CALLPEAK(SAMPLE_CONTROL_PAIRING.out.paired_ch,1)

}















//-------------------------------------------------------------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------------------------------------------------
//
// Bits of code that were not used in the final workflow but might be useful for future reference or debugging
//





    // input_Reads_Channel = Channel
    //                         .fromPath(params.metadata)
    //                         .splitCsv(header: true)
    //                         .map { row ->
    //                             def meta = [id: row.id] // Assuming one of the columns is named 'id'
    //                             def reads = [file(row.fastq_path_1), file(row.fastq_path_2)] // Assuming columns are named 'fastq_path_1' and 'fastq_path_2'
    //                             [meta, reads]
    //                         }

    
    // ####################################################################################
    // PICARD : Remove duplicates
    // ####################################################################################

    //Prepare fai fasta channel for picard
    // fasta_fai_channel = fasta_channel.map { _meta, fasta -> [[id : params.genome], file("${fasta}.fai")] }

    
    // PICARD_ADDORREPLACEREADGROUPS(        
    //     SAMTOOLS_SORT.out.bam, //BAM file to feed picard
    //     fasta_channel, // Reference genome fasta
    //     fasta_fai_channel)

    // PICARD_MARKDUPLICATES (
    //     PICARD_ADDORREPLACEREADGROUPS.out.bam, //BAM file to feed picard
    //     fasta_channel, // Reference genome fasta
    //     fasta_fai_channel
    // )
// BOWTIE2_TARGET_ALIGN.out.

//  ch_markduplicates_metrics = Channel.empty()
//     if (params.run_mark_dups) {
//         PICARD_MARKDUPLICATES (
//             BAM_SORT_STATS_SAMTOOLS.out.bam, //BAM file to feed picard
//             BAM_SORT_STATS_SAMTOOLS.out.bai,
//             true,
//             PREPARE_GENOME.out.fasta.collect(),
//             PREPARE_GENOME.out.fasta_index.collect()
//         )


//         ch_samtools_bam           = MARK_DUPLICATES_PICARD.out.bam
//         ch_samtools_bai           = MARK_DUPLICATES_PICARD.out.bai
//         ch_samtools_stats         = MARK_DUPLICATES_PICARD.out.stats
//         ch_samtools_flagstat      = MARK_DUPLICATES_PICARD.out.flagstat
//         ch_samtools_idxstats      = MARK_DUPLICATES_PICARD.out.idxstats
//         ch_markduplicates_metrics = MARK_DUPLICATES_PICARD.out.metrics
//         ch_software_versions      = ch_software_versions.mix(MARK_DUPLICATES_PICARD.out.versions)
//     }




