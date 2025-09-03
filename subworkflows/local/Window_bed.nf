include { BEDTOOLS_COVERAGE } from '../../modules/nf-core/bedtools/coverage/main'   

workflow SPLIT_WINDOWS {
    take:
    deeptools_ch                // channel of tuples: [meta, bam, bai]
    fasta_fai_channel           // fasta fai channel

    main:



    if(params.split){

    //Formate the chanel to include bam
    bed_bam_ch = deeptools_ch.map{meta, bam, _bai -> [meta, params.splitwindows, bam]}


    //Compute the covarage file per window with bedtools
    BEDTOOLS_COVERAGE(bed_bam_ch, fasta_fai_channel)
    }

    // if(params.split){
    //     // List all BED files in folder
    //     bed_files_ch = Channel.fromPath("${params.splitwindows}/*.bed")

    //     // Format BAM channel
    //     bed_bam_ch = deeptools_ch.map { meta, bam, _bai -> [meta, bam] }

    //     // Cross each BAM with each BED
    //     bam_bed_ready_ch = bed_bam_ch.cross(bed_files_ch)
    //                                 .map { (meta, bam), bed -> [meta, bed, bam] }

    //     // Compute coverage
    //     BEDTOOLS_COVERAGE(bam_bed_ready_ch, fasta_fai_channel)
    // }

   }
    


