
include { DEEPTOOLS_BAMCOVERAGE } from '../../modules/nf-core/deeptools/bamcoverage/main'
include { BEDTOOLS_GENOMECOV } from '../../modules/nf-core/bedtools/genomecov/main'


workflow BAM_NORM {

    take:
    deeptools_ch                // channel of tuples: [meta, bam, bai]
    fasta_only_channel_dp       // fasta only channel
    fasta_fai_channel           // fasta fai channel
    spikein_bowtie2_log         // channel of tuples: [meta, log_file]
    sizes_ch                     // channel of chromosome sizes   [meta, sizes]


    main:


    //Make the Coverage file binned with bedtools



    if (params.normalisation_to_spikein) {

        //
        // Compute the Scale factor and join in a channel for bedtools
        // 

        bedtools_ch = deeptools_ch
            .map { meta, bam, bai -> [meta, bam] }  // left: [id, [meta, bam]]
            .join(
                spikein_bowtie2_log.map { meta, log_file ->
                    def aligned_count = log_file.text
                        .split('\n')
                        .findAll { it =~ /aligned concordantly (exactly 1 time|>1 times)/ }
                        .sum { it.split()[0].toInteger() }
                    def scale_factor = params.normalisation_c / (aligned_count != 0 ? aligned_count : params.normalisation_c)
                    [meta, scale_factor]  // right: [id, scale_factor]
                }
            )

        //  
        // Output the scale factors as a CSV file
        // 

        ch_scale_factors_csv = bedtools_ch
            .map { meta, bam, scale_factor -> "${meta.id},${scale_factor}" }
            .toList()
       scale_factors_csv =  WriteScaleFactors(ch_scale_factors_csv)

        //
        // Run BEDTOOLS
        // 

        sizes_only = sizes_ch.map{ meta, sizes -> [sizes]}
        sizes_only.view()
    
        BEDTOOLS_GENOMECOV(
            bedtools_ch,
            sizes_only,
            "bedgraph",
            true
         )
    

        // Collect scale factors into a CSV

    } 
        
    // Call DEEPTOOLS_BAMCOVERAGE without the scale factor
    // This is ALWAYS RUN
    DEEPTOOLS_BAMCOVERAGE(
        deeptools_ch,
        fasta_only_channel_dp,
        fasta_fai_channel
    )
    
}


process WriteScaleFactors {
    tag "scale_factors_csv"

    input:
    val list

    output:
    path "scale_factors.csv"  // File must be created in the work dir

   // publishDir:
   // path "${params.outdir}", mode: 'copy'  // Copy it to the final output dir

    script:
    """
    echo "sample_id,scale_factor" > scale_factors.csv
    printf "%s\\n" ${list.join(" ")} >> scale_factors.csv
    """
}



