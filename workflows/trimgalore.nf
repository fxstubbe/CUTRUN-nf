


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


// Imports the TRIMGALORE process
include { TRIMGALORE } from '../modules/nf-core/trimgalore/main'
include { FASTQC } from '../modules/nf-core/fastqc/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Define the workflow
workflow CleanReads {


            input_Reads_Channel = Channel
            .fromPath(params.metadata)
            .splitCsv(header: true)
            .map { row ->
                def meta = [id: row.id] // Assuming one of the columns is named 'id'
                def reads = [file(row.fastq_path_1), file(row.fastq_path_2)] // Assuming columns are named 'fastq_path_1' and 'fastq_path_2'
                [meta, reads]
            }

   // take:
    //Channel of the medata file
   // Reads_ch

    //main:
    // Print the input data
   // Reads_ch.view()

    // Pass the input data to TRIMGALORE
    TRIMGALORE(input_Reads_Channel)

    // Print the output
    TRIMGALORE.out.reads.view()


    // Connect the output of TRIMGALORE to SECOND_PROCESS
    //SECOND_PROCESS(TRIM.out.reads)

    // Print the output of SECOND_PROCESS
    //SECOND_PROCESS.out.processed_reads.view()

}