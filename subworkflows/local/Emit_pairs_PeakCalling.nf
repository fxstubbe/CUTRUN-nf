


workflow SAMPLE_CONTROL_PAIRING{
    take:
    bedgraph_ch  // channel of tuples: [meta, bedgraph]

    main:
    paired_ch = bedgraph_ch
        .toList()  // collect all items so we can build a lookup map
        .map { all_bedgraphs ->

            // Build a lookup map: id -> [meta, bedgraph]
            def bedgraph_map = [:]
            all_bedgraphs.each { meta, bedgraph ->
                bedgraph_map[meta.id] = [meta: meta, bedgraph: bedgraph]
            }

            // Pair samples with their controls
            all_bedgraphs.collect { meta, bedgraph ->
                if (meta.control) {
                    def ctrl_data = bedgraph_map[meta.control]
                    if (ctrl_data) {
                        [
                            meta,
                            bedgraph,
                            ctrl_data.bedgraph,
                            meta.threshold ?: 0.01  // default threshold
                        ]
                    } else {
                        null // control not found, skip
                    }
                } else null
            }.findAll()  // remove nulls
        }
        .flatMap { it }  // flatten the list of pairs into the channel

    emit:
    paired_ch  // channel of tuples: [meta, sample_bedgraph, control_bedgraph, threshold]
}
