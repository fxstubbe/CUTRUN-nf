/*
 * Pair each target bedGraph with the bedGraph of its control (samplesheet `control` column)
 *
 * emits [ meta, target.bedgraph, control.bedgraph ]  when a control is declared and use_control
 *       [ meta, target.bedgraph, [] ]                otherwise (SEACR numeric-threshold mode)
 * Libraries that serve as somebody's control are not peak-called themselves.
 */

workflow SAMPLE_CONTROL_PAIRING {

    take:
    bedgraph_ch     // channel: [ meta, bedgraph ]

    main:

    ch_targets  = bedgraph_ch.filter { meta, _bg -> !meta.is_control }
    ch_controls = bedgraph_ch
        .filter { meta, _bg -> meta.is_control }
        .map { meta, bg -> [ meta.id, bg ] }

    paired_ch = Channel.empty()
    if (params.use_control.toString().toBoolean()) {
        ch_with_control = ch_targets
            .filter { meta, _bg -> meta.control }
            .map { meta, bg -> [ meta.control, meta, bg ] }
            .combine(ch_controls, by: 0)
            .map { _cid, meta, bg, ctrl -> [ meta, bg, ctrl ] }
        ch_no_control = ch_targets
            .filter { meta, _bg -> !meta.control }
            .map { meta, bg -> [ meta, bg, [] ] }
        paired_ch = ch_with_control.mix(ch_no_control)
    } else {
        paired_ch = ch_targets.map { meta, bg -> [ meta, bg, [] ] }
    }

    emit:
    paired_ch       // channel: [ meta, bedgraph, control bedgraph or [] ]
}
