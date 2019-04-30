/*!
  * runQualityCardContent component v1.0.0
  * Copyright 2018 IUCT-O
  * Author Frederic Escudie
  * Licensed under GNU General Public License
  */


const rquality_spec_error = {
    'whole_phase_error': { 1.5:"error", 1:"warning", 0.5:"success", 0:"good" },
    'block_error': { 2:"error", 1.5:"warning", 1:"success", 0:"good" },
    'ponctual_error_few_reads': { 0:"good", 1:"warning", 3:"error" },
    'ponctual_error_lot_reads': { 0:"good", 1:"error" }
}


Vue.component('run-quality-card-content', {
    props: {
        'interop_summary': {
            default: null,
            type: Object
        },
        'run_phases': {
            default: null,
            type: Array
        },
        'block_size': {
            default: 20,
            type: Number
        },
        'lot_reads': {
            default: 10,
            type: Number
        },
        'few_reads': {
            default: 5,
            type: Number
        },
        'ponctual_error_threshold': {
            default: 3,
            type: Number
        },
        'alerts_thresholds': {
            default: function(){ return rquality_spec_error },
            type: Object
        }
    },
    computed: {
        thrld_whole_phase_error: function(){
            return new ThresholdAlert(
                this.alerts_thresholds.whole_phase_error,
                "Error rate on whole phase:",
                "",
                "%"
            )
        },
        thrld_block_error: function(){
            return new ThresholdAlert(
                this.alerts_thresholds.block_error,
                "Error rate by block (" + this.block_size + " cycles):",
                "",
                "%"
            )
        },
        thrld_ponctual_error_few_reads: function(){
            return new ThresholdAlert(
                this.alerts_thresholds.ponctual_error_few_reads,
                "",
                "cycles with strong error rate (" + this.ponctual_error_threshold + "%) on " + this.few_reads + "% of reads."
            )
        },
        thrld_ponctual_error_lot_reads: function(){
            return new ThresholdAlert(
                this.alerts_thresholds.ponctual_error_lot_reads,
                "",
                "cycles with strong error rate (" + this.ponctual_error_threshold + "%) on " + this.lot_reads + "% of reads."
            )
        },
        metrics: function(){
            const self = this
            let metrics = {
                'avg_error_in_phase': 0,
                'avg_error_in_block': 0,
                'nb_ponctual_error_few_reads': 0,
                'nb_ponctual_error_lot_reads': 0,
            }
            let average_error_blocks = []
            const few_percentile = (100 - this.lot_reads) + "_percentile"
            const lot_percentile = (100 - this.few_reads) + "_percentile"
            let idx_cycle = 0
            this.run_phases.forEach(function(curr_phase, idx_phase){
                if(!curr_phase.is_index){
                    const nb_cycles_in_phase = curr_phase.nb_cycles
                    let sum_error_upper_quartile = 0
                    let nb_error_upper_quartile = 0
                    let errors_stack = []
                    for(let idx=0 ; idx < nb_cycles_in_phase - 1 ; idx++){
                        // Error on whole read
                        sum_error_upper_quartile += self.interop_summary.error_rate["75_percentile"][idx_cycle]
                        nb_error_upper_quartile++
                        // Error on block
                        errors_stack.push(self.interop_summary.error_rate["75_percentile"][idx_cycle])
                        if(idx + 1 >= self.block_size){
                            let curr_sum = 0
                            for(let j=0 ; j < self.block_size ; j++){
                                curr_sum += errors_stack[j]
                            }
                            average_error_blocks.push(curr_sum/self.block_size)
                            errors_stack.shift()
                        }
                        // Error on single cycle
                        if(self.interop_summary.error_rate[few_percentile][idx] >= self.ponctual_error_thrld){
                            metrics.nb_ponctual_error_lot_reads += 1
                        }
                        if(self.interop_summary.error_rate[lot_percentile][idx] >= self.ponctual_error_thrld){
                            metrics.nb_ponctual_error_few_reads += 1
                        }
                        idx_cycle++
                    }
                    metrics.avg_error_in_phase = Math.max(metrics.avg_error_in_phase, sum_error_upper_quartile/nb_error_upper_quartile).toFixed(2)
                }
            })
            metrics.avg_error_in_block = Math.max(...average_error_blocks).toFixed(2)
            return metrics
        }
    },
    template:
        `<div>
            <div class="row" v-if="interop_summary !== null">
                <div class="col-md-10">
                    <run-error-graph v-if="interop_summary !== null"
                        :run_phases=run_phases
                        :analysis=interop_summary>
                    </run-error-graph>
                </div>
                <div class="col-md-2" if="metrics !== null">
                    <threshold-alert
                        :threshold=thrld_whole_phase_error
                        :value=metrics.avg_error_in_phase>
                    </threshold-alert>
                    <threshold-alert
                        :threshold=thrld_block_error
                        :value=metrics.avg_error_in_block>
                    </threshold-alert>
                    <threshold-alert
                        :threshold=thrld_ponctual_error_few_reads
                        :value=metrics.nb_ponctual_error_few_reads>
                    </threshold-alert>
                    <threshold-alert
                        :threshold=thrld_ponctual_error_lot_reads
                        :value=metrics.nb_ponctual_error_lot_reads>
                    </threshold-alert>
                </div>
            </div>
        </div>`
})
