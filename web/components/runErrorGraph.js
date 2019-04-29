/*!
  * runErrorGraph component v1.0.0
  * Copyright 2018 IUCT-O
  * Author Frederic Escudie
  * Licensed under GNU General Public License
  */


Vue.component('run-error-graph', {
    props: {
        'title': {
            type: String,
            default: "Control sequence error rate"
        },
		'analysis': {
            type: Object,
            required: true
        },
        'run_phases': {
            type: Array,
            required: true
        }
	},
    computed: {
        categories: function(){
            let x_ticks = []
            let x_max = 0
            this.chart_series.forEach(function(curr_series){
                x_max = Math.max(x_max, curr_series.data.length)
            })
            return Array.from({length: x_max}, (val, key) => key+1)
        },
        chart_series: function(){
            const self = this
            let series = []
            this.run_phases.forEach(function(curr_phase){
                if(curr_phase.title.startsWith("R")){
                    series.push({
                        "name": curr_phase.title,
                        "data": []
                    })
                }
            })
            let idx_cycle = 0
            let idx_category = 0
            this.run_phases.forEach(function(curr_phase, idx_phase){
                if(curr_phase.title.startsWith("R")){
                    const nb_cycles_in_phase = self.run_phases[idx_phase].nb_cycles
                    for(let idx=0 ; idx < nb_cycles_in_phase - 1 ; idx++){
                        series[idx_category].data.push([
                            self.analysis.error_rate.min[idx_cycle],
                            self.analysis.error_rate["25_percentile"][idx_cycle],
                            self.analysis.error_rate["50_percentile"][idx_cycle],
                            self.analysis.error_rate["75_percentile"][idx_cycle],
                            self.analysis.error_rate.max[idx_cycle]
                        ])
                        idx_cycle++
                    }
                    idx_category++
                }
            })
            return series
        }
    },
    template:
        `<tiny-highchart
            chart_type="boxplot"
            :series="chart_series"
            :title="title"
            y_title="Error rate (%)"
            :y_min=0
            x_title="Sequencing cycle"
            :x_tick_categories=categories>
        </tiny-highchart>`
})
