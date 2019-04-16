/*!
  * runDensityGraph component v1.0.0
  * Copyright 2018 IUCT-O
  * Author Frederic Escudie
  * Licensed under GNU General Public License
  */


Vue.component('run-density-graph', {
    props: {
		analysis: {
            type: Object,
            required: true
        },
	},
    data: function(){
        return {"categories": ["Run"]}
    },
    computed: {
        chart_series: function(){
            return [
                {
                    "name": "Density",
                    "data": [[
                        this.analysis.density.min,
                        this.analysis.density["25_percentile"],
                        this.analysis.density["50_percentile"],
                        this.analysis.density["75_percentile"],
                        this.analysis.density.max
                    ]]
                }, {
                    "name": "Density PF",
                    "data": [[
                        this.analysis.density_PF.min,
                        this.analysis.density_PF["25_percentile"],
                        this.analysis.density_PF["50_percentile"],
                        this.analysis.density_PF["75_percentile"],
                        this.analysis.density_PF.max
                    ]]
                }
            ]
        }
    },
    template:
        `<tiny-highchart
            chart_type="boxplot"
            :series="chart_series"
            title="FlowCell clusters density"
            y_title="Clusters/mm2"
            :x_tick_categories=categories>
        </tiny-highchart>`
})
