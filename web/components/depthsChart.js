/*!
  * depthsChart component v1.0.0
  * Copyright 2019 IUCT-O
  * Author Frederic Escudie
  * Licensed under GNU General Public License
  */


Vue.component('depths-chart', {
    props: {
		analysis: {
            type: Object,
            required: true
        },
	},
    computed: {
        chart_series: function(){
            const self = this
            let chart_series = []
            Object.keys(self.analysis.depths_classes.count_by_spl).forEach(function( curr_spl ){
                let curr_series = {
                    "data": [],
                    "name": curr_spl
                }
                const counts = self.analysis.depths_classes.count_by_spl[curr_spl]
                self.analysis.depths_classes.depths_list.forEach(function(curr_depth, curr_idx){
                    curr_series["data"].push([
                        curr_depth,
                        counts[curr_idx]
                    ])
                })
                chart_series.push(curr_series)
            })
            return chart_series
        },
    },
	template:
        `<tiny-highchart
            chart_type="area"
            :series=chart_series
            x_title="Depth"
            y_title="Number of nucleotids with this depth"
            title="Depths distribution"
            :tick_interval=1>
        </tiny-highchart>`
})
