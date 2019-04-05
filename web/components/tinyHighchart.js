/*!
  * tinyHighchart component v1.1.0
  * Copyright 2018 IUCT-O
  * Author Frederic Escudie
  * Licensed under GNU General Public License
  */


function getGraph( pChart_type, pSeries, pX_title, pY_title, pTitle=null, pTickInterval=1 ){
	let chart_param = {
        chart: {
            type: pChart_type,
            zoomType: "x",
        },
        xAxis: {
            tickInterval: pTickInterval,
            title: {
                text: pX_title
            }
        },
        yAxis: {
            title: {
                text: pY_title
            }
        },
        tooltip: {
            crosshairs: [true],
            shared: true,
			formatter: function(){
				const self = this
				let lines = [this.points[0].series.xAxis.axisTitle.textStr + ' <b>' + this.points[0].x + '</b>:']
                this.points.forEach(function(point, idx) {
                    lines.push(
                        '<b><span style="color:' + point.series.color + '">' + point.series.name + '</b></span>' +
                        ': ' +
                        '<b>' + point.y + '</b>'
                    )
                });
                return lines.join('<br />')
            }
        },
        credits: {
            enabled: false
        },
        series: pSeries
    }
    if( pTitle != null ){
		chart_param["title"] = { "text": pTitle }
	}
    return chart_param
}


Vue.component('tiny-highchart', {
    props: {
		series: {
            type: Array,
            required: true
        },
		chart_type: {
			type: String,
            default: "line"
		},
        x_title: {
            type: String,
            required: true
        },
        y_title: {
            type: String,
            required: true
        },
		title: {
            type: String,
            default: null
        },
        tick_interval: {
            type: Number,
            default: null
        }
	},
    mounted: function() {
		this.drawGraph()
	},
    watch: {
        count_by_cat: function(val){
            this.drawGraph()
        }
    },
    computed: {
        graph_parameters: function(){
            return getGraph( this.chart_type, this.series, this.x_title, this.y_title, this.title, this.tick_interval )
        }
    },
    methods: {
        drawGraph: function(){
		    const ctx = this.$el
            const scroll_pos = $(window).scrollTop()
        	Highcharts.chart( ctx, this.graph_parameters )
            $(window).scrollTop(scroll_pos)
        }
    },
	template: `<div></div>`
})
