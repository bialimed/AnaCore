/*!
  * posCtrlTable component v1.0.0
  * Copyright 2019 IUCT-O
  * Author Frederic Escudie
  * Licensed under GNU General Public License
  */
const pctrl_header = [
    {"title": "Position", "value": function(entry, col){ return entry.chrom + ":" + entry.pos }},
    {"title": "Change", "value": function(entry, col){ return entry.ref + " / " + entry.alt }},
    {"title": "Expected frequency (%)", "value": function(entry, col){ return (entry.expected * 100).toFixed(1) }},
    {"title": "Detected frequency (%)", "value": function(entry, col){ return (entry.detected * 100).toFixed(1) }},
    {
        "title": "Frequency error",
        "value": function(entry, col){
            return (Math.abs(entry.expected - entry.detected) * 100).toFixed(2)
        }
    },
    {
        "title": "Frequency error rate",
        "value": function(entry, col){
            const valid_rate = Math.min(entry.expected, entry.detected) / Math.max(entry.expected, entry.detected)
            return ((1 - valid_rate) * 100).toFixed(2)
        },
        "sort": "desc"
    },
    {
        "title": "Status",
        "value": function(entry, col, max_error_rate){
            let status = 'pass'
            if(entry.detected == 0){
                status = '<span class="pctrltable-label pctrltable-label-danger">lost</span>'
            } else {
                const valid_rate = Math.min(entry.expected, entry.detected) / Math.max(entry.expected, entry.detected)
                if((1 - valid_rate) > max_error_rate){
                    status = '<span class="pctrltable-label pctrltable-label-warning">out of threshold</span>'
                }
            }
            return status
        },
        "is_html": true
    }
]


Vue.component('pos-ctrl-table', {
    extends: DynamicTable,
    props: {
        data: {
            type: Array,
            required: true
        },
        max_error_rate: {
            default: 0.33,
            type: Number
        },
        title: {
            default: "Expected variants",
            type: String
        },
        // Immutable property
        header: {
			default: function(){ return pctrl_header },
			type: Array
		}
    },
    methods: {
		getValue: function( entry, col ){
			return this.$options.filters.getValue( entry, col, this.max_error_rate  )
		},
	},
    filters: {
		getValue: function( entry, col, max_error_rate ){
			let formatted_entry = null
			if( col.hasOwnProperty("value") ){
                if( col["title"] == "Status" ){
                    formatted_entry = col["value"]( entry, col, max_error_rate )
                } else {
				    formatted_entry = col["value"]( entry, col )
                }
			} else if( entry.hasOwnProperty(col['title']) ){
				formatted_entry = entry[col['title']]
			}
			return formatted_entry
		}
    }
})
