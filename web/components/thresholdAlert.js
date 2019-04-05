class ThresholdAlert {
    constructor(rules, text_pre, text_post="", value_suffix="") {
        this.text_pre = text_pre
		this.rules = rules
        this.text_post = text_post
        this.value_suffix = value_suffix
	}
}

Vue.component('threshold-alert', {
	props: {
		'threshold': {
            required: true,
            type: Object
        },
        'value': {
            required: true,
            type: Number
        }
	},
	methods: {
		getThresholdCategory: function (thresholds, eval_value) {
			let value_category = null ;
			const asc_thresholds = Object.keys(thresholds).sort(function (a, b) {  return parseFloat(a) - parseFloat(b)  })
			asc_thresholds.forEach(function(curr_threshold){
				if( parseFloat(curr_threshold) <= eval_value ){
					value_category = thresholds[curr_threshold]
				}
			});
			return( value_category )
		}
	},
	computed: {
		category_classes: function(){
			const category_classes = "alert alert-" + this.getThresholdCategory( this.threshold.rules, this.value )
			return( category_classes )
		}
	},
	template:
		`<div :class="category_classes" role="alert">
			{{threshold.text_pre}} <span class="right-badge">{{value}}{{threshold.value_suffix}}</span> {{threshold.text_post}}
        </div>`
})
