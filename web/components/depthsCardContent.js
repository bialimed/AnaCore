/*!
  * depthsCardContent component v1.1.0
  * Copyright 2018 IUCT-O
  * Author Frederic Escudie
  * Licensed under GNU General Public License
  */

Vue.component('depths-card-content', {
	props: {
        'shallowsAreas': {
            default: null,
            type: Object
        },
        'depthsMetrics': {
            default: null,
            type: Object
        }
    },
    data: function(){
        return {
            "threshold_shallows_prct": new ThresholdAlert(
                {1:"danger", 0.5:"warning", 0.001:"success", 0:"good"},
                "",
                "",
                "%"
            ),
            "threshold_shallows_with_variants": new ThresholdAlert(
                {1:"warning", 0:"good"},
                "",
                " shallows areas contains known variants potentially missed"
            )
        }
    },
    computed: {
        values_shallows_with_variants: function(){
            let nb_shallows_with_variants = null
            if(this.shallowsAreas !== null){
                nb_shallows_with_variants = 0
                this.shallowsAreas.results.forEach(function(curr_area){
                    if(curr_area.known_variants.length != 0){
                        nb_shallows_with_variants += 1
                    }
                })
            }
            return nb_shallows_with_variants
        },
        values_shallows_prct: function(){
            let prct_shallows_nt = null
            if(this.depthsMetrics !== null){
                this.threshold_shallows_prct.text_post = " of targeted nucleotids have a depth <= " +  this.depthsMetrics.under_threshold.threshold
                const spl_name = Object.keys(this.depthsMetrics.under_threshold.count_by_spl)[0]
                prct_shallows_nt = this.depthsMetrics.under_threshold.count_by_spl[spl_name].rate * 100
            }
            return prct_shallows_nt
        }
    },
	template:
        `<div>
            <div class="row">
                <div class="col-md-10">
                    <template v-if="depthsMetrics !== null">
                        <depths-chart
                            :analysis="depthsMetrics">
                        </depths-chart>
                    </template>
                </div>
                <div class="col-md-2">
                    <template v-if="values_shallows_prct !== null">
                        <threshold-alert
                            :threshold=threshold_shallows_prct
                            :value=values_shallows_prct
                        >
                        </threshold-alert>
                    </template>
                </div>
            </div>
            <div class="row">
                <div class="col-md-10">
                    <template v-if="shallowsAreas !== null">
                        <shallows-analysis-table
                            :analysis="shallowsAreas">
                        </shallows-analysis-table>
                    </template>
                </div>
                <div class="col-md-2">
                    <template v-if="values_shallows_with_variants !== null">
                        <threshold-alert
                            :threshold=threshold_shallows_with_variants
                            :value=values_shallows_with_variants
                        >
                        </threshold-alert>
                    </template>
                </div>
            </div>
        </div>`
})
