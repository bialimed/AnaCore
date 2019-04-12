/*!
  * depthsCardContent component v1.2.0
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
            "threshold_missing_prct": new ThresholdAlert(
                {0.5:"danger", 0.001:"success", 0:"good"},
                "",
                "",
                "%"
            ),
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
                const spl_name = Object.keys(this.depthsMetrics.under_threshold.count_by_spl)[0]
                prct_shallows_nt = this.depthsMetrics.under_threshold.count_by_spl[spl_name].rate * 100
                this.threshold_shallows_prct.text_post = " (" + this.depthsMetrics.under_threshold.count_by_spl[spl_name].count + " nt) of targeted nucleotids have a depth <= " +  this.depthsMetrics.under_threshold.threshold
            }
            return prct_shallows_nt
        },
        values_missing_prct: function(){
            let prct_missing_nt = null
            if(this.depthsMetrics !== null){
                const spl_name = Object.keys(this.depthsMetrics.under_threshold.count_by_spl)[0]
                let nb_missing_nt = 0
                const idx_missing = this.depthsMetrics.depths_classes.depths_list.indexOf(0)
                if(idx_missing != -1){
                    nb_missing_nt = this.depthsMetrics.depths_classes.count_by_spl[spl_name][idx_missing]
                }
                const nb_targeted_nt = this.depthsMetrics.sequencing_by_spl[spl_name].nt_targeted
                prct_missing_nt = (nb_missing_nt * 100 / nb_targeted_nt).toFixed(2)
                this.threshold_missing_prct.text_post = " (" + nb_missing_nt + " nt) of targeted nucleotids are not covered"
            }
            return prct_missing_nt
        }
    },
    template:
        `<div>
            <div class="row">
                <div class="col-md-10">
                    <depths-chart v-if="depthsMetrics !== null"
                        :analysis="depthsMetrics">
                    </depths-chart>
                </div>
                <div class="col-md-2">
                    <threshold-alert v-if="values_missing_prct !== null"
                        :threshold=threshold_missing_prct
                        :value=values_missing_prct>
                    </threshold-alert>
                    <threshold-alert v-if="values_shallows_prct !== null"
                        :threshold=threshold_shallows_prct
                        :value=values_shallows_prct>
                    </threshold-alert>
                </div>
            </div>
            <div class="row">
                <div class="col-md-10">
                    <shallows-analysis-table v-if="shallowsAreas !== null"
                        :analysis="shallowsAreas">
                    </shallows-analysis-table>
                </div>
                <div class="col-md-2">
                    <threshold-alert v-if="values_shallows_with_variants !== null"
                        :threshold=threshold_shallows_with_variants
                        :value=values_shallows_with_variants>
                    </threshold-alert>
                </div>
            </div>
        </div>`
})
