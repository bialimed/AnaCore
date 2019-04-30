/*!
  * runDensityCardContent component v1.0.0
  * Copyright 2018 IUCT-O
  * Author Frederic Escudie
  * Licensed under GNU General Public License
  */


const rdensity_spec_by_instrument = {
    /*
     * sources:
     *   https:///www.illumina.com/content/dam/illumina-marketing/documents/products/other/miseq-overclustering-primer-770-2014-038.pdf
     *   https:///www.illumina.com/bulletins/2016/10/cluster-density-guidelines-for-illumina-sequencing-platforms-.html
     */
    "MiniSeq": {220:"error", 170:"good", 119:"warning", 0:"error"},
    "MiSeq_v2": {1200:"error", 1000:"good", 700:"warning", 0:"error"},
    "MiSeq_v3": {1400:"error", 1200:"good", 840:"warning", 0:"error"},
    "NextSeq": {220:"error", 170:"good", 119:"warning", 0:"error"},
    "HiSeq_v2": {1000:"error", 850:"good", 595:"warning", 0:"error"},
    "HiSeq_v3": {850:"error", 750:"good", 525:"warning", 0:"error"},
    "HiSeq_v4": {1050:"error", 950:"good", 665:"warning", 0:"error"}
}


Vue.component('run-density-card-content', {
    props: {
        'interop_summary': {
            default: null,
            type: Object
        },
        'run_info': {
            default: null,
            type: Object
        },
        'density_spec': {
            default: function(){ return rdensity_spec_by_instrument },
            type: Object
        }
    },
    data: function(){
        return {
            "threshold_ratio_PF": new ThresholdAlert(
                {20:"danger", 10:"warning", 5:"success", 0:"good"},
                "",
                "of clusters do not pass filters",
                "%"
            )
        }
    },
    computed: {
        threshold_density: function(){
            let threshold_density = null
            if( this.run_info !== null ){
                let spec_id = this.run_info.instrument.platform
                if( spec_id == "MiSeq" || spec_id == "HiSeq" ){
                    spec_id += "_" +  this.run_info.kit.version
                }
                if( this.density_spec.hasOwnProperty(spec_id) ){
                    threshold_density = new ThresholdAlert(
                        this.density_spec[spec_id],
                        "Run density",
                        "",
                        " KC/mm2"
                    )
                }
            }
            return threshold_density
        },
        value_density: function(){
            return (this.interop_summary.density["50_percentile"] / 1000).toFixed(0)
        },
        value_ratio_PF: function(){
            return (
                (1 - (this.interop_summary.density_PF["50_percentile"]/this.interop_summary.density["50_percentile"])) * 100
            ).toFixed(1)
        }
    },
    template:
        `<div>
            <div class="row" v-if="interop_summary !== null">
                <div class="col-md-2"></div>
                <div class="col-md-6">
                    <run-density-graph
                        :analysis=interop_summary>
                    </run-density-graph>
                </div>
                <div class="col-md-2">
                    <threshold-alert v-if="threshold_density !== null"
                        :threshold=threshold_density
                        :value=value_density>
                    </threshold-alert>
                    <threshold-alert
                        :threshold=threshold_ratio_PF
                        :value=value_ratio_PF>
                    </threshold-alert>
                </div>
                <div class="col-md-2"></div>
            </div>
        </div>`
})
