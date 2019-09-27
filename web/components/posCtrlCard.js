/*!
  * posCtrlCard component v1.0.0
  * Copyright 2019 IUCT-O
  * Author Frederic Escudie
  * Licensed under GNU General Public License
  */
Vue.component('pos-ctrl-card-content', {
    props: {
        'pos_ctrl': {
            default: null,
            type: Array
        },
        'max_error_rate': {
            default: 0.33,
            type: Number
        },
        'max_sensitivity': {
            default: 0.05,
            type: Number
        },
        'max_error_value': {
            default: 0.03,
            type: Number
        },
    },
    computed: {
        threshold_lost_variants: function(){
            return new ThresholdAlert(
                {1:"warning", 0:"good"},
                "",
                "lost variants",
                "/" + this.nb_expected_variants
            )
        },
        threshold_lost_strong_variants: function(){
            return new ThresholdAlert(
                {1:"danger", 0:"good"},
                "",
                "lost variants with an expected AF >= " + this.max_sensitivity * 100 + "%",
                "/" + this.nb_expected_strong_variants
            )
        },
        threshold_error_strong_value: function(){
            return new ThresholdAlert(
                {0:"good", 1:"warning", 3:"danger"},
                "",
                "variants have an allele frequency error > " + this.max_error_value * 100 + "%"
            )
        },
        nb_expected_variants: function(){
            let nb_expected_vartiants = null
            if(this.pos_ctrl !== null){
                nb_lost_variants = 0
                this.pos_ctrl.forEach(function(curr_var){
                    if(curr_var.expected != 0){
                        nb_expected_vartiants += 1
                    }
                })
            }
            return nb_expected_vartiants
        },
        nb_lost_variants: function(){
            let nb_lost_variants = null
            if(this.pos_ctrl !== null){
                nb_lost_variants = 0
                this.pos_ctrl.forEach(function(curr_var){
                    if(curr_var.expected != 0 && curr_var.detected == 0){
                        nb_lost_variants += 1
                    }
                })
            }
            return nb_lost_variants
        },
        nb_expected_strong_variants: function(){
            let nb_expected_strong_variants = null
            if(this.pos_ctrl !== null){
                const self = this
                nb_expected_strong_variants = 0
                this.pos_ctrl.forEach(function(curr_var){
                    if(curr_var.expected >= self.max_sensitivity){
                        nb_expected_strong_variants += 1
                    }
                })
            }
            return nb_expected_strong_variants
        },
        nb_lost_strong_variants: function(){
            let nb_lost_strong_variants = null
            if(this.pos_ctrl !== null){
                const self = this
                nb_lost_strong_variants = 0
                this.pos_ctrl.forEach(function(curr_var){
                    if(curr_var.expected >= self.max_sensitivity && curr_var.detected == 0){
                        nb_lost_strong_variants += 1
                    }
                })
            }
            return nb_lost_strong_variants
        },
        nb_error_strong_value: function(){
            let nb_error_strong_value = null
            if(this.pos_ctrl !== null){
                const self = this
                nb_error_strong_value = 0
                this.pos_ctrl.forEach(function(curr_var){
                    if(Math.abs(curr_var.expected - curr_var.detected) > self.max_error_value){
                        nb_error_strong_value += 1
                    }
                })
            }
            return nb_error_strong_value
        }
    },
    template:
        `<div>
            <div class="row">
                <div class="col-md-10">
                    <pos-ctrl-table v-if="pos_ctrl !== null"
                        :data=pos_ctrl
                        :max_error_rate=max_error_rate
                        :scroll_x=true>
                    </pos-ctrl-table>
                </div>
                <div class="col-md-2" v-if="pos_ctrl !== null">
                    <threshold-alert
                        :threshold=threshold_lost_variants
                        :value=nb_lost_variants>
                    </threshold-alert>
                    <threshold-alert
                        :threshold=threshold_lost_strong_variants
                        :value=nb_lost_strong_variants>
                    </threshold-alert>
                    <threshold-alert
                        :threshold=threshold_error_strong_value
                        :value=nb_error_strong_value>
                    </threshold-alert>
                </div>
            </div>
        </div>`
})
