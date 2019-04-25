/*!
  * variantsTable component v2.1.0
  * Copyright 2018 IUCT-O
  * Author Frederic Escudie
  * Licensed under GNU General Public License
  */


Vue.filter('vtablePrct', function (value, precision=2) {
    return value.toFixed(precision)
})


Vue.filter('vtableFormatLocation', function(value){
    return value.split(":")[0] + ": " + value.split(":")[1]
})


Vue.filter('vtableFormatHGVS', function(value){
    let cleaned_value = ""
    if( value !== null && value != "" ){
        cleaned_value = value.split(":")[0] + ": " + value.split(":")[1]
    }
    return cleaned_value
})


function vtablesGetStrThresholdCategory(thresholds, eval_value) {
    let value_category = [] ;
    if(thresholds.hasOwnProperty(eval_value)){
        value_category = [
            "vtable-label",
            "vtable-label-" + thresholds[eval_value]
        ]
    } else {
        value_category = [
            "vtable-label",
            "vtable-label-na"
        ]
    }
    return value_category
}


function vtablesGetNumThresholdCategory(thresholds, eval_value) {
    let value_category = [] ;
    const asc_thresholds = Object.keys(thresholds).sort(function (a, b) {  return parseFloat(a) - parseFloat(b)  })
    asc_thresholds.forEach(function(curr_threshold){
        if( parseFloat(curr_threshold) <= parseFloat(eval_value) ){
            value_category = [
                "vtable-label",
                "vtable-label-" + thresholds[curr_threshold]
            ]
        }
    })
    return value_category
}


const vtable_pred_thresholds = {
    "CADD_phred": {"type": "num", "thresholds": {0:"danger", 15:"warning", 25:"good"}},
    "MetaLR": {"type": "num", "thresholds": {0:"danger", 0.80:"warning", 0.90:"good"}},
    "VEST3": {"type": "num", "thresholds": {0:"danger", 0.80:"warning", 0.90:"good"}},
    "ClinVar": {
        "type": "str", "thresholds": {
            "benign": "danger",
            "likely_benign": "danger",
            "pathogenic": "good",
            "likely_pathogenic": "good",
            "drug_response": "good",
            "association": "warning",
            "risk_factor": "warning",
            "protective": "warning",
            "uncertain_significance": "na",
            "not_provided": "na",
            "other": "na",
            "conflicting_data_from_submitters": "na",
            "affects": "na"
        }
    }
}


const vtable_header = [
    {'title': 'Location', '_val': function(variant, annot){return variant.coord.cpltPos}},
    {'title': 'Gene', '_val': function(variant, annot){return annot.subject.symbol}},
    {'title': 'HGVSc', '_val': function(variant, annot){return annot.changes.HGVSc}},
    {'title': 'HGVSp', '_val': function(variant, annot){return annot.changes.HGVSp}},
    {'title': 'Source', '_val': function(variant, annot){return variant.sources}},
    {
        'title': 'AF',
        'type_fct': parseFloat,
        '_val': function(variant, annot, default_source=null){
            support = null
            if(default_source !== null && variant.supports.hasOwnProperty(default_source)){
                support = variant.supports[default_source]
            } else {
                support = variant.supports[variant.sources[0]]
            }
            return support.AF * 100
        }
    },
    {
        'title': 'Depth',
        'type_fct': parseInt,
        '_val': function(variant, annot, default_source=null){
            support = null
            if(default_source !== null && variant.supports.hasOwnProperty(default_source)){
                support = variant.supports[default_source]
            } else {
                support = variant.supports[variant.sources[0]]
            }
            return support.depth
        }
    },
    {
        'title': 'Ratio lib',
        '_val': function(variant, annot, default_source=null){
            let value = []
            support = variant.supports[variant.sources[0]]
            support.libraries.forEach(function(curr_lib) {
                value.push(
                    curr_lib.alt_depth + "/" + curr_lib.depth
                )
            })
            return value.join(", ")
        }
    },
    {'title': 'Filters', '_val': function(variant, annot){return variant.filters}},
    {
        'title': 'Max frequency',
        'type_fct': parseFloat,
        '_val': function(variant, annot){
            let MAF = variant.MAF_pop
            if(MAF !== null){
                MAF = MAF.AF * 100
            }
            return MAF
        }
    },
    {'title': 'Conseq', '_val': function(variant, annot){return annot.conseq.split("&")}},
    {'title': 'Known', '_val': function(variant, annot){return variant.xref}},
    {
        'title': 'Transcript',
        '_val': function(variant, annot){
            if(annot.changes.HGVSc != null){
                return annot.changes.HGVSc.split(".")[0]
            }
            return ""
        }
    },
    {
        'title': 'Pathogenicity',
        '_val': function(variant, annot){
            return annot.pathogenicity
        }
    },
    {'title': 'Chromosome', '_val': function(variant, annot){return variant.coord.region}},
    {'title': 'RefStart', 'type_fct': parseFloat, '_val': function(variant, annot){return variant.coord.refStart()}},
    {'title': 'RefEnd', 'type_fct': parseFloat, '_val': function(variant, annot){return variant.coord.refEnd()}},
]


VariantsTable = Vue.component('variants-table', {
    extends: DynamicTable,
    props: {
        data: {
            default: null,
            type: Array, // List of Variants objects
        },
        filters: {
            default: null,
            type: Object, // An instance of FiltersCombiner
        },
        pred_thresholds: {
            type: Object,
            default: function(){ return vtable_pred_thresholds }
        },
        default_source: {
            type: String,
            default: "unknown"
        },
        header: { // Immutable property
            type: Array,
            default: function(){ return vtable_header }
        },
        menu: {
            default: function(){ return {"export": true, "filter": true, "pagination": true} },
            type: Object
        },
        scroll_x: {
            default: true,
            type: Boolean
        },
        title: {
            default: "Variants found",
            type: String
        },
        export_title: String
    },
    computed: {
        filteredData: function () {
            const self = this
            let filtered_data = []
            this.data.forEach(function( curr_variant ){
                //let unique_val = {}
                const variant_rows = self.variantToRows(curr_variant)
                variant_rows.forEach(function( curr_row ){
                      if( self.filters == null || self.filters.eval(curr_row) ){
                        const row_id = JSON.stringify(Object.values(curr_row))
                        //if(!unique_val.hasOwnProperty(row_id)){
                            filtered_data.push(curr_row)
                        //    unique_val[row_id] = 1
                        //}
                      }
                })
              })
            return filtered_data
        },
        datatableParams: function(){
            order_rule = [[ 0, "asc" ]]
            this.columns.forEach(function( curr_col, col_idx ){
                if( curr_col.hasOwnProperty("sort") ){
                    order_rule = [[ col_idx, curr_col["sort"] ]]
                }
            })
            let datatable_params = {
                dom: '',
                order: order_rule
            }
            if( this.scroll_x ){
                datatable_params["scrollX"] = true
                datatable_params["fixedColumns"] = {
                    "leftColumns": 2
                }
            }
            if( this.menu.pagination ){
                if( this.menu.filter ){
                    datatable_params["dom"] = 'f'
                }
                if(this.menu.export){
                    datatable_params["dom"] = "<'d-md-none'<'row'<'col-sm-12'B><'col-sm-12'l><'col-sm-12'" + datatable_params["dom"] + "r>>>" +
						"<'d-none d-md-block'<'row'<'col-md-6'<'row'<'dtable-head-btn'B>l>><'col-md-6'" + datatable_params["dom"] + "r>>>" +
                        "<'row'<'col-sm-12't>>" +
                        "<'row'<'col-sm-12 col-md-5'i><'col-sm-12 col-md-7'p>>"
                } else {
                    datatable_params["dom"] = "<'row'<'col-sm-12 col-md-6'l><'col-sm-12 col-md-6'" + datatable_params["dom"] + "r>>" +
                        "<'row'<'col-sm-12't>>" +
                        "<'row'<'col-sm-12 col-md-5'i><'col-sm-12 col-md-7'p>>"
                }
            }
            if( this.menu.export ){
                if(this.export_title === null){
                    datatable_params["buttons"] = ['csv', 'excel']
                } else {
                    datatable_params["buttons"] = [{
                        extend: 'csv',
                        title: this.export_title
                    },{
                        extend: 'excel',
                        title: this.export_title
                    }]
                }
            }
            return datatable_params
        },
    },
    methods: {
        variantToRows: function(variant){
            let entries = []
            if( variant.annot.length == 0 ){ // The variant has no annotation
                let curr_entry = {}
                const fake_annot = {
                    "subject": {"symbol": ""},
                    "changes": {"HGVSc": "", "HGVSp": ""},
                    "conseq": ""
                }
                this.columns.forEach(function( curr_col ){
                    curr_entry[curr_col.title.split(' ').join('_')] = curr_col._val(variant, fake_annot, this.default_source)
                })
                curr_entry._variant = variant
                entries = [curr_entry]
            } else { // The variant has annotation(s)
                const self = this
                variant.annot.forEach(function( curr_annot ){
                    let curr_entry = {}
                    self.columns.forEach(function( curr_col ){
                        curr_entry[curr_col.title.split(' ').join('_')] = curr_col._val(variant, curr_annot, self.default_source)
                    })
                    curr_entry._variant = variant
                    entries.push(curr_entry)
                })
            }
            return entries
        },
        predScoreClasses: function(predictor, score){
            let classes = []
            if(score !== null){
                classes = ["vtable-label", "vtable-label-na"]
                if(this.pred_thresholds.hasOwnProperty(predictor)){
                    const pred_param = this.pred_thresholds[predictor]
                    if(pred_param.type == "str"){
                        classes = vtablesGetStrThresholdCategory(pred_param.thresholds, score)
                    } else {
                        classes = vtablesGetNumThresholdCategory(pred_param.thresholds, score)
                    }
                }
            }
            return classes
        },
        split: function(val, separator="&"){
            let splitted = [null]
            if(val !== null){
                splitted = val.split(separator)
            }
            return splitted
        }
    },
    template:
        `<table class="table" :style="datatableStyle">
            <caption>
                {{title}}
            </caption>
            <thead>
                <tr>
                    <th>Location</th>
                    <th>Gene</th>
                    <th>HGVSc</th>
                    <th>HGVSp</th>
                    <th>Source</th>
                    <th>AF</th>
                    <th>Depth</th>
                    <th>Ratio in libraries</th>
                    <th>Filters</th>
                    <th>Max frequency (%)</th>
                    <th>Conseq</th>
                    <th>Pathogenicity</th>
                    <th>Known</th>
                </tr>
            </thead>
            <tbody>
                <tr v-for="entry in filteredData">
                    <td>{{entry.Location | vtableFormatLocation}}</td>
                    <td>{{entry.Gene}}</td>
                    <td>{{entry.HGVSc | vtableFormatHGVS}}</td>
                    <td>{{entry.HGVSp | vtableFormatHGVS}}</td>
                    <td>
                        <template v-for="curr_source in entry.Source">
                            {{curr_source}}
                        </template>
                    </td>
                    <td :class="{'vtable-noStdCaller': !entry._variant.supports.hasOwnProperty(default_source)}">{{entry.AF | vtablePrct}}</td>
                    <td :class="{'vtable-noStdCaller': !entry._variant.supports.hasOwnProperty(default_source)}">{{entry.Depth}}</td>
                    <td :class="{'vtable-noStdCaller': !entry._variant.supports.hasOwnProperty(default_source)}">{{entry.Ratio_lib}}</td>
                    <td>
                        <span v-for="curr_filter in entry.Filters" class="block-list">
                            <span v-if='curr_filter == "PASS"' class="vtable-label vtable-label-good">{{curr_filter}}</span>
                            <span v-else-if='curr_filter == "popAF"' class="vtable-label vtable-label-warning">{{curr_filter}}</span>
                            <span v-else class='vtable-label vtable-label-danger'>{{curr_filter}}</span>
                        </span>
                    </td>
                    <td v-if="entry.Max_frequency !== null">{{entry.Max_frequency | vtablePrct(3)}}</td><td v-else></td>
                    <td>{{entry.Conseq.join(" ")}}</td>
                    <td>
                        <span v-for="(score, pred) in entry.Pathogenicity" class="block-list">
                            <span>
                                {{pred}}:
                                <template v-for="curr_score in split(score)">
                                    <span :class="predScoreClasses(pred, curr_score)">{{curr_score}}</span>
                                </template>
                            </span><br />
                        </span>
                    </td>
                    <td>
                        <template v-for="(xref_ids, db) in entry.Known">
                            <template v-if="db == 'cosmic'">
                                <a v-for="id in xref_ids" target="_blank" :href="'https://cancer.sanger.ac.uk/cosmic/search?q=' + id">{{id}} </a>
                            </template>
                            <template v-else-if="db == 'dbSNP'">
                                <a v-for="id in xref_ids" target="_blank" :href="'https://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=' + id">{{id}} </a>
                            </template>
                            <template v-else>
                                <span v-for="id in xref_ids">{{id}} </span>
                            </template>
                        </template>
                    </td>
                </tr>
            </tbody>
        </table>`
})
