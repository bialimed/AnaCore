/*!
  * dynamicTable component v1.1.0
  * Copyright 2018 IUCT-O
  * Author Frederic Escudie
  * Licensed under GNU General Public License
  */


DynamicTable = Vue.component('dynamic-table', {
	props: {
		data: Array,
		filters: Object,  // An instance of FiltersCombiner
		header: {
			required: true,
			type: Array
		},
		// [  // Header row
			// {  // First column
			// 	'title': "str",
			// 	'is_html': "bool",
			// 	'class': "str",
			// 	'entryClass': fct(entry, col),
			// 	'sort': ["asc", "desc"],  // only on columns without colpsan above
			// 	'href': fct(entry, col),
			// 	'value': fct(entry, col),
			// 	'colspan': "int",  // only in multiple row header and not in last row
			// 	'rospan': "int",  // only in multiple row header and not in last row
			// },
		// 	{  // Second column
		// 		...
		// 	}
		// ]
		menu: {
			default: function(){ return {"export": true, "filter": true, "pagination": true} },
			type: Object
		},
		export_title:  {
			default: null,
			type: String
		},
		scroll_x: {
			default: false,  // This option is incompatible with header containing rowspan or colspan
			type: Boolean
		},
		title: String,
	},
	mounted: function() {
		$(this.$el).DataTable(this.datatableParams)
	},
	computed: {
		is_multi_rows_header: function(){
			return this.header.length !== 0 && Array.isArray(this.header[0])
		},
		nb_columns: function(){
			let nb_col = 0
			let first_row = this.header
			if( this.is_multi_rows_header ){
				first_row = this.header[0]
			}
			first_row.forEach(function(curr_col){
				if( !curr_col.hasOwnProperty("colspan") ){
					nb_col += 1
				} else {
					nb_col += curr_col.colspan
				}
			})
			return nb_col
		},
		header_mask: function(){
			let mask = this.header
			if( this.is_multi_rows_header ){
				const nb_row = this.header.length
				const nb_col = this.nb_columns
				mask = new Array(nb_row).fill([])
				mask.forEach(function(curr_row, row_idx){
					mask[row_idx] = new Array(nb_col).fill(null)
				})
				this.header.forEach(function(curr_row, row_idx){
					col_idx = 0
					curr_row.forEach(function(curr_title){
						while( mask[row_idx][col_idx] !== null ){
							col_idx += 1
						}
						colspan = (curr_title.hasOwnProperty("colspan") ? curr_title.colspan : 1)
						rowspan = (curr_title.hasOwnProperty("rowspan") ? curr_title.rowspan : 1)
						for(let sub_row_idx = row_idx ; sub_row_idx < row_idx + rowspan ; sub_row_idx++){
							for(let sub_col_idx = col_idx ; sub_col_idx < col_idx + colspan ; sub_col_idx++){
								mask[sub_row_idx][sub_col_idx] = curr_title
							}
						}
						col_idx += 1
					})
				})
			}
			return mask
		},
		columns: function(){
			let returned_cols = this.header_mask
			if( this.is_multi_rows_header ){
				returned_cols = returned_cols[returned_cols.length - 1]
			}
			return returned_cols
		},
		filteredData: function () {
			let filtered_data = this.data
			if( this.filters != null ){
				const self = this
				filtered_data = filtered_data.filter(function (row) {
					return self.filters.eval( row )
				})
			}
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
			}
			if( this.menu.filter ){
				datatable_params["dom"] = 'f'
			}
			if( this.menu.pagination ){
				datatable_params["dom"] = 'l' + datatable_params["dom"] + 'rtip'
			}
			if( this.menu.export ){
				if(this.export_title === null){
					datatable_params["buttons"] = ['csv', 'excel']
				} else {
					datatable_params["buttons"] = [
						{extend: 'csv', title: this.export_title},
						{extend: 'excel', title: this.export_title}
					]
				}
				datatable_params["dom"] = 'B' + datatable_params["dom"]
			}
			return datatable_params
		},
		_datatableStyle: function(){
			let style = {}
			if( this.scroll_x ){
				style = {
					"width": "100%",
					"border-bottom": "none"
				}
			}
			return style
		}
	},
	beforeUpdate: function() {
		$(this.$el).DataTable().destroy()
	},
	updated: function() {
		$(this.$el).DataTable(this.datatableParams)
	},
	methods: {
		getValue: function( entry, col ){
			return this.$options.filters.getValue( entry, col )
		},
	},
	filters: {
		getValue: function( entry, col ){
			let formatted_entry = null
			if( col.hasOwnProperty("value") ){
				formatted_entry = col["value"]( entry )
			} else if( entry.hasOwnProperty(col['title']) ){
				formatted_entry = entry[col['title']]
			}
			return formatted_entry
		},
		getHref: function( entry, col ){
			return col["href"](entry)
		},
		getClass: function( entry, col ){
			let classes = null
			if( col.hasOwnProperty("class") ){
				classes = [col["class"]]
			}
			if( col.hasOwnProperty("entryClass") ){
				added_class = col["entryClass"](entry)
				if( classes === null ){
					classes = [added_class]
				} else {
					classes.push(added_class)
				}
			}
			return classes
		},
		getRowSpan: function(col){
			return col.hasOwnProperty("rowspan") ? col.rowspan : 1
		},
		getColSpan: function(col){
			return col.hasOwnProperty("colspan") ? col.colspan : 1
		}
	},
	template:
		`<table class="table table-striped" :style="_datatableStyle">
			<caption>
				{{title}}
			</caption>
			<thead>
				<template v-if="!is_multi_rows_header">
					<tr>
						<th v-for="col in columns">{{col.title}}</th>
					</tr>
				</template>
				<template v-else>
					<tr v-for="row in header">
						<th v-for="col in row"
							:rowspan="col | getRowSpan"
							:colspan="col | getColSpan">
							{{col.title}}
						</th>
					</tr>
				</template>
			</thead>
			<tbody>
				<tr v-for="entry in filteredData">
					<td v-for="(col, col_idx) in columns" :class="entry | getClass(col)">
						<template v-if="col.is_html">
							<span v-html="getValue(entry, col)"></span>
						</template>
						<template v-else>
							<template v-if="col.href">
								<template v-if="col.href_target">
									<a :target="col.href_target" :href="entry | getHref(col)">{{entry | getValue(col)}}</a>
								</template>
								<template v-else>
									<a :href="entry | getHref(col)">{{entry | getValue(col)}}</a>
								</template>
							</template>
							<template v-else>
								{{entry | getValue(col)}}
							</template>
						</template>
					</td>
				</tr>
			</tbody>
		</table>`
})
