/*!
  * MSI v1.0.0
  * Copyright 2018 IUCT-O
  * Author Frederic Escudie
  * Licensed under GNU General Public License
  */
class MSIReport {
	constructor(samples=null) {
		this.samples = (samples === null ? [] : samples)
	}

	getMethods(){
		let methods = new Set()
		this.samples.forEach(function(curr_spl){
			curr_spl.getMethods().forEach(function(curr_method){
				methods.add(curr_method)
			})
		})
		methods = Array.from(methods).sort()
		return methods
	}

	getLocusNameByLocusID(){
		let name_by_ID = {}
		this.samples.forEach(function(curr_spl){
			for(let locus_id in curr_spl.loci){
				const curr_locus = curr_spl.loci[locus_id]
				name_by_ID[locus_id] = curr_locus.name
			}
		})
		return name_by_ID
	}

	getLoci(){
		let loci = new Set()
		this.samples.forEach(function(curr_spl){
			curr_spl.getLoci().forEach(function(curr_locus){
				loci.add(curr_locus)
			})
		})
		loci = Array.from(loci).sort()
		return loci
	}

	getResByLocus(){
		const self = this
		let res_by_locus = {}
		self.getLoci().forEach(function(locus_id) {
			let locus_res = []
			self.samples.forEach(function(curr_spl) {
				locus_res.push(curr_spl.loci[locus_id])
			})
			res_by_locus[locus_id] = locus_res
		})
		return res_by_locus
	}

	getMostRepresentedByLocus(method, status){
		let peaks_by_locus = {}
		const results_by_locus = this.getResByLocus()
		for(let locus_id in results_by_locus){
			const locus_samples = results_by_locus[locus_id]
			let peaks = []
			locus_samples.forEach(function (curr_spl) {
				if( curr_spl.results.hasOwnProperty(method) ){
					const curr_locus_res = curr_spl.results[method]
					if( status.indexOf(curr_locus_res.status) != -1 ){
						const nb_by_len = curr_locus_res.data["nb_by_length"]
						let most_represented_len = null
						let max_count = -1
						for(let curr_len in nb_by_len){
							if(nb_by_len[curr_len] >= max_count){
								max_count = nb_by_len[curr_len]
								most_represented_len = curr_len
							}
						}
						peaks.push(most_represented_len)
					}
				}
			})
			peaks_by_locus[locus_id] = peaks
		}
		return peaks_by_locus
	}

	static fromJSON(list) {
		let obj = new MSIReport()
		list.forEach(function (spl) {
			obj.samples.push(
				MSISample.fromJSON(spl)
			)
		});
		return obj
	}
}


class MSISample {
	constructor(name, loci=null, results=null) {
		this.name = name
		this.loci = (loci === null ? {} : loci)
		this.results = (results === null ? {} : results)
	}

	getMajorityStatus(){
		let nb_methods = 0
		let count_by_status = {"MSI": 0, "MSS": 0, "Undetermined": 0}
		// Get most-represented status
		for(let curr_method in this.results){
			nb_methods += 1
			const curr_status = this.results[curr_method].status
			count_by_status[curr_status] += 1
		}
		const desc_status = Object.keys(count_by_status).sort(function(a, b){
			return count_by_status[a] - count_by_status[b]
		}).reverse()
		const most_represented_status = desc_status[0]
		// Get consensus status
		let consensus_status = null
		if( count_by_status[most_represented_status] > 0.5*nb_methods ){
			consensus_status = most_represented_status
		}
		return consensus_status
	}

	getLoci(){
		let loci = new Set()
		for(let curr_locus in this.loci){
			loci.add(curr_locus)
		}
		loci = Array.from(loci).sort()
		return loci
	}

	getMethods(){
		let methods = new Set()
		for(let curr_method in this.results){
			methods.add(curr_method)
		}
		methods = Array.from(methods).sort()
		return methods
	}

	static fromJSON(hash) {
		let obj = new MSISample(null)
		Object.keys(hash).forEach(function (key) {
			obj[key] = hash[key]
		});
		return obj
	}
}
