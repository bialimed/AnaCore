import sys
import json

depth_distrib_file = sys.argv[1]
strand = sys.argv[2]
breakpoint = sys.argv[3]

cov_by_tr = dict()
with open(depth_distrib_file) as FH_in:
    for line in FH_in:
        if not line.startswith("#"):
            region, start, end, name, sample, cov_min, cov_lower_quartile, cov_median, cov_upper_quartile, cov_max = [field.strip() for field in line.split("\t")]
            transcript = name.split(".")[0]
            if transcript not in cov_by_tr:
                cov_by_tr[transcript] = list()
            cov_by_tr[transcript].append({
                "start": start,
                "end": end,
                "depths": [float(cov_min), float(cov_lower_quartile), float(cov_median), float(cov_upper_quartile), float(cov_max)]
            })

depths_series = list()
for tr_id in sorted(cov_by_tr):
    tr_data = list()
    exons_pos = list()
    idx_bp_categories = None
    sorted_exons = sorted(cov_by_tr[tr_id], key=lambda x: (x["start"], x["end"]))
    if strand == "-":
        sorted_exons = sorted_exons[::-1]
    for idx_exon, exon in enumerate(sorted_exons):
        tr_data.append( exon["depths"] )
        if strand == "+":
            exons_pos.append( str(exon["start"]) + " - " + str(exon["end"]) )
        else:
            exons_pos.append( str(exon["end"]) + " - " + str(exon["start"]) )
        if breakpoint is not None:
            if int(breakpoint) == int(exon["start"]):
                idx_bp_categories = idx_exon - 0.25
                if strand == "-":
                    idx_bp_categories = idx_exon + 0.25
            elif int(breakpoint) == int(exon["end"]):
                idx_bp_categories = idx_exon + 0.25
                if strand == "-":
                    idx_bp_categories = idx_exon - 0.25
            elif int(breakpoint) > int(exon["start"]) and  int(breakpoint) < int(exon["end"]):
                idx_bp_categories = idx_exon
            elif int(breakpoint) < int(exon["start"]) and idx_bp_categories is None:
                idx_bp_categories = idx_exon - 0.5
    if breakpoint is not None and idx_bp_categories is None:
        idx_bp_categories = len(sorted_exons) - 0.5
    depths_series.append({
        "name": tr_id,
        "data": tr_data,
        "categories": exons_pos,
        "vertical_line": idx_bp_categories
    })

print( json.dumps(depths_series) )
