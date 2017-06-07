import sys
import itertools
from SVIO import *

def getOverlap( CNV_list ):
	support = list()
	contradict = list()
	prev_start = -1
	prev_end = 99999999999999999999999999999999999999999999999999
	for current_CNV in CNV_list:
		change_start = False
		change_end = False
		if current_CNV["start"] >= prev_start and current_CNV["start"] <= prev_end:
			change_start = True	
		if current_CNV["end"] <= prev_end and current_CNV["end"] >= prev_start:
			change_end = True
		if not change_start and not change_end:
			if current_CNV["start"] < prev_start and current_CNV["end"] > prev_end:
				support.append( current_CNV )
			else:
				contradict.append( current_CNV )
		elif change_start:
			prev_start = current_CNV["start"]
		elif change_end:
			prev_end = current_CNV["end"]
		if change_start or change_end:
			support.append( current_CNV )
	return({
		"overlap_start": prev_start,
		"overlap_end": prev_start,
		"support": support,
		"contradict": contradict
	})

for in_path in sys.argv[1:]:
	CNV = list()
	FH_in = SVIO(in_path, "r", ",")
	try:
		for record in FH_in:
			if record["CNV Id"].strip() != "":
				gene = "" if "Gene" not in record else record["Gene"]
				region, positions = record["CN Segment Posn."].split(":")
				start, end = positions.split("..")
				start = int(start)
				end = int(end)
				if start > end:
					old_start = start
					end = start
					start = old_start
				CNV.append({
					"region": region,
					"gene": gene,
					"start": start,
					"end": end
				})
	finally:
		FH_in.close()
    # Find best overlap
	best_overlap = {
		"overlap_start": None,
		"overlap_end": None,
		"support": list(),
		"contradict": list()
	}
	for current_combi in itertools.combinations( CNV, len(CNV) ):
		overlap_data = getOverlap( current_combi )
		if len(overlap_data["support"]) > len(best_overlap["support"]):
			best_overlap = overlap_data
		elif len(overlap_data["support"]) == len(best_overlap["support"]) and (overlap_data["overlap_end"] - overlap_data["overlap_start"]) > (best_overlap["overlap_end"] - best_overlap["overlap_start"]):
			best_overlap = overlap_data
	# Display
	print( in_path )
	print( "------------------------------------------------------------" )
	# Overlap
	support_ratio = str(len(best_overlap["support"])) + "/" + str(len(best_overlap["contradict"]) + len(best_overlap["support"]))
	print( "Common overlap (" + support_ratio + "):" )
	print( "", "Gene", "Chr", "Overlap_start", "Overlap_end", "Support_ratio", sep="\t" )
	print(
		"",
		best_overlap["support"][0]["gene"],
		best_overlap["support"][0]["region"],
		best_overlap["overlap_start"],
		best_overlap["overlap_end"],
		support_ratio,
		sep="\t"
	)
	# Excluded
	contradict_ratio = str(len(best_overlap["contradict"])) + "/" + str(len(best_overlap["contradict"]) + len(best_overlap["support"]))
	print( "Excluded (" + contradict_ratio + "):" )
	print( "", "start", "end", sep="\t" )
	for excluded in best_overlap["contradict"]:
		print( "", excluded["start"], excluded["end"], sep="\t" )
	print("")
