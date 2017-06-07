import sys
from GTFI import *
min_support = 200

def getGenesByChrom( gtf_file ):
    genes_by_chr = dict()
    FH_gtf = GTFI( gtf_file )
    for record in FH_gtf:
        if record["feature"] == "gene":
            if record["region"] not in genes_by_chr:
                genes_by_chr[record["region"]] = list()
            genes_by_chr[record["region"]].append( record )
    FH_gtf.close()
    for chrom in genes_by_chr:
        genes_by_chr[chrom] = sorted(genes_by_chr[chrom], key=lambda x: (x["start"], x["end"]))
    return( genes_by_chr )

def getOverlappingGene( genes_by_chr, chrom, pos ):
    overlapping_genes = list()
    if chrom in genes_by_chr:
        for gene in genes_by_chr[chrom]:
            if gene["start"] <= pos and gene["end"] >= pos:
                overlapping_genes.append( gene )
    return( overlapping_genes )


genes_by_chr = getGenesByChrom( sys.argv[1] )
print(
    "#Region_A_name",
    "Break_A_position",
    "Genes_A_overlapping_break",
    "Region_B_name",
    "Break_B_position",
    "Genes_B_overlapping_break",
    "Fusion_orientation",
    "Nb_spanning_reads",
    "Nb_spanning_pairs",
    "Nb_spanning_reads_and_pairs",
    "Nb_contradict_reads",
    sep="\t"
)
nb_fusions = 0
nb_filter_support = 0
nb_filter_chrom = 0
nb_filter_gene = 0
with open(sys.argv[2]) as FH_fusion_out:
    for line in FH_fusion_out:
        nb_fusions += 1
        # Parse
        fusion, trash_1, contig_a, contig_b, depth_a, depth_b, mate_distances = [field.strip() for field in line.strip().split("@")]
        chrom, break_a, break_b, orientation, nb_splitted_reads, nb_splitted_pairs, nb_pairs_splitted_reads, nb_contradict, base_cover_left, base_cover_right, trash_1 = [field.strip() for field in fusion.split("\t")]
        chrom_a, chrom_b = chrom.split("-")
        break_a = int(break_a)
        break_b = int(break_b)
        nb_splitted_reads = int(nb_splitted_reads)
        nb_splitted_pairs = int(nb_splitted_pairs)
        nb_pairs_splitted_reads = int(nb_pairs_splitted_reads)
        nb_support = nb_splitted_reads + nb_splitted_pairs + nb_pairs_splitted_reads
        nb_contradict = int(nb_contradict)
        base_cover_left = int(base_cover_left)
        base_cover_right = int(base_cover_right)
        # Filters
        kept = True
        #~ if nb_support < min_support:
            #~ nb_filter_support += 1
            #~ kept = False
        if chrom_a == chrom_b:
            nb_filter_chrom += 1
            kept = False
        #~ elif nb_splitted_reads < 10 or nb_splitted_pairs < 10 or nb_pairs_splitted_reads < 10:
            #~ kept = False
        genes_a = getOverlappingGene( genes_by_chr, chrom_a, break_a )
        genes_b = getOverlappingGene( genes_by_chr, chrom_b, break_b )
        if len(genes_a) < 1 and len(genes_b) < 1:
            nb_filter_gene += 1
            kept = False
        if kept:
            print( 
                chrom_a,
                break_a,
                ",".join([gene["attr"]["gene_name"] for gene in genes_a]),
                chrom_b,
                break_b,
                ",".join([gene["attr"]["gene_name"] for gene in genes_b]),
                orientation,
                nb_splitted_reads,
                nb_splitted_pairs,
                nb_pairs_splitted_reads,
                nb_contradict,
                sep="\t"
            )
#~ print(
    #~ nb_fusions,
    #~ nb_filter_support,
    #~ nb_filter_chrom,
    #~ nb_filter_gene,
    #~ sep="\n"
#~ )
