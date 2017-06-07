import os
import sys

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
LIB_DIR = os.path.abspath(os.path.join(os.path.dirname(CURRENT_DIR), "lib"))
sys.path.append(LIB_DIR)
if os.getenv('PYTHONPATH') is None: os.environ['PYTHONPATH'] = LIB_DIR
else: os.environ['PYTHONPATH'] = os.environ['PYTHONPATH'] + os.pathsep + LIB_DIR

from GTFI import *
from VEPvcf import VEPVCFIO, getAlleleRecord

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


min_threshold = 0.05

genes_map = getGenesByChrom( sys.argv[1] )

samples = list()
genes_with_variants = dict()
for idx_vcf, current_vcf in enumerate(sys.argv[2:]):
    with VEPVCFIO(current_vcf) as FH_vcf:
        for idx_spl, spl_name in enumerate(FH_vcf.samples):
            samples.append( spl_name )
            for record in FH_vcf:
                overlapping_genes = [gene["attr"]["gene_id"] for gene in getOverlappingGene(genes_map, record.chrom, record.pos)]
                for idx_alt, alt in enumerate(record.alt):
                    allele_record = getAlleleRecord( FH_vcf, record, idx_alt )
                    # Get sample depth
                    spl_depth = None
                    if len(record.samples) == 1 and "DP" in allele_record.info:
                        spl_depth = allele_record.info["DP"]
                    elif "DP" in allele_record.samples[spl_name]:
                        spl_depth = allele_record.samples[spl_name]["DP"]
                    else:
                        raise Exception( 'The sample depth cannot be retrieved for variant "' + record.chrom + ":" + str(record.pos) + '" in file "' + current_vcf + '".' )
                    # Get allele frequency
                    allele_freq = None
                    if "AF" in allele_record.samples[spl_name]:
                        frequencies = allele_record.samples[spl_name]["AF"]
                        print(frequencies)
                        if FH_vcf.info["AF"]["number_tag"] == "R": 
                            allele_freq = frequencies[1]
                        else:
                            allele_freq = frequencies[0]
                    elif "AD" in allele_record.samples[spl_name]:
                        frequencies = allele_record.samples[spl_name]["AD"]
                        if FH_vcf.format["AD"]["number_tag"] == "R":
                            allele_freq = frequencies[1]/float(spl_depth)
                        else:
                            allele_freq = frequencies[0]/float(spl_depth)
                    elif len(record.samples) <= 1 and "AF" in record.info:
                        allele_freq = record.info["AF"]
                        if isinstance(record.info["AF"], list) and len(record.info["AF"]) == 1:
                            allele_freq = record.info["AF"][0]
                        else:
                            raise Exception( 'The allele frequency cannot be retrieved for variant "' + record.chrom + ":" + str(record.pos) + '" in file "' + current_vcf + '".' )
                    else:
                        raise Exception( 'The allele frequency cannot be retrieved for variant "' + record.chrom + ":" + str(record.pos) + '" in file "' + current_vcf + '".' )
                    # Store
                    if allele_freq >= min_threshold:
                        for gene_id in overlapping_genes:
                            if gene_id not in genes_with_variants:
                                genes_with_variants[gene_id] = dict()
                            if spl_name not in genes_with_variants[gene_id]:
                                genes_with_variants[gene_id][spl_name] = 0
                            genes_with_variants[gene_id][spl_name] += 1

print( "Gene_ID", "\t".join(samples), sep="\t" )
for gene_id in genes_with_variants:
    counts = list()
    for spl_name in samples:
        if spl_name in genes_with_variants[gene_id]:
            counts.append( genes_with_variants[gene_id][spl_name] )
        else:
            counts.append( 0 )
    print(
        gene_id,
        "\t".join(map(str, counts)),
        sep="\t"
    )
