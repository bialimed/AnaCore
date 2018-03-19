import os
import sys

CURRENT_DIR = os.path.dirname(__file__)
LIB_DIR = os.path.abspath(os.path.join(os.path.dirname(CURRENT_DIR), "lib"))
sys.path.append(LIB_DIR)
if os.getenv('PYTHONPATH') is None: os.environ['PYTHONPATH'] = LIB_DIR
else: os.environ['PYTHONPATH'] = os.environ['PYTHONPATH'] + os.pathsep + LIB_DIR

from VEPvcf import VEPVCFIO, getAlleleRecord

samples = list()
variants = dict()
for idx_vcf, current_vcf in enumerate(sys.argv[1:]):
    with VEPVCFIO(current_vcf) as FH_vcf:
        for idx_spl, spl_name in enumerate(FH_vcf.samples):
            samples.append( spl_name )
            for record in FH_vcf:
                for idx_alt, alt in enumerate(record.alt):
                    allele_record = getAlleleRecord( FH_vcf, record, idx_alt )
                    variant_id = record.chrom + ":" + str(record.pos) + "=" + record.ref + "/" + alt
                    if variant_id not in variants:
                        variants[variant_id] = {
                            "chrom": record.chrom,
                            "pos": record.pos,
                            "ref": record.ref,
                            "alt": alt,
                            "freq": dict(),
                            "AD": dict()
                        }
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
                    variants[variant_id]["freq"][spl_name] = allele_freq
                    variants[variant_id]["AD"][spl_name] = spl_depth * allele_freq

print( "Variant_ID", "\t".join(samples), sep="\t" )
for variant_id in variants:
    counts = list()
    for spl_name in samples:
        if spl_name in variants[variant_id]["freq"]:
            counts.append( variants[variant_id]["freq"][spl_name] )
        else:
            counts.append( 0 )
    print(
        variant_id,
        "\t".join(map(str, counts)),
        sep="\t"
    )
