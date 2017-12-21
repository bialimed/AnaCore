import os
import sys
import uuid
import argparse
import subprocess

CURRENT_DIR = os.path.dirname(__file__)
LIB_DIR = os.path.abspath(os.path.join(os.path.dirname(CURRENT_DIR), "lib"))
sys.path.append(LIB_DIR)
if os.getenv('PYTHONPATH') is None: os.environ['PYTHONPATH'] = LIB_DIR
else: os.environ['PYTHONPATH'] = os.environ['PYTHONPATH'] + os.pathsep + LIB_DIR

from VEPvcf import VEPVCFIO, getAlleleRecord


if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser( description='**********************************.' )
    group_input = parser.add_argument_group( 'Inputs' ) # Inputs
    group_input.add_argument( '--aln-files', nargs='+', required=True, help='**************************** (format: BAM).' )
    group_input.add_argument( '--ids-files', nargs='+', required=True, help='****************************.' )
    group_output = parser.add_argument_group( 'Outputs' ) # Outputs
    group_output.add_argument( '--out-dir', default=os.getcwd(), help='****************************.' )
    args = parser.parse_args()

    # Get variants
    variants = dict()
    for current_ids in args.ids_files:
        with open(current_ids) as FH_ids:
            for line in FH_ids:
                pos_id, allele_id = line.strip().split("=")
                chrom, pos = pos_id.split(":")
                if pos_id not in variants:
                    variants[pos_id] = {
                        "chrom": chrom,
                        "pos": pos,
                        "alleles": [],
                        "coverage": dict()                   
                    }
                variants[pos_id]["alleles"].append( allele_id )

    # Variants to BED
    tmp_bed = os.path.join( args.out_dir, "tmp_" + str(uuid.uuid1()) + ".bed" )
    with open(tmp_bed, "w") as FH_bed:
        for pos_id in variants:
            current_var = variants[pos_id]
            FH_bed.write(
                "\t".join([
                    current_var["chrom"],
                    str(int(current_var["pos"]) - 1),
                    current_var["pos"]
                ]) + "\n"
            )

    # Positions coverage
    samples = list()
    for spl_aln in args.aln_files:
        spl = os.path.basename(spl_aln).split("_")[0]
        samples.append(spl)
        tmp_coverage = os.path.join( args.out_dir, "tmp_" + str(uuid.uuid1()) + ".tsv" )
        # Process depth
        cmd = "/save/fescudie/softwares/samtools-1.3.1/samtools depth -aa -b " + tmp_bed + " -m 1000000 " + spl_aln + " > " + tmp_coverage
        subprocess.check_call(cmd, shell=True)
        # Store depth
        with open(tmp_coverage) as FH_coverage:
            for line in FH_coverage:
                chrom, pos, depth = line.strip().split("\t")
                pos_id = chrom + ":" + pos
                if spl in variants[pos_id]["coverage"]:
                    raise Exception("The depth for variant '" + pos_id + "' already exists.")
                variants[pos_id]["coverage"][spl] = depth
        os.remove( tmp_coverage )
    os.remove( tmp_bed )

    # Display results
    print( "#Chromosome", "Position", "ID", "\t".join(samples), sep="\t" )
    for pos_id in variants:
        current_var = variants[pos_id]
        print(
            current_var["chrom"],
            current_var["pos"],
            ",".join(current_var["alleles"]),
            "\t".join([current_var["coverage"][spl] for spl in samples]),
            sep="\t"
        )
