import os
import sys
from SVIO import *

spl_B_path = sys.argv[2]
spl_A_path = sys.argv[1]

TN_genes = set()
FH_TN = SVIO( sys.argv[1] )
try:
    for record in FH_TN:
        genes_A = record["Genes_A_overlapping_break"] ########################### problem same region but different gene (downstream one gene is not overlapping)
        genes_B = record["Genes_B_overlapping_break"]
        if genes_A <= genes_B:
            TN_genes.add( genes_A + " @@-@-@@ " + genes_B )
        else:
            TN_genes.add( genes_B + " @@-@-@@ " + genes_A )
finally:
    FH_TN.close()

TT_genes = set()
FH_TT = SVIO( sys.argv[2] )
try:
    for record in FH_TT:
        genes_A = record["Genes_A_overlapping_break"]
        genes_B = record["Genes_B_overlapping_break"]
        if genes_A <= genes_B:
            TT_genes.add( genes_A + " @@-@-@@ " + genes_B )
        else:
            TT_genes.add( genes_B + " @@-@-@@ " + genes_A )
finally:
    FH_TT.close()


TT_TN_common = TT_genes.intersection(TN_genes)
TT_specific = TT_genes.difference(TN_genes)
TN_specific = TN_genes.difference(TT_genes)
print("#Sample\tGene_A\tGene_B")
for genes in sorted(TT_specific):
	print( os.path.basename( spl_B_path ) + "\t" + "\t".join(genes.split(' @@-@-@@ ')) )
for genes in sorted(TN_specific):
	print( os.path.basename( spl_A_path ) + "\t" + "\t".join(genes.split(' @@-@-@@ ')) )
for genes in sorted(TT_TN_common):
	print( "All\t" +  "\t".join(genes.split(' @@-@-@@ ')) )
