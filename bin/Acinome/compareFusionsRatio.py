import os
import sys
from SVIO import *

spl_TT_path = sys.argv[2]
spl_TN_path = sys.argv[1]

TN_genes = dict()
FH_TN = SVIO( sys.argv[1] )
try:
    for record in FH_TN:
        genes_A = record["Genes_A_overlapping_break"] ########################### problem same region but different gene (downstream one gene is not overlapping)
        genes_B = record["Genes_B_overlapping_break"]
        key = None
        if genes_A <= genes_B:
            key = genes_A + " @@-@-@@ " + genes_B
        else:
            key = genes_B + " @@-@-@@ " + genes_A
        if key not in TN_genes:
            TN_genes[key] = 0
        TN_genes[key] += int(record["Nb_spanning_reads"]) + int(record["Nb_spanning_pairs"]) + int(record["Nb_spanning_reads_and_pairs"])
finally:
    FH_TN.close()

TT_genes = dict()
FH_TT = SVIO( sys.argv[2] )
try:
    for record in FH_TT:
        genes_A = record["Genes_A_overlapping_break"] ########################### problem same region but different gene (downstream one gene is not overlapping)
        genes_B = record["Genes_B_overlapping_break"]
        key = None
        if genes_A <= genes_B:
            key = genes_A + " @@-@-@@ " + genes_B
        else:
            key = genes_B + " @@-@-@@ " + genes_A
        if key not in TT_genes:
            TT_genes[key] = 0
        TT_genes[key] += int(record["Nb_spanning_reads"]) + int(record["Nb_spanning_pairs"]) + int(record["Nb_spanning_reads_and_pairs"])
finally:
    FH_TT.close()

print("Gene_A\tGene_B\t" + os.path.basename(spl_TT_path) + "\t" + os.path.basename(spl_TN_path) )
for genes in sorted(TT_genes):
    nb_TT = TT_genes[genes]
    nb_TN = 1
    real_nb_TN = 0
    if genes in TN_genes:
        nb_TN = TT_genes[genes]
        real_nb_TN = nb_TN
    if (nb_TT > nb_TN) and (nb_TT/nb_TN > 3) and nb_TT > 200:
        print( "\t".join(genes.split(' @@-@-@@ ')) + "\t" + str(nb_TT) + "\t" + str(real_nb_TN) )
