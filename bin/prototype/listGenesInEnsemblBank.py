#!/usr/bin/env python3
#
# Copyright (C) 2017 IUCT-O
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2017 IUCT-O'
__license__ = 'GNU General Public License'
__version__ = '1.0.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import re
import sys

in_ensembl_fasta = sys.argv[1]

# Get nb transcripts by gene
count_by_gene_id = dict()
with open(in_ensembl_fasta) as FH:
    for line in FH:
        if line.startswith(">"):
            line = line.strip()
            match = re.search("gene:(ENSG[^\s]+)", line)
            gene_id = match.group(1).split(".")[0]
            if gene_id not in count_by_gene_id:
                count_by_gene_id[gene_id] = 0
            count_by_gene_id[gene_id] += 1

# Display genes
print("#Gene_stable_id", "Nb_transcripts", sep="\t")
for gene_id in sorted(count_by_gene_id):
    print(gene_id, count_by_gene_id[gene_id], sep="\t")
