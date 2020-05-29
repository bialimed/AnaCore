#!/usr/bin/env python3

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2018 IUCT-O'
__license__ = 'GNU General Public License'
__version__ = '1.2.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import os
import sys
import uuid
import tempfile
import unittest

TEST_DIR = os.path.dirname(os.path.abspath(__file__))
PACKAGE_DIR = os.path.dirname(TEST_DIR)
sys.path.append(PACKAGE_DIR)

from anacore.gtf import GTFIO, loadModel
from anacore.region import RegionTree
from anacore.genomicRegion import Gene, Transcript, Protein, Exon, CDS


########################################################################
#
# FUNCTIONS
#
########################################################################
class TestGTFIO(unittest.TestCase):
    def setUp(self):
        tmp_folder = tempfile.gettempdir()
        unique_id = str(uuid.uuid1())

        # Temporary files
        self.tmp_ensembl_in_gtf = os.path.join(tmp_folder, unique_id + "_ensembl_in.gtf")
        self.tmp_ensembl_out_gtf = os.path.join(tmp_folder, unique_id + "_ensembl_out.gtf")
        self.tmp_ncbi_in_gtf = os.path.join(tmp_folder, unique_id + "_ncbi_in.gtf")
        self.tmp_ncbi_out_gtf = os.path.join(tmp_folder, unique_id + "_ncbi_out.gtf")

        # Ensembl GTF
        self.ensembl_expected = [
            Gene(54770583, 54771134, "+", "6", "KRASP1", {"id": "ENSG00000220635"}, None, [
                Transcript(54770583, 54771134, "+", "6", "KRASP1-201", {"id": "ENST00000407852"}, None, [
                    Exon(54770583, 54771134, "+", "6", "ENST00000407852_e1")
                ])
            ]),
            Gene(25204789, 25250936, "-", "12", "KRAS", {"id": "ENSG00000133703"}, None, [
                Transcript(25204789, 25250931, "-", "12", "KRAS-202", {"id": "ENST00000311936"}, None,
                    [
                        Exon(25250751, 25250931, "-", "12", "ENST00000311936_e1"),
                        Exon(25245274, 25245395, "-", "12", "ENST00000311936_e2"),
                        Exon(25227234, 25227412, "-", "12", "ENST00000311936_e3"),
                        Exon(25225614, 25225773, "-", "12", "ENST00000311936_e5"),
                        Exon(25204789, 25209911, "-", "12", "ENST00000311936_e6")
                    ], [
                        Protein(25209798, 25245384, "-", "12", None, None, None,
                            [
                                CDS(25245274, 25245384, "-", "12", ""),
                                CDS(25227234, 25227412, "-", "12", ""),
                                CDS(25225614, 25225773, "-", "12", ""),
                                CDS(25209798, 25209911, "-", "12", "")
                            ]
                        )
                    ]
                ),
                Transcript(25209168, 25250936, "-", "12", "KRAS-204", {"id": "ENST00000557334"}, None,
                    [
                        Exon(25250751, 25250936, "-", "12", "ENST00000557334_e1"),
                        Exon(25245274, 25245395, "-", "12", "ENST00000557334_e2"),
                        Exon(25209168, 25209911, "-", "12", "ENST00000557334_e3")
                    ], [
                        Protein(25209798, 25245384, "-", "12", None, None, None,
                            [
                                CDS(25245274, 25245384, "-", "12", ""),
                                CDS(25209798, 25209911, "-", "12", "")
                            ]
                        )
                    ]
                ),
                Transcript(25209431, 25250803, "-", "12", "KRAS-201", {"id": "ENST00000256078"}, None,
                    [
                        Exon(25250751, 25250803, "-", "12", "ENST00000256078_e1"),
                        Exon(25245274, 25245395, "-", "12", "ENST00000256078_e2"),
                        Exon(25227234, 25227412, "-", "12", "ENST00000256078_e3"),
                        Exon(25225614, 25225773, "-", "12", "ENST00000256078_e4"),
                        Exon(25215437, 25215560, "-", "12", "ENST00000256078_e5"),
                        Exon(25209431, 25209911, "-", "12", "ENST00000256078_e6")
                    ], [
                        Protein(25215444, 25245384, "-", "12", None, None, None,
                            [
                                CDS(25245274, 25245384, "-", "12", ""),
                                CDS(25227234, 25227412, "-", "12", ""),
                                CDS(25225614, 25225773, "-", "12", ""),
                                CDS(25215444, 25215560, "-", "12", "")
                            ]
                        )
                    ]
                ),
                Transcript(25233819, 25250929, "-", "12", "KRAS-203", {"id": "ENST00000556131"}, None,
                    [
                        Exon(25250764, 25250929, "-", "12", "ENST00000556131_e1"),
                        Exon(25245274, 25245395, "-", "12", "ENST00000556131_e2"),
                        Exon(25233819, 25235226, "-", "12", "ENST00000556131_e3")
                    ], [
                        Protein(25235209, 25245384, "-", "12", None, None, None,
                            [
                                CDS(25245274, 25245384, "-", "12", ""),
                                CDS(25235209, 25235226, "-", "12", "")
                            ]
                        )
                    ]
                )
            ])
        ]
        with open(self.tmp_ensembl_in_gtf, "w") as FH_gtf:
            FH_gtf.write("""#!genome-build GRCh38.p12
#!genome-version GRCh38
#!genome-date 2013-12
#!genome-build-accession NCBI:GCA_000001405.27
#!genebuild-last-updated 2018-07
6	havana	gene	54770583	54771134	.	+	.	gene_id \"ENSG00000220635\"; gene_version \"2\"; gene_name \"KRASP1\"; gene_source \"havana\"; gene_biotype \"processed_pseudogene\";
6	havana	transcript	54770583	54771134	.	+	.	gene_id \"ENSG00000220635\"; gene_version \"2\"; transcript_id \"ENST00000407852\"; transcript_version \"2\"; gene_name \"KRASP1\"; gene_source \"havana\"; gene_biotype \"processed_pseudogene\"; transcript_name \"KRASP1-201\"; transcript_source \"havana\"; transcript_biotype \"processed_pseudogene\"; tag \"basic\"; transcript_support_level \"NA\";
6	havana	exon	54770583	54771134	.	+	.	gene_id \"ENSG00000220635\"; gene_version \"2\"; transcript_id \"ENST00000407852\"; transcript_version \"2\"; exon_number \"1\"; gene_name \"KRASP1\"; gene_source \"havana\"; gene_biotype \"processed_pseudogene\"; transcript_name \"KRASP1-201\"; transcript_source \"havana\"; transcript_biotype \"processed_pseudogene\"; exon_id \"ENSE00001550689\"; exon_version \"2\"; tag \"basic\"; transcript_support_level \"NA\";
12	ensembl_havana	gene	25204789	25250936	.	-	.	gene_id \"ENSG00000133703\"; gene_version \"11\"; gene_name \"KRAS\"; gene_source \"ensembl_havana\"; gene_biotype \"protein_coding\";
12	ensembl_havana	transcript	25204789	25250931	.	-	.	gene_id \"ENSG00000133703\"; gene_version \"11\"; transcript_id \"ENST00000311936\"; transcript_version \"7\"; gene_name \"KRAS\"; gene_source \"ensembl_havana\"; gene_biotype \"protein_coding\"; transcript_name \"KRAS-202\"; transcript_source \"ensembl_havana\"; transcript_biotype \"protein_coding\"; tag \"CCDS\"; ccds_id \"CCDS8702\"; tag \"basic\"; transcript_support_level \"1\";
12	ensembl_havana	exon	25250751	25250931	.	-	.	gene_id \"ENSG00000133703\"; gene_version \"11\"; transcript_id \"ENST00000311936\"; transcript_version \"7\"; exon_number \"1\"; gene_name \"KRAS\"; gene_source \"ensembl_havana\"; gene_biotype \"protein_coding\"; transcript_name \"KRAS-202\"; transcript_source \"ensembl_havana\"; transcript_biotype \"protein_coding\"; tag \"CCDS\"; ccds_id \"CCDS8702\"; exon_id \"ENSE00001189804\"; exon_version \"4\"; tag \"basic\"; transcript_support_level \"1\";
12	ensembl_havana	exon	25245274	25245395	.	-	.	gene_id \"ENSG00000133703\"; gene_version \"11\"; transcript_id \"ENST00000311936\"; transcript_version \"7\"; exon_number \"2\"; gene_name \"KRAS\"; gene_source \"ensembl_havana\"; gene_biotype \"protein_coding\"; transcript_name \"KRAS-202\"; transcript_source \"ensembl_havana\"; transcript_biotype \"protein_coding\"; tag \"CCDS\"; ccds_id \"CCDS8702\"; exon_id \"ENSE00000936617\"; exon_version \"1\"; tag \"basic\"; transcript_support_level \"1\";
12	ensembl_havana	CDS	25245274	25245384	.	-	0	gene_id \"ENSG00000133703\"; gene_version \"11\"; transcript_id \"ENST00000311936\"; transcript_version \"7\"; exon_number \"2\"; gene_name \"KRAS\"; gene_source \"ensembl_havana\"; gene_biotype \"protein_coding\"; transcript_name \"KRAS-202\"; transcript_source \"ensembl_havana\"; transcript_biotype \"protein_coding\"; tag \"CCDS\"; ccds_id \"CCDS8702\"; protein_id \"ENSP00000308495\"; protein_version \"3\"; tag \"basic\"; transcript_support_level \"1\";
12	ensembl_havana	start_codon	25245382	25245384	.	-	0	gene_id \"ENSG00000133703\"; gene_version \"11\"; transcript_id \"ENST00000311936\"; transcript_version \"7\"; exon_number \"2\"; gene_name \"KRAS\"; gene_source \"ensembl_havana\"; gene_biotype \"protein_coding\"; transcript_name \"KRAS-202\"; transcript_source \"ensembl_havana\"; transcript_biotype \"protein_coding\"; tag \"CCDS\"; ccds_id \"CCDS8702\"; tag \"basic\"; transcript_support_level \"1\";
12	ensembl_havana	exon	25227234	25227412	.	-	.	gene_id \"ENSG00000133703\"; gene_version \"11\"; transcript_id \"ENST00000311936\"; transcript_version \"7\"; exon_number \"3\"; gene_name \"KRAS\"; gene_source \"ensembl_havana\"; gene_biotype \"protein_coding\"; transcript_name \"KRAS-202\"; transcript_source \"ensembl_havana\"; transcript_biotype \"protein_coding\"; tag \"CCDS\"; ccds_id \"CCDS8702\"; exon_id \"ENSE00001719809\"; exon_version \"1\"; tag \"basic\"; transcript_support_level \"1\";
12	ensembl_havana	CDS	25227234	25227412	.	-	0	gene_id \"ENSG00000133703\"; gene_version \"11\"; transcript_id \"ENST00000311936\"; transcript_version \"7\"; exon_number \"3\"; gene_name \"KRAS\"; gene_source \"ensembl_havana\"; gene_biotype \"protein_coding\"; transcript_name \"KRAS-202\"; transcript_source \"ensembl_havana\"; transcript_biotype \"protein_coding\"; tag \"CCDS\"; ccds_id \"CCDS8702\"; protein_id \"ENSP00000308495\"; protein_version \"3\"; tag \"basic\"; transcript_support_level \"1\";
12	ensembl_havana	exon	25225614	25225773	.	-	.	gene_id \"ENSG00000133703\"; gene_version \"11\"; transcript_id \"ENST00000311936\"; transcript_version \"7\"; exon_number \"4\"; gene_name \"KRAS\"; gene_source \"ensembl_havana\"; gene_biotype \"protein_coding\"; transcript_name \"KRAS-202\"; transcript_source \"ensembl_havana\"; transcript_biotype \"protein_coding\"; tag \"CCDS\"; ccds_id \"CCDS8702\"; exon_id \"ENSE00001644818\"; exon_version \"1\"; tag \"basic\"; transcript_support_level \"1\";
12	ensembl_havana	CDS	25225614	25225773	.	-	1	gene_id \"ENSG00000133703\"; gene_version \"11\"; transcript_id \"ENST00000311936\"; transcript_version \"7\"; exon_number \"4\"; gene_name \"KRAS\"; gene_source \"ensembl_havana\"; gene_biotype \"protein_coding\"; transcript_name \"KRAS-202\"; transcript_source \"ensembl_havana\"; transcript_biotype \"protein_coding\"; tag \"CCDS\"; ccds_id \"CCDS8702\"; protein_id \"ENSP00000308495\"; protein_version \"3\"; tag \"basic\"; transcript_support_level \"1\";
12	ensembl_havana	exon	25204789	25209911	.	-	.	gene_id \"ENSG00000133703\"; gene_version \"11\"; transcript_id \"ENST00000311936\"; transcript_version \"7\"; exon_number \"5\"; gene_name \"KRAS\"; gene_source \"ensembl_havana\"; gene_biotype \"protein_coding\"; transcript_name \"KRAS-202\"; transcript_source \"ensembl_havana\"; transcript_biotype \"protein_coding\"; tag \"CCDS\"; ccds_id \"CCDS8702\"; exon_id \"ENSE00002456976\"; exon_version \"1\"; tag \"basic\"; transcript_support_level \"1\";
12	ensembl_havana	CDS	25209798	25209911	.	-	0	gene_id \"ENSG00000133703\"; gene_version \"11\"; transcript_id \"ENST00000311936\"; transcript_version \"7\"; exon_number \"5\"; gene_name \"KRAS\"; gene_source \"ensembl_havana\"; gene_biotype \"protein_coding\"; transcript_name \"KRAS-202\"; transcript_source \"ensembl_havana\"; transcript_biotype \"protein_coding\"; tag \"CCDS\"; ccds_id \"CCDS8702\"; protein_id \"ENSP00000308495\"; protein_version \"3\"; tag \"basic\"; transcript_support_level \"1\";
12	ensembl_havana	stop_codon	25209795	25209797	.	-	0	gene_id \"ENSG00000133703\"; gene_version \"11\"; transcript_id \"ENST00000311936\"; transcript_version \"7\"; exon_number \"5\"; gene_name \"KRAS\"; gene_source \"ensembl_havana\"; gene_biotype \"protein_coding\"; transcript_name \"KRAS-202\"; transcript_source \"ensembl_havana\"; transcript_biotype \"protein_coding\"; tag \"CCDS\"; ccds_id \"CCDS8702\"; tag \"basic\"; transcript_support_level \"1\";
12	ensembl_havana	five_prime_utr	25250751	25250931	.	-	.	gene_id \"ENSG00000133703\"; gene_version \"11\"; transcript_id \"ENST00000311936\"; transcript_version \"7\"; gene_name \"KRAS\"; gene_source \"ensembl_havana\"; gene_biotype \"protein_coding\"; transcript_name \"KRAS-202\"; transcript_source \"ensembl_havana\"; transcript_biotype \"protein_coding\"; tag \"CCDS\"; ccds_id \"CCDS8702\"; tag \"basic\"; transcript_support_level \"1\";
12	ensembl_havana	five_prime_utr	25245385	25245395	.	-	.	gene_id \"ENSG00000133703\"; gene_version \"11\"; transcript_id \"ENST00000311936\"; transcript_version \"7\"; gene_name \"KRAS\"; gene_source \"ensembl_havana\"; gene_biotype \"protein_coding\"; transcript_name \"KRAS-202\"; transcript_source \"ensembl_havana\"; transcript_biotype \"protein_coding\"; tag \"CCDS\"; ccds_id \"CCDS8702\"; tag \"basic\"; transcript_support_level \"1\";
12	ensembl_havana	three_prime_utr	25204789	25209794	.	-	.	gene_id \"ENSG00000133703\"; gene_version \"11\"; transcript_id \"ENST00000311936\"; transcript_version \"7\"; gene_name \"KRAS\"; gene_source \"ensembl_havana\"; gene_biotype \"protein_coding\"; transcript_name \"KRAS-202\"; transcript_source \"ensembl_havana\"; transcript_biotype \"protein_coding\"; tag \"CCDS\"; ccds_id \"CCDS8702\"; tag \"basic\"; transcript_support_level \"1\";
12	havana	transcript	25209168	25250936	.	-	.	gene_id \"ENSG00000133703\"; gene_version \"11\"; transcript_id \"ENST00000557334\"; transcript_version \"5\"; gene_name \"KRAS\"; gene_source \"ensembl_havana\"; gene_biotype \"protein_coding\"; transcript_name \"KRAS-204\"; transcript_source \"havana\"; transcript_biotype \"protein_coding\"; tag \"basic\"; transcript_support_level \"5\";
12	havana	exon	25250751	25250936	.	-	.	gene_id \"ENSG00000133703\"; gene_version \"11\"; transcript_id \"ENST00000557334\"; transcript_version \"5\"; exon_number \"1\"; gene_name \"KRAS\"; gene_source \"ensembl_havana\"; gene_biotype \"protein_coding\"; transcript_name \"KRAS-204\"; transcript_source \"havana\"; transcript_biotype \"protein_coding\"; exon_id \"ENSE00002446502\"; exon_version \"1\"; tag \"basic\"; transcript_support_level \"5\";
12	havana	exon	25245274	25245395	.	-	.	gene_id \"ENSG00000133703\"; gene_version \"11\"; transcript_id \"ENST00000557334\"; transcript_version \"5\"; exon_number \"2\"; gene_name \"KRAS\"; gene_source \"ensembl_havana\"; gene_biotype \"protein_coding\"; transcript_name \"KRAS-204\"; transcript_source \"havana\"; transcript_biotype \"protein_coding\"; exon_id \"ENSE00000936617\"; exon_version \"1\"; tag \"basic\"; transcript_support_level \"5\";
12	havana	CDS	25245274	25245384	.	-	0	gene_id \"ENSG00000133703\"; gene_version \"11\"; transcript_id \"ENST00000557334\"; transcript_version \"5\"; exon_number \"2\"; gene_name \"KRAS\"; gene_source \"ensembl_havana\"; gene_biotype \"protein_coding\"; transcript_name \"KRAS-204\"; transcript_source \"havana\"; transcript_biotype \"protein_coding\"; protein_id \"ENSP00000452512\"; protein_version \"1\"; tag \"basic\"; transcript_support_level \"5\";
12	havana	start_codon	25245382	25245384	.	-	0	gene_id \"ENSG00000133703\"; gene_version \"11\"; transcript_id \"ENST00000557334\"; transcript_version \"5\"; exon_number \"2\"; gene_name \"KRAS\"; gene_source \"ensembl_havana\"; gene_biotype \"protein_coding\"; transcript_name \"KRAS-204\"; transcript_source \"havana\"; transcript_biotype \"protein_coding\"; tag \"basic\"; transcript_support_level \"5\";
12	havana	exon	25209168	25209911	.	-	.	gene_id \"ENSG00000133703\"; gene_version \"11\"; transcript_id \"ENST00000557334\"; transcript_version \"5\"; exon_number \"3\"; gene_name \"KRAS\"; gene_source \"ensembl_havana\"; gene_biotype \"protein_coding\"; transcript_name \"KRAS-204\"; transcript_source \"havana\"; transcript_biotype \"protein_coding\"; exon_id \"ENSE00002464674\"; exon_version \"1\"; tag \"basic\"; transcript_support_level \"5\";
12	havana	CDS	25209798	25209911	.	-	0	gene_id \"ENSG00000133703\"; gene_version \"11\"; transcript_id \"ENST00000557334\"; transcript_version \"5\"; exon_number \"3\"; gene_name \"KRAS\"; gene_source \"ensembl_havana\"; gene_biotype \"protein_coding\"; transcript_name \"KRAS-204\"; transcript_source \"havana\"; transcript_biotype \"protein_coding\"; protein_id \"ENSP00000452512\"; protein_version \"1\"; tag \"basic\"; transcript_support_level \"5\";
12	havana	stop_codon	25209795	25209797	.	-	0	gene_id \"ENSG00000133703\"; gene_version \"11\"; transcript_id \"ENST00000557334\"; transcript_version \"5\"; exon_number \"3\"; gene_name \"KRAS\"; gene_source \"ensembl_havana\"; gene_biotype \"protein_coding\"; transcript_name \"KRAS-204\"; transcript_source \"havana\"; transcript_biotype \"protein_coding\"; tag \"basic\"; transcript_support_level \"5\";
12	havana	five_prime_utr	25250751	25250936	.	-	.	gene_id \"ENSG00000133703\"; gene_version \"11\"; transcript_id \"ENST00000557334\"; transcript_version \"5\"; gene_name \"KRAS\"; gene_source \"ensembl_havana\"; gene_biotype \"protein_coding\"; transcript_name \"KRAS-204\"; transcript_source \"havana\"; transcript_biotype \"protein_coding\"; tag \"basic\"; transcript_support_level \"5\";
12	havana	five_prime_utr	25245385	25245395	.	-	.	gene_id \"ENSG00000133703\"; gene_version \"11\"; transcript_id \"ENST00000557334\"; transcript_version \"5\"; gene_name \"KRAS\"; gene_source \"ensembl_havana\"; gene_biotype \"protein_coding\"; transcript_name \"KRAS-204\"; transcript_source \"havana\"; transcript_biotype \"protein_coding\"; tag \"basic\"; transcript_support_level \"5\";
12	havana	three_prime_utr	25209168	25209794	.	-	.	gene_id \"ENSG00000133703\"; gene_version \"11\"; transcript_id \"ENST00000557334\"; transcript_version \"5\"; gene_name \"KRAS\"; gene_source \"ensembl_havana\"; gene_biotype \"protein_coding\"; transcript_name \"KRAS-204\"; transcript_source \"havana\"; transcript_biotype \"protein_coding\"; tag \"basic\"; transcript_support_level \"5\";
12	ensembl_havana	transcript	25209431	25250803	.	-	.	gene_id \"ENSG00000133703\"; gene_version \"11\"; transcript_id \"ENST00000256078\"; transcript_version \"8\"; gene_name \"KRAS\"; gene_source \"ensembl_havana\"; gene_biotype \"protein_coding\"; transcript_name \"KRAS-201\"; transcript_source \"ensembl_havana\"; transcript_biotype \"protein_coding\"; tag \"CCDS\"; ccds_id \"CCDS8703\"; tag \"basic\"; transcript_support_level \"1\";
12	ensembl_havana	exon	25250751	25250803	.	-	.	gene_id \"ENSG00000133703\"; gene_version \"11\"; transcript_id \"ENST00000256078\"; transcript_version \"8\"; exon_number \"1\"; gene_name \"KRAS\"; gene_source \"ensembl_havana\"; gene_biotype \"protein_coding\"; transcript_name \"KRAS-201\"; transcript_source \"ensembl_havana\"; transcript_biotype \"protein_coding\"; tag \"CCDS\"; ccds_id \"CCDS8703\"; exon_id \"ENSE00002513959\"; exon_version \"1\"; tag \"basic\"; transcript_support_level \"1\";
12	ensembl_havana	exon	25245274	25245395	.	-	.	gene_id \"ENSG00000133703\"; gene_version \"11\"; transcript_id \"ENST00000256078\"; transcript_version \"8\"; exon_number \"2\"; gene_name \"KRAS\"; gene_source \"ensembl_havana\"; gene_biotype \"protein_coding\"; transcript_name \"KRAS-201\"; transcript_source \"ensembl_havana\"; transcript_biotype \"protein_coding\"; tag \"CCDS\"; ccds_id \"CCDS8703\"; exon_id \"ENSE00000936617\"; exon_version \"1\"; tag \"basic\"; transcript_support_level \"1\";
12	ensembl_havana	CDS	25245274	25245384	.	-	0	gene_id \"ENSG00000133703\"; gene_version \"11\"; transcript_id \"ENST00000256078\"; transcript_version \"8\"; exon_number \"2\"; gene_name \"KRAS\"; gene_source \"ensembl_havana\"; gene_biotype \"protein_coding\"; transcript_name \"KRAS-201\"; transcript_source \"ensembl_havana\"; transcript_biotype \"protein_coding\"; tag \"CCDS\"; ccds_id \"CCDS8703\"; protein_id \"ENSP00000256078\"; protein_version \"4\"; tag \"basic\"; transcript_support_level \"1\";
12	ensembl_havana	start_codon	25245382	25245384	.	-	0	gene_id \"ENSG00000133703\"; gene_version \"11\"; transcript_id \"ENST00000256078\"; transcript_version \"8\"; exon_number \"2\"; gene_name \"KRAS\"; gene_source \"ensembl_havana\"; gene_biotype \"protein_coding\"; transcript_name \"KRAS-201\"; transcript_source \"ensembl_havana\"; transcript_biotype \"protein_coding\"; tag \"CCDS\"; ccds_id \"CCDS8703\"; tag \"basic\"; transcript_support_level \"1\";
12	ensembl_havana	exon	25227234	25227412	.	-	.	gene_id \"ENSG00000133703\"; gene_version \"11\"; transcript_id \"ENST00000256078\"; transcript_version \"8\"; exon_number \"3\"; gene_name \"KRAS\"; gene_source \"ensembl_havana\"; gene_biotype \"protein_coding\"; transcript_name \"KRAS-201\"; transcript_source \"ensembl_havana\"; transcript_biotype \"protein_coding\"; tag \"CCDS\"; ccds_id \"CCDS8703\"; exon_id \"ENSE00001719809\"; exon_version \"1\"; tag \"basic\"; transcript_support_level \"1\";
12	ensembl_havana	CDS	25227234	25227412	.	-	0	gene_id \"ENSG00000133703\"; gene_version \"11\"; transcript_id \"ENST00000256078\"; transcript_version \"8\"; exon_number \"3\"; gene_name \"KRAS\"; gene_source \"ensembl_havana\"; gene_biotype \"protein_coding\"; transcript_name \"KRAS-201\"; transcript_source \"ensembl_havana\"; transcript_biotype \"protein_coding\"; tag \"CCDS\"; ccds_id \"CCDS8703\"; protein_id \"ENSP00000256078\"; protein_version \"4\"; tag \"basic\"; transcript_support_level \"1\";
12	ensembl_havana	exon	25225614	25225773	.	-	.	gene_id \"ENSG00000133703\"; gene_version \"11\"; transcript_id \"ENST00000256078\"; transcript_version \"8\"; exon_number \"4\"; gene_name \"KRAS\"; gene_source \"ensembl_havana\"; gene_biotype \"protein_coding\"; transcript_name \"KRAS-201\"; transcript_source \"ensembl_havana\"; transcript_biotype \"protein_coding\"; tag \"CCDS\"; ccds_id \"CCDS8703\"; exon_id \"ENSE00001644818\"; exon_version \"1\"; tag \"basic\"; transcript_support_level \"1\";
12	ensembl_havana	CDS	25225614	25225773	.	-	1	gene_id \"ENSG00000133703\"; gene_version \"11\"; transcript_id \"ENST00000256078\"; transcript_version \"8\"; exon_number \"4\"; gene_name \"KRAS\"; gene_source \"ensembl_havana\"; gene_biotype \"protein_coding\"; transcript_name \"KRAS-201\"; transcript_source \"ensembl_havana\"; transcript_biotype \"protein_coding\"; tag \"CCDS\"; ccds_id \"CCDS8703\"; protein_id \"ENSP00000256078\"; protein_version \"4\"; tag \"basic\"; transcript_support_level \"1\";
12	ensembl_havana	exon	25215437	25215560	.	-	.	gene_id \"ENSG00000133703\"; gene_version \"11\"; transcript_id \"ENST00000256078\"; transcript_version \"8\"; exon_number \"5\"; gene_name \"KRAS\"; gene_source \"ensembl_havana\"; gene_biotype \"protein_coding\"; transcript_name \"KRAS-201\"; transcript_source \"ensembl_havana\"; transcript_biotype \"protein_coding\"; tag \"CCDS\"; ccds_id \"CCDS8703\"; exon_id \"ENSE00001189807\"; exon_version \"5\"; tag \"basic\"; transcript_support_level \"1\";
12	ensembl_havana	CDS	25215444	25215560	.	-	0	gene_id \"ENSG00000133703\"; gene_version \"11\"; transcript_id \"ENST00000256078\"; transcript_version \"8\"; exon_number \"5\"; gene_name \"KRAS\"; gene_source \"ensembl_havana\"; gene_biotype \"protein_coding\"; transcript_name \"KRAS-201\"; transcript_source \"ensembl_havana\"; transcript_biotype \"protein_coding\"; tag \"CCDS\"; ccds_id \"CCDS8703\"; protein_id \"ENSP00000256078\"; protein_version \"4\"; tag \"basic\"; transcript_support_level \"1\";
12	ensembl_havana	stop_codon	25215441	25215443	.	-	0	gene_id \"ENSG00000133703\"; gene_version \"11\"; transcript_id \"ENST00000256078\"; transcript_version \"8\"; exon_number \"5\"; gene_name \"KRAS\"; gene_source \"ensembl_havana\"; gene_biotype \"protein_coding\"; transcript_name \"KRAS-201\"; transcript_source \"ensembl_havana\"; transcript_biotype \"protein_coding\"; tag \"CCDS\"; ccds_id \"CCDS8703\"; tag \"basic\"; transcript_support_level \"1\";
12	ensembl_havana	exon	25209431	25209911	.	-	.	gene_id \"ENSG00000133703\"; gene_version \"11\"; transcript_id \"ENST00000256078\"; transcript_version \"8\"; exon_number \"6\"; gene_name \"KRAS\"; gene_source \"ensembl_havana\"; gene_biotype \"protein_coding\"; transcript_name \"KRAS-201\"; transcript_source \"ensembl_havana\"; transcript_biotype \"protein_coding\"; tag \"CCDS\"; ccds_id \"CCDS8703\"; exon_id \"ENSE00002477035\"; exon_version \"1\"; tag \"basic\"; transcript_support_level \"1\";
12	ensembl_havana	five_prime_utr	25250751	25250803	.	-	.	gene_id \"ENSG00000133703\"; gene_version \"11\"; transcript_id \"ENST00000256078\"; transcript_version \"8\"; gene_name \"KRAS\"; gene_source \"ensembl_havana\"; gene_biotype \"protein_coding\"; transcript_name \"KRAS-201\"; transcript_source \"ensembl_havana\"; transcript_biotype \"protein_coding\"; tag \"CCDS\"; ccds_id \"CCDS8703\"; tag \"basic\"; transcript_support_level \"1\";
12	ensembl_havana	five_prime_utr	25245385	25245395	.	-	.	gene_id \"ENSG00000133703\"; gene_version \"11\"; transcript_id \"ENST00000256078\"; transcript_version \"8\"; gene_name \"KRAS\"; gene_source \"ensembl_havana\"; gene_biotype \"protein_coding\"; transcript_name \"KRAS-201\"; transcript_source \"ensembl_havana\"; transcript_biotype \"protein_coding\"; tag \"CCDS\"; ccds_id \"CCDS8703\"; tag \"basic\"; transcript_support_level \"1\";
12	ensembl_havana	three_prime_utr	25215437	25215440	.	-	.	gene_id \"ENSG00000133703\"; gene_version \"11\"; transcript_id \"ENST00000256078\"; transcript_version \"8\"; gene_name \"KRAS\"; gene_source \"ensembl_havana\"; gene_biotype \"protein_coding\"; transcript_name \"KRAS-201\"; transcript_source \"ensembl_havana\"; transcript_biotype \"protein_coding\"; tag \"CCDS\"; ccds_id \"CCDS8703\"; tag \"basic\"; transcript_support_level \"1\";
12	ensembl_havana	three_prime_utr	25209431	25209911	.	-	.	gene_id \"ENSG00000133703\"; gene_version \"11\"; transcript_id \"ENST00000256078\"; transcript_version \"8\"; gene_name \"KRAS\"; gene_source \"ensembl_havana\"; gene_biotype \"protein_coding\"; transcript_name \"KRAS-201\"; transcript_source \"ensembl_havana\"; transcript_biotype \"protein_coding\"; tag \"CCDS\"; ccds_id \"CCDS8703\"; tag \"basic\"; transcript_support_level \"1\";
12	havana	transcript	25233819	25250929	.	-	.	gene_id \"ENSG00000133703\"; gene_version \"11\"; transcript_id \"ENST00000556131\"; transcript_version \"1\"; gene_name \"KRAS\"; gene_source \"ensembl_havana\"; gene_biotype \"protein_coding\"; transcript_name \"KRAS-203\"; transcript_source \"havana\"; transcript_biotype \"protein_coding\"; tag \"basic\"; transcript_support_level \"1\";
12	havana	exon	25250764	25250929	.	-	.	gene_id \"ENSG00000133703\"; gene_version \"11\"; transcript_id \"ENST00000556131\"; transcript_version \"1\"; exon_number \"1\"; gene_name \"KRAS\"; gene_source \"ensembl_havana\"; gene_biotype \"protein_coding\"; transcript_name \"KRAS-203\"; transcript_source \"havana\"; transcript_biotype \"protein_coding\"; exon_id \"ENSE00002530521\"; exon_version \"1\"; tag \"basic\"; transcript_support_level \"1\";
12	havana	exon	25245274	25245395	.	-	.	gene_id \"ENSG00000133703\"; gene_version \"11\"; transcript_id \"ENST00000556131\"; transcript_version \"1\"; exon_number \"2\"; gene_name \"KRAS\"; gene_source \"ensembl_havana\"; gene_biotype \"protein_coding\"; transcript_name \"KRAS-203\"; transcript_source \"havana\"; transcript_biotype \"protein_coding\"; exon_id \"ENSE00000936617\"; exon_version \"1\"; tag \"basic\"; transcript_support_level \"1\";
12	havana	CDS	25245274	25245384	.	-	0	gene_id \"ENSG00000133703\"; gene_version \"11\"; transcript_id \"ENST00000556131\"; transcript_version \"1\"; exon_number \"2\"; gene_name \"KRAS\"; gene_source \"ensembl_havana\"; gene_biotype \"protein_coding\"; transcript_name \"KRAS-203\"; transcript_source \"havana\"; transcript_biotype \"protein_coding\"; protein_id \"ENSP00000451856\"; protein_version \"1\"; tag \"basic\"; transcript_support_level \"1\";
12	havana	start_codon	25245382	25245384	.	-	0	gene_id \"ENSG00000133703\"; gene_version \"11\"; transcript_id \"ENST00000556131\"; transcript_version \"1\"; exon_number \"2\"; gene_name \"KRAS\"; gene_source \"ensembl_havana\"; gene_biotype \"protein_coding\"; transcript_name \"KRAS-203\"; transcript_source \"havana\"; transcript_biotype \"protein_coding\"; tag \"basic\"; transcript_support_level \"1\";
12	havana	exon	25233819	25235226	.	-	.	gene_id \"ENSG00000133703\"; gene_version \"11\"; transcript_id \"ENST00000556131\"; transcript_version \"1\"; exon_number \"3\"; gene_name \"KRAS\"; gene_source \"ensembl_havana\"; gene_biotype \"protein_coding\"; transcript_name \"KRAS-203\"; transcript_source \"havana\"; transcript_biotype \"protein_coding\"; exon_id \"ENSE00002478081\"; exon_version \"1\"; tag \"basic\"; transcript_support_level \"1\";
12	havana	CDS	25235209	25235226	.	-	0	gene_id \"ENSG00000133703\"; gene_version \"11\"; transcript_id \"ENST00000556131\"; transcript_version \"1\"; exon_number \"3\"; gene_name \"KRAS\"; gene_source \"ensembl_havana\"; gene_biotype \"protein_coding\"; transcript_name \"KRAS-203\"; transcript_source \"havana\"; transcript_biotype \"protein_coding\"; protein_id \"ENSP00000451856\"; protein_version \"1\"; tag \"basic\"; transcript_support_level \"1\";
12	havana	stop_codon	25235206	25235208	.	-	0	gene_id \"ENSG00000133703\"; gene_version \"11\"; transcript_id \"ENST00000556131\"; transcript_version \"1\"; exon_number \"3\"; gene_name \"KRAS\"; gene_source \"ensembl_havana\"; gene_biotype \"protein_coding\"; transcript_name \"KRAS-203\"; transcript_source \"havana\"; transcript_biotype \"protein_coding\"; tag \"basic\"; transcript_support_level \"1\";
12	havana	five_prime_utr	25250764	25250929	.	-	.	gene_id \"ENSG00000133703\"; gene_version \"11\"; transcript_id \"ENST00000556131\"; transcript_version \"1\"; gene_name \"KRAS\"; gene_source \"ensembl_havana\"; gene_biotype \"protein_coding\"; transcript_name \"KRAS-203\"; transcript_source \"havana\"; transcript_biotype \"protein_coding\"; tag \"basic\"; transcript_support_level \"1\";
12	havana	five_prime_utr	25245385	25245395	.	-	.	gene_id \"ENSG00000133703\"; gene_version \"11\"; transcript_id \"ENST00000556131\"; transcript_version \"1\"; gene_name \"KRAS\"; gene_source \"ensembl_havana\"; gene_biotype \"protein_coding\"; transcript_name \"KRAS-203\"; transcript_source \"havana\"; transcript_biotype \"protein_coding\"; tag \"basic; other\"; transcript_support_level \"1\";
12	havana	three_prime_utr	25233819	25235205	.	-	.	gene_id \"ENSG00000133703\"; gene_version \"11\"; transcript_id \"ENST00000556131\"; transcript_version \"1\"; gene_name \"KRAS\"; gene_source \"ensembl_havana\"; gene_biotype \"protein_coding\"; transcript_name \"KRAS-203\"; transcript_source \"havana\"; transcript_biotype \"protein_coding\"; tag \"basic;\"; transcript_support_level \"1\";""")


        # NCBI GTF
        self.ncbi_expected = [
            Gene(54635272, 54640529, "+", "6", "KRASP1", {"id": "3844"}, None, [
                Transcript(54635272, 54640529, "+", "6", "gene14201", {"id": "gene14201"}, None, [
                    Exon(54635272, 54640529, "+", "6", "gene14201_e1")
                ])
            ]),
            Gene(25357723, 25403865, "-", "12", "KRAS", {"id": "3845"}, None, [
                Transcript(25357723, 25403865, "-", "12", "rna36549", {"id": "rna36549"}, None,
                    [
                        Exon(25357723, 25362845, "-", "12", "rna36549_e1"),
                        Exon(25378548, 25378707, "-", "12", "rna36549_e2"),
                        Exon(25380168, 25380346, "-", "12", "rna36549_e3"),
                        Exon(25398208, 25398329, "-", "12", "rna36549_e5"),
                        Exon(25403685, 25403865, "-", "12", "rna36549_e6")
                    ], [
                        Protein(25362729, 25398318, "-", "12", None, None, None,
                            [
                                CDS(25362729, 25362845, "-", "12", ""),
                                CDS(25378548, 25378707, "-", "12", ""),
                                CDS(25380168, 25380346, "-", "12", ""),
                                CDS(25398208, 25398318, "-", "12", "")
                            ]
                        )
                    ]
                ),
                Transcript(25357723, 25403865, "-", "12", "rna36550", {"id": "rna36550"}, None,
                    [
                        Exon(25357723, 25362845, "-", "12", "rna36550_e1"),
                        Exon(25368371, 25368494, "-", "12", "rna36550_e2"),
                        Exon(25378548, 25378707, "-", "12", "rna36550_e3"),
                        Exon(25380168, 25380346, "-", "12", "rna36550_e4"),
                        Exon(25398208, 25398329, "-", "12", "rna36550_e5"),
                        Exon(25403685, 25403865, "-", "12", "rna36550_e6")
                    ], [
                        Protein(25368375, 25398318, "-", "12", None, None, None,
                            [
                                CDS(25368375, 25368494, "-", "12", ""),
                                CDS(25378548, 25378707, "-", "12", ""),
                                CDS(25380168, 25380346, "-", "12", ""),
                                CDS(25398208, 25398318, "-", "12", ""),
                            ]
                        )
                    ]
                )
            ]),
            Gene(2527306, 2529079, "+", "X", "CD99P1", {"id": "401577"}, None, [
                Transcript(2527306, 2529079, "+", "X", "rna58916", {"id": "rna58916"}, None, [
                    Exon(2527306, 2527522, "+", "X", "rna58916_e1"),
                    Exon(2529037, 2529079, "+", "X", "rna58916_e2")
                ])
            ]),
            Gene(2477306, 2479079, "+", "Y", "CD99P1", {"id": "401577"}, None, [
                Transcript(2477306, 2479079, "+", "Y", "rna61353", {"id": "rna61353"}, None, [
                    Exon(2477306, 2477522, "+", "Y", "rna61353_e1"),
                    Exon(2479037, 2479079, "+", "Y", "rna61353_e2")
                ])
            ])
        ]
        with open(self.tmp_ncbi_in_gtf, "w") as FH_gtf:
            FH_gtf.write("""6	Curated Genomic	exon	54635272	54640529	.	+	.	transcript_id \"gene14201\"; gene_id \"3844\"; gene_name \"KRASP1\";
12	BestRefSeq	exon	25357723	25362845	.	-	.	transcript_id \"rna36549\"; gene_id \"3845\"; gene_name \"KRAS\";
12	BestRefSeq	exon	25378548	25378707	.	-	.	transcript_id \"rna36549\"; gene_id \"3845\"; gene_name \"KRAS\";
12	BestRefSeq	exon	25380168	25380346	.	-	.	transcript_id \"rna36549\"; gene_id \"3845\"; gene_name \"KRAS\";
12	BestRefSeq	exon	25398208	25398329	.	-	.	transcript_id \"rna36549\"; gene_id \"3845\"; gene_name \"KRAS\";
12	BestRefSeq	exon	25403685	25403865	.	-	.	transcript_id \"rna36549\"; gene_id \"3845\"; gene_name \"KRAS\";
12	BestRefSeq	CDS	25362729	25362845	.	-	0	transcript_id \"rna36549\"; gene_id \"3845\"; gene_name \"KRAS\";
12	BestRefSeq	CDS	25378548	25378707	.	-	1	transcript_id \"rna36549\"; gene_id \"3845\"; gene_name \"KRAS\";
12	BestRefSeq	CDS	25380168	25380346	.	-	0	transcript_id \"rna36549\"; gene_id \"3845\"; gene_name \"KRAS\";
12	BestRefSeq	CDS	25398208	25398318	.	-	0	transcript_id \"rna36549\"; gene_id \"3845\"; gene_name \"KRAS\";
12	BestRefSeq	exon	25357723	25362845	.	-	.	transcript_id \"rna36550\"; gene_id \"3845\"; gene_name \"KRAS\";
12	BestRefSeq	exon	25368371	25368494	.	-	.	transcript_id \"rna36550\"; gene_id \"3845\"; gene_name \"KRAS\";
12	BestRefSeq	exon	25378548	25378707	.	-	.	transcript_id \"rna36550\"; gene_id \"3845\"; gene_name \"KRAS\";
12	BestRefSeq	exon	25380168	25380346	.	-	.	transcript_id \"rna36550\"; gene_id \"3845\"; gene_name \"KRAS\";
12	BestRefSeq	exon	25398208	25398329	.	-	.	transcript_id \"rna36550\"; gene_id \"3845\"; gene_name \"KRAS\";
12	BestRefSeq	exon	25403685	25403865	.	-	.	transcript_id \"rna36550\"; gene_id \"3845\"; gene_name \"KRAS\";
12	BestRefSeq	CDS	25368375	25368494	.	-	0	transcript_id \"rna36550\"; gene_id \"3845\"; gene_name \"KRAS\";
12	BestRefSeq	CDS	25378548	25378707	.	-	1	transcript_id \"rna36550\"; gene_id \"3845\"; gene_name \"KRAS\";
12	BestRefSeq	CDS	25380168	25380346	.	-	0	transcript_id \"rna36550\"; gene_id \"3845\"; gene_name \"KRAS\";
12	BestRefSeq	CDS	25398208	25398318	.	-	0	transcript_id \"rna36550\"; gene_id \"3845\"; gene_name \"KRAS\";
X	BestRefSeq	exon	2527306	2527522	.	+	.	transcript_id \"rna58916\"; gene_id \"401577\"; gene_name \"CD99P1\";
X	BestRefSeq	exon	2529037	2529079	.	+	.	transcript_id \"rna58916\"; gene_id \"401577\"; gene_name \"CD99P1\";
Y	BestRefSeq	exon	2477306	2477522	.	+	.	transcript_id \"rna61353\"; gene_id \"401577\"; gene_name \"CD99P1\";
Y	BestRefSeq	exon	2479037	2479079	.	+	.	transcript_id \"rna61353\"; gene_id \"401577\"; gene_name \"CD99P1\";""")


    def testLoadModelEnsembl(self):
        ensembl_genes = loadModel(self.tmp_ensembl_in_gtf, "genes")
        self.assertEqual(
            toBracketTree(self.ensembl_expected),
            toBracketTree(ensembl_genes)
        )

    def testLoadModelNCBI(self):
        ncbi_genes = loadModel(self.tmp_ncbi_in_gtf, "genes")
        self.assertEqual(
            toBracketTree(self.ncbi_expected),
            toBracketTree(ncbi_genes)
        )

    def testWriteEnsembl(self):
        # Read expected
        in_records = None
        with GTFIO(self.tmp_ensembl_in_gtf) as FH_in:
            in_records = FH_in.read()
        # Write
        with GTFIO(self.tmp_ensembl_out_gtf, "w") as FH_out:
            for record in in_records:
                FH_out.write(record)
        # Read writted
        out_records = None
        with GTFIO(self.tmp_ensembl_in_gtf) as FH_in:
            out_records = FH_in.read()
        # Compare expected and observed
        for exp_rec, obs_rec in zip(in_records, out_records):
            self.assertEqual(
                recordToStr(exp_rec),
                recordToStr(obs_rec)
            )

    def testWriteNCBI(self):
        # Read expected
        in_records = None
        with GTFIO(self.tmp_ncbi_in_gtf) as FH_in:
            in_records = FH_in.read()
        # Write
        with GTFIO(self.tmp_ncbi_out_gtf, "w") as FH_out:
            for record in in_records:
                FH_out.write(record)
        # Read writted
        out_records = None
        with GTFIO(self.tmp_ncbi_in_gtf) as FH_in:
            out_records = FH_in.read()
        # Compare expected and observed
        for exp_rec, obs_rec in zip(in_records, out_records):
            self.assertEqual(
                recordToStr(exp_rec),
                recordToStr(obs_rec)
            )

    def tearDown(self):
        # Clean temporary files
        for curr_file in [self.tmp_ensembl_in_gtf, self.tmp_ensembl_out_gtf, self.tmp_ncbi_in_gtf, self.tmp_ncbi_out_gtf]:
            if os.path.exists(curr_file):
                os.remove(curr_file)


def recordToStr(record):
    attributes = []
    for key, val in sorted(record.annot.items()):
        if key not in ["source", "feature", "score", "frame"]:
            attributes.append('{} "{}"'.format(key, val))
    return "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(
        ("." if record.reference is None else record.reference.name),
        ("." if "source" not in record.annot else record.annot["source"]),
        ("." if "feature" not in record.annot else record.annot["feature"]),
        record.start,
        record.end,
        ("." if "score" not in record.annot else record.annot["score"]),
        ("." if record.strand is None else record.strand),
        ("." if "frame" not in record.annot else record.annot["frame"]),
        "; ".join(attributes)
    )

def toBracketTree(node):
    node_btree = ""
    if issubclass(node.__class__, RegionTree):
        node_btree = "{}:{}".format(node.__class__.__name__, node.getCoordinatesStr())
        if len(node.children) > 0:
            children_btree = []
            for child in node.children:
                children_btree.append(toBracketTree(child))
            node_btree += "({})".format(",".join(children_btree))
        if issubclass(node.__class__, Transcript) and len(node.proteins) > 0:
            prot_btree = []
            for prot in node.proteins:
                prot_btree.append(toBracketTree(prot))
            node_btree += "@({})".format(",".join(prot_btree))
    else:
        children_btree = []
        for child in node:
            children_btree.append(toBracketTree(child))
        node_btree = "({})".format(",".join(children_btree))
    return node_btree


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    unittest.main()
    ###### tester UTRs
    ###### tester prot inclus dans autre
