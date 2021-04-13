#!/usr/bin/env python3

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2020 IUCT-O'
__license__ = 'GNU General Public License'
__version__ = '1.1.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import os
import sys
import unittest

TEST_DIR = os.path.dirname(os.path.abspath(__file__))
PACKAGE_DIR = os.path.dirname(TEST_DIR)
sys.path.append(PACKAGE_DIR)

from anacore.hgvs import HGVSProtChange


########################################################################
#
# FUNCTIONS
#
########################################################################
class TestHGVSProtChange(unittest.TestCase):
    def testInsCouldBeIdentical(self):
        dataset = {
            "Ala3_Tyr4insAla": {
                "Ala3Dup": True
            },
            "T16_I17insGTTT": {
                "L16V": False,
                "G13_T16dup": True,
                "G13_L16dup": False,
                "G13_L16[1]": False,
                "G13_L16[2]": False
            },
            "L16_I17insGTTL": {
                "L16V": False,
                "G13_T16dup": False,
                "G13_L16dup": True,
                "G13_L16[1]": True,
                "G13_L16[2]": False
            },
            "L16_I17insGTTLGTTL": {
                "L16V": False,
                "G13_T16dup": False,
                "G13_L16dup": False,
                "G13_L16[1]": False,
                "G13_L16[2]": True
            },
            "L16_I17insGTTLGTLL": {
                "L16V": False,
                "G13_T16dup": False,
                "G13_L16dup": False,
                "G13_L16[1]": False,
                "G13_L16[2]": False
            }
        }
        for hgvs_ins, expected_by_cmp in dataset.items():
            hgvs_ins_obj = HGVSProtChange.fromStr(hgvs_ins)
            for hgvs_repeat, expected in expected_by_cmp.items():
                hgvs_repeat_obj = HGVSProtChange.fromStr(hgvs_repeat)
                self.assertEqual(HGVSProtChange.insCouldBeIdentical(hgvs_ins_obj, hgvs_repeat_obj), expected)

    def testIsInFrameIns(self):
        dataset = {
            # Special cases
            "?": False,
            "0": False,
            # Substitution
            "Val600=": False,
            "Val600Glu": False,
            "(G56A^S^C)": False,  # Uncertain
            "(Gly56Ala^Ser^Cys)": False,  # Uncertain
            "W24=/C": False,  # Mosaic
            "Trp24=/Cys": False,  # Mosaic
            # Deletion
            "(Val7del)": False,
            "Trp4del": False,
            "Lys23_Val25del": False,
            "(Pro458_Gly460del)": False,
            "Gly2_Met46del": False,
            "Trp26*": False,
            "Pro551_Glu554delProMetTyrGlu": False,
            "Trp557_Lys558DelTrpLys": False,
            "Asp842_His845DelAspIleMetHis": False,
            "Asp842DelAsp": False,
            # DelIns
            "Cys28delinsTrpVal": False,
            # Duplication
            "Ala3dup": True,
            "Ala3_Ser5dup": True,
            "A767_V769dupASV": True,
            # Insertion
            "His4_Gln5insAla": True,
            "Lys2_Gly3insGlnSerLys": True,
            "(Met3_His4insGlyTer)": True,
            "Arg78_Gly79ins23": True,
            "Gln746_Lys747ins*63": True,
            "(Ser332_Ser333ins(1))": True,
            "(Ser332_Ser333insX)": True,
            "(Val582_Asn583ins(5))": True,
            "(Val582_Asn583insXaaXaaXaaXaaXaa)": True,
            # Frameshift
            "Arg97ProfsTer23": False,
            "Arg97fs": False,
            "(Tyr4*)": False,
            "Ile327Argfs*?": False,
            "Gln151Thrfs*9": False,
            "*757Leu": False,
            "Val600Ter": False,
            # Extension
            "Met1ext-5": False,
            "Ter110GlnextTer17": False,
            "(Ter315TyrextAsnLysGlyThrTer)": False,
            "*315TyrextAsnLysGlyThr*": False,
            "Ter327Argext*?": False,
            # Reapeated sequences
            "Ala2[10]": True,
            "Ala2[1]": True,
            "Ala2_Pro5[10]": True,
            "Ala2_Pro5[1]": True,
        }
        for hgvs, expected in dataset.items():
            hgvs_obj = HGVSProtChange.fromStr(hgvs)
            self.assertEqual(hgvs_obj.isInFrameIns(), expected)

    def testParse(self):
        expected = HGVSProtChange("Ter", 757, None, None, None, None, ["Leu"], False)
        observed = HGVSProtChange.fromStr("*757L")
        self.assertEqual(observed, expected)
        observed = HGVSProtChange.fromStr("*757Leu")
        self.assertEqual(observed, expected)
        expected = HGVSProtChange("Val", 600, None, None, None, None, ["="], False)
        observed = HGVSProtChange.fromStr("V600=")
        self.assertEqual(observed, expected)
        observed = HGVSProtChange.fromStr("Val600=")
        self.assertEqual(observed, expected)
        expected = HGVSProtChange("Val", 600, None, None, None, None, ["Glu"], False)
        observed = HGVSProtChange.fromStr("V600E")
        self.assertEqual(observed, expected)
        observed = HGVSProtChange.fromStr("Val600Glu")
        self.assertEqual(observed, expected)
        expected = HGVSProtChange("Val", 600, None, None, None, None, ["Ter"], False)
        observed = HGVSProtChange.fromStr("V600*")
        self.assertEqual(observed, expected)
        observed = HGVSProtChange.fromStr("Val600Ter")
        self.assertEqual(observed, expected)
        expected = HGVSProtChange("Val", 600, None, None, None, "del", [], False)
        observed = HGVSProtChange.fromStr("V600del")
        self.assertEqual(observed, expected)
        observed = HGVSProtChange.fromStr("Val600del")
        self.assertEqual(observed, expected)
        expected = HGVSProtChange("Val", 600, "Ala", 603, None, "del", [], False)
        observed = HGVSProtChange.fromStr("V600_A603del")
        self.assertEqual(observed, expected)
        observed = HGVSProtChange.fromStr("Val600_Ala603del")
        self.assertEqual(observed, expected)
        expected = HGVSProtChange("Val", 600, "Ala", 603, None, "dup", [], False)
        observed = HGVSProtChange.fromStr("V600_A603dup")
        self.assertEqual(observed, expected)
        observed = HGVSProtChange.fromStr("Val600_Ala603dup")
        self.assertEqual(observed, expected)
        expected = HGVSProtChange("Val", 600, "Val", 601, None, "ins", ["(15)"], False)
        observed = HGVSProtChange.fromStr("Val600_Val601ins(15)")
        self.assertEqual(observed, expected)
        expected = HGVSProtChange("Val", 600, "Val", 601, None, "ins", ["15"], False)
        observed = HGVSProtChange.fromStr("Val600_Val601ins15")
        self.assertEqual(observed, expected)
        expected = HGVSProtChange("Cys", 28, None, None, None, "delins", ["Trp", "Val"], False)
        observed = HGVSProtChange.fromStr("C28delinsWV")
        self.assertEqual(observed, expected)
        observed = HGVSProtChange.fromStr("Cys28delinsTrpVal")
        self.assertEqual(observed, expected)
        expected = HGVSProtChange("Cys", 28, "Lys", 29, None, "delins", ["Trp", "Val", "Ter"], False)
        observed = HGVSProtChange.fromStr("C28_K29delinsWV*")
        self.assertEqual(observed, expected)
        observed = HGVSProtChange.fromStr("Cys28_Lys29delinsTrpValTer")
        self.assertEqual(observed, expected)
        # Special cases
        expected = HGVSProtChange(None, None, None, None, None, None, None, False)
        observed = HGVSProtChange.fromStr("?")  # Unknown
        self.assertEqual(observed, expected)
        expected = HGVSProtChange(None, "0", None, None, None, None, None, False)
        observed = HGVSProtChange.fromStr("0")  # No protein
        self.assertEqual(observed, expected)
        # Substitution
        expected = HGVSProtChange("Trp", 24, None, None, None, None, ["Cys"], False)
        observed = HGVSProtChange.fromStr("W24C")  # Missense
        self.assertEqual(observed, expected)
        observed = HGVSProtChange.fromStr("Trp24Cys")  # Missense
        self.assertEqual(observed, expected)
        expected = HGVSProtChange("Trp", 24, None, None, None, None, ["Cys"], True)
        observed = HGVSProtChange.fromStr("(Trp24Cys)")  # Missense
        self.assertEqual(observed, expected)
        expected = HGVSProtChange("Trp", 24, None, None, None, None, ["Ter"], False)
        observed = HGVSProtChange.fromStr("W24*")  # Nonsense
        self.assertEqual(observed, expected)
        observed = HGVSProtChange.fromStr("Trp24*")  # Nonsense
        self.assertEqual(observed, expected)
        observed = HGVSProtChange.fromStr("Trp24Ter")  # Nonsense
        self.assertEqual(observed, expected)
        expected = HGVSProtChange("Cys", 188, None, None, None, None, ["="], False)
        observed = HGVSProtChange.fromStr("C188=")  # No change
        self.assertEqual(observed, expected)
        observed = HGVSProtChange.fromStr("Cys188=")  # No change
        self.assertEqual(observed, expected)
        expected = HGVSProtChange("Gly", 56, None, None, None, None, ["Ala", "^", "Ser", "^", "Cys"], True)
        observed = HGVSProtChange.fromStr("(G56A^S^C)")  # Uncertain
        self.assertEqual(observed, expected)
        observed = HGVSProtChange.fromStr("(Gly56Ala^Ser^Cys)")  # Uncertain
        self.assertEqual(observed, expected)
        expected = HGVSProtChange("Trp", 24, None, None, None, None, ["=", "/", "Cys"], False)
        observed = HGVSProtChange.fromStr("W24=/C")  # Mosaic
        self.assertEqual(observed, expected)
        observed = HGVSProtChange.fromStr("Trp24=/Cys")  # Mosaic
        self.assertEqual(observed, expected)
        # Deletion
        expected = HGVSProtChange("Val", 7, None, None, None, "del", [], False)
        observed = HGVSProtChange.fromStr("V7del")
        self.assertEqual(observed, expected)
        observed = HGVSProtChange.fromStr("Val7del")
        self.assertEqual(observed, expected)
        expected = HGVSProtChange("Val", 7, None, None, None, "del", [], True)
        observed = HGVSProtChange.fromStr("(V7del)")
        self.assertEqual(observed, expected)
        observed = HGVSProtChange.fromStr("(Val7del)")
        self.assertEqual(observed, expected)
        expected = HGVSProtChange("Trp", 4, None, None, None, "del", [], False)
        observed = HGVSProtChange.fromStr("Trp4del")
        self.assertEqual(observed, expected)
        expected = HGVSProtChange("Lys", 23, "Val", 25, None, "del", [], False)
        observed = HGVSProtChange.fromStr("K23_V25del")
        self.assertEqual(observed, expected)
        observed = HGVSProtChange.fromStr("Lys23_Val25del")
        self.assertEqual(observed, expected)
        expected = HGVSProtChange("Pro", 458, "Gly", 460, None, "del", [], True)
        observed = HGVSProtChange.fromStr("(P458_G460del)")
        self.assertEqual(observed, expected)
        observed = HGVSProtChange.fromStr("(Pro458_Gly460del)")
        self.assertEqual(observed, expected)
        expected = HGVSProtChange("Gly", 2, "Met", 46, None, "del", [], False)
        observed = HGVSProtChange.fromStr("G2_M46del")
        self.assertEqual(observed, expected)
        observed = HGVSProtChange.fromStr("Gly2_Met46del")
        self.assertEqual(observed, expected)
        expected = HGVSProtChange("Trp", 26, None, None, None, None, ["Ter"], False)
        observed = HGVSProtChange.fromStr("W26*")
        self.assertEqual(observed, expected)
        observed = HGVSProtChange.fromStr("Trp26*")
        self.assertEqual(observed, expected)
        expected = HGVSProtChange("Pro", 551, "Glu", 554, None, "del", [], False)
        observed = HGVSProtChange.fromStr("P551_E554delPMYE")
        self.assertEqual(observed, expected)
        observed = HGVSProtChange.fromStr("Pro551_Glu554delProMetTyrGlu")
        self.assertEqual(observed, expected)
        expected = HGVSProtChange("Trp", 557, "Lys", 558, None, "del", [], False)
        observed = HGVSProtChange.fromStr("W557_K558DELWK")
        self.assertEqual(observed, expected)
        observed = HGVSProtChange.fromStr("Trp557_Lys558DelTrpLys")
        self.assertEqual(observed, expected)
        expected = HGVSProtChange("Asp", 842, "His", 845, None, "del", [], False)
        observed = HGVSProtChange.fromStr("D842_H845DelDIMH")
        self.assertEqual(observed, expected)
        observed = HGVSProtChange.fromStr("Asp842_His845DelAspIleMetHis")
        self.assertEqual(observed, expected)
        expected = HGVSProtChange("Asp", 842, None, None, None, "del", [], False)
        observed = HGVSProtChange.fromStr("D842DELD")
        self.assertEqual(observed, expected)
        observed = HGVSProtChange.fromStr("Asp842DelAsp")
        self.assertEqual(observed, expected)
        # observed = HGVSProtChange.fromStr("Val7=/del")
        # Duplication
        expected = HGVSProtChange("Ala", 767, "Val", 769, None, "dup", [], False)
        observed = HGVSProtChange.fromStr("A767_V769dupASV")
        self.assertEqual(observed, expected)
        expected = HGVSProtChange("Ala", 3, None, None, None, "dup", [], False)
        observed = HGVSProtChange.fromStr("Ala3dup")
        self.assertEqual(observed, expected)
        observed = HGVSProtChange.fromStr("A3dup")
        self.assertEqual(observed, expected)
        expected = HGVSProtChange("Ala", 3, "Ser", 5, None, "dup", [], False)
        observed = HGVSProtChange.fromStr("A3_S5dup")
        self.assertEqual(observed, expected)
        observed = HGVSProtChange.fromStr("Ala3_Ser5dup")
        self.assertEqual(observed, expected)
        expected = HGVSProtChange("Ser", 6, None, None, None, "dup", [], False)
        observed = HGVSProtChange.fromStr("S6dup")
        self.assertEqual(observed, expected)
        observed = HGVSProtChange.fromStr("Ser6dup")
        self.assertEqual(observed, expected)
        # Insertion
        expected = HGVSProtChange("His", 4, "Gln", 5, None, "ins", ["Ala"], False)
        observed = HGVSProtChange.fromStr("H4_Q5insA")
        self.assertEqual(observed, expected)
        observed = HGVSProtChange.fromStr("His4_Gln5insAla")
        self.assertEqual(observed, expected)
        expected = HGVSProtChange("Lys", 2, "Gly", 3, None, "ins", ["Gln", "Ser", "Lys"], False)
        observed = HGVSProtChange.fromStr("K2_G3insQSK")
        self.assertEqual(observed, expected)
        observed = HGVSProtChange.fromStr("Lys2_Gly3insGlnSerLys")
        self.assertEqual(observed, expected)
        expected = HGVSProtChange("Met", 3, "His", 4, None, "ins", ["Gly", "Ter"], True)
        observed = HGVSProtChange.fromStr("(M3_H4insG*)")
        self.assertEqual(observed, expected)
        observed = HGVSProtChange.fromStr("(Met3_His4insGlyTer)")
        self.assertEqual(observed, expected)
        expected = HGVSProtChange("Arg", 78, "Gly", 79, None, "ins", ["23"], False)
        observed = HGVSProtChange.fromStr("R78_G79ins23")
        self.assertEqual(observed, expected)
        observed = HGVSProtChange.fromStr("Arg78_Gly79ins23")
        self.assertEqual(observed, expected)
        expected = HGVSProtChange("Gln", 746, "Lys", 747, None, "ins", ["Ter", "63"], False)
        observed = HGVSProtChange.fromStr("Q746_K747ins*63")
        self.assertEqual(observed, expected)
        observed = HGVSProtChange.fromStr("Gln746_Lys747ins*63")
        self.assertEqual(observed, expected)
        expected = HGVSProtChange("Ser", 332, "Ser", 333, None, "ins", ["(1)"], True)
        observed = HGVSProtChange.fromStr("(S332_S333ins(1))")
        self.assertEqual(observed, expected)
        observed = HGVSProtChange.fromStr("(Ser332_Ser333ins(1))")
        self.assertEqual(observed, expected)
        expected = HGVSProtChange("Ser", 332, "Ser", 333, None, "ins", ["Xaa"], True)
        observed = HGVSProtChange.fromStr("(S332_S333insX)")
        self.assertEqual(observed, expected)
        observed = HGVSProtChange.fromStr("(Ser332_Ser333insX)")
        self.assertEqual(observed, expected)
        expected = HGVSProtChange("Val", 582, "Asn", 583, None, "ins", ["(5)"], True)
        observed = HGVSProtChange.fromStr("(V582_N583ins(5))")
        self.assertEqual(observed, expected)
        observed = HGVSProtChange.fromStr("(Val582_Asn583ins(5))")
        self.assertEqual(observed, expected)
        expected = HGVSProtChange("Val", 582, "Asn", 583, None, "ins", ["Xaa", "Xaa", "Xaa", "Xaa", "Xaa"], True)
        observed = HGVSProtChange.fromStr("(Val582_Asn583insXXXXX)")
        self.assertEqual(observed, expected)
        observed = HGVSProtChange.fromStr("(Val582_Asn583insXaaXaaXaaXaaXaa)")
        self.assertEqual(observed, expected)
        # Frameshift
        expected = HGVSProtChange("Val", 2288, None, None, None, "fs", ["Ter", "1"], False)
        observed = HGVSProtChange.fromStr("V2288FS*1")
        self.assertEqual(observed, expected)
        expected = HGVSProtChange("Arg", 97, None, None, "Pro", "fs", ["Ter", "23"], False)
        observed = HGVSProtChange.fromStr("R97Pfs*23")
        self.assertEqual(observed, expected)
        observed = HGVSProtChange.fromStr("Arg97ProfsTer23")
        self.assertEqual(observed, expected)
        expected = HGVSProtChange("Arg", 97, None, None, None, "fs", [], False)
        observed = HGVSProtChange.fromStr("R97FS")
        self.assertEqual(observed, expected)
        observed = HGVSProtChange.fromStr("Arg97fs")
        self.assertEqual(observed, expected)
        expected = HGVSProtChange("Tyr", 4, None, None, None, None, ["Ter"], True)
        observed = HGVSProtChange.fromStr("(Y4*)")
        self.assertEqual(observed, expected)
        observed = HGVSProtChange.fromStr("(Tyr4*)")
        self.assertEqual(observed, expected)
        expected = HGVSProtChange("Ile", 327, None, None, "Arg", "fs", ["Ter", "?"], False)
        observed = HGVSProtChange.fromStr("I327Rfs*?")
        self.assertEqual(observed, expected)
        observed = HGVSProtChange.fromStr("Ile327Argfs*?")
        self.assertEqual(observed, expected)
        expected = HGVSProtChange("Gln", 151, None, None, "Thr", "fs", ["Ter", "9"], False)
        observed = HGVSProtChange.fromStr("Q151TFS*9")
        self.assertEqual(observed, expected)
        observed = HGVSProtChange.fromStr("Gln151Thrfs*9")
        self.assertEqual(observed, expected)
        # Extension
        expected = HGVSProtChange("Met", 1, None, None, None, "ext", ["-5"], False)
        observed = HGVSProtChange.fromStr("M1ext-5")
        self.assertEqual(observed, expected)
        observed = HGVSProtChange.fromStr("Met1ext-5")
        self.assertEqual(observed, expected)
        expected = HGVSProtChange("Ter", 110, None, None, "Gln", "ext", ["Ter", "17"], False)
        observed = HGVSProtChange.fromStr("*110Qext*17")
        self.assertEqual(observed, expected)
        observed = HGVSProtChange.fromStr("Ter110GlnextTer17")
        self.assertEqual(observed, expected)
        observed = HGVSProtChange.fromStr("*110Glnext*17")
        self.assertEqual(observed, expected)
        expected = HGVSProtChange("Ter", 315, None, None, "Tyr", "ext", ["Asn", "Lys", "Gly", "Thr", "Ter"], True)
        observed = HGVSProtChange.fromStr("(*315YextNKGT*)")
        self.assertEqual(observed, expected)
        observed = HGVSProtChange.fromStr("(Ter315TyrextAsnLysGlyThrTer)")
        self.assertEqual(observed, expected)
        expected = HGVSProtChange("Ter", 315, None, None, "Tyr", "ext", ["Asn", "Lys", "Gly", "Thr", "Ter"], False)
        observed = HGVSProtChange.fromStr("*315TyrextAsnLysGlyThr*")
        self.assertEqual(observed, expected)
        expected = HGVSProtChange("Ter", 327, None, None, "Arg", "ext", ["Ter", "?"], False)
        observed = HGVSProtChange.fromStr("*327Rext*?")
        self.assertEqual(observed, expected)
        observed = HGVSProtChange.fromStr("Ter327Argext*?")
        self.assertEqual(observed, expected)
        # Reapeated sequences
        expected = HGVSProtChange("Ala", 2, None, None, None, None, ["[", "10", "]"], False)
        observed = HGVSProtChange.fromStr("A2[10]")
        self.assertEqual(observed, expected)
        expected = HGVSProtChange("Ala", 2, None, None, None, None, ["[", "10", "]"], False)
        observed = HGVSProtChange.fromStr("Ala2[10]")
        self.assertEqual(observed, expected)
        expected = HGVSProtChange("Ala", 2, None, None, None, "dup", [], False)
        observed = HGVSProtChange.fromStr("A2[1]")
        self.assertEqual(observed, expected)
        expected = HGVSProtChange("Ala", 2, None, None, None, "dup", [], False)
        observed = HGVSProtChange.fromStr("Ala2[1]")
        self.assertEqual(observed, expected)
        expected = HGVSProtChange("Ala", 2, "Pro", 5, None, None, ["[", "10", "]"], False)
        observed = HGVSProtChange.fromStr("A2_P5[10]")
        self.assertEqual(observed, expected)
        expected = HGVSProtChange("Ala", 2, "Pro", 5, None, None, ["[", "10", "]"], False)
        observed = HGVSProtChange.fromStr("Ala2_Pro5[10]")
        self.assertEqual(observed, expected)
        expected = HGVSProtChange("Ala", 2, "Pro", 5, None, "dup", [], False)
        observed = HGVSProtChange.fromStr("A2_P5[1]")
        self.assertEqual(observed, expected)
        expected = HGVSProtChange("Ala", 2, "Pro", 5, None, "dup", [], False)
        observed = HGVSProtChange.fromStr("Ala2_Pro5[1]")
        self.assertEqual(observed, expected)


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    unittest.main()
