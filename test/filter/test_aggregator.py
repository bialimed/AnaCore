#!/usr/bin/env python3

__author__ = 'Rémi THEVENOUX'
__copyright__ = 'Copyright (C) 2023 CHU Toulouse'
__license__ = 'GNU General Public License'
__version__ = '1.0.0'
__status__ = 'prod'

import os
import sys
import unittest
from dict_util import flattenDict

TEST_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PACKAGE_DIR = os.path.dirname(TEST_DIR)
sys.path.append(PACKAGE_DIR)


#
# Temporay class (must be implemented)
#
class AggregatorFunction:
    @staticmethod
    def build(input):
        return lambda positive, total: True


########################################################################
#
# Test cases
#
########################################################################
class TestAggregator(unittest.TestCase):

    def testInvalidAggregator(self):
        suite = {
            "invalid type": "foo:0.5",
            "nb": {
                "no colon": "nb",
                "double colon 1": "nb::2",
                "double colon 2": "nb:2:",
                "empty": "nb:",
                "not a number 1": "nb:foo",
                "not a number 2": "nb:4foo",
                "not an integer": "nb:4.2",
                "zero": "nb:0",
                "negative": "nb:-3",
            },
            "ratio": {
                "no colon": "ratio",
                "double colon 1": "ratio::0.5",
                "double colon 2": "ratio:0.5:",
                "empty": "ratio:",
                "not a number-1": "ratio:foo",
                "not a number-2": "ratio:0.foo",
                "negative": "ratio:-0.5",
                "zero": "ratio:0.0",
                ">1": "ratio:3.14",
            }
        }

        case_list = flattenDict(suite)
        for caseName, caseValue in case_list.items():
            with self.subTest(caseName):
                try:
                    AggregatorFunction.build(caseValue)
                except Exception:
                    return  # fail as expected
                raise AssertionError("Aggregator build does not fail")

    def testAggregator(self):
        case_dict = {
            "nb": {
                "0": ["nb:2", 0, 5, False],
                "<": ["nb:2", 1, 5, False],
                "=": ["nb:2", 2, 5, True],
                ">": ["nb:2", 4, 5, True],
                "0/0": ["nb:2", 0, 0, False],
            },
            "ratio": {
                "0": ["ratio:0.5", 0, 6, False],
                "<": ["ratio:0.5", 2, 6, False],
                "=": ["ratio:0.5", 3, 6, True],
                ">": ["ratio:0.5", 5, 6, True],
                "0/0": {
                    "50%": ["ratio:0.5", 0, 0, True],
                    "all": ["ratio:1", 0, 0, True],
                }
            }
        }

        case_list = flattenDict(case_dict)
        for caseName, caseValues in case_list.items():
            with self.subTest(caseName):
                aggregator = AggregatorFunction.build(caseValues[0])
                result = aggregator(caseValues[1], caseValues[2])
                self.assertEqual(result, caseValues[3])


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    unittest.main()
