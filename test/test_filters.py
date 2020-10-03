#!/usr/bin/env python3

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2017 IUCT-O'
__license__ = 'GNU General Public License'
__version__ = '1.1.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import os
import sys
import warnings
import unittest

TEST_DIR = os.path.dirname(os.path.abspath(__file__))
PACKAGE_DIR = os.path.dirname(TEST_DIR)
sys.path.append(PACKAGE_DIR)

from anacore.filters import Filter, FiltersCombiner, filtersFromDict


########################################################################
#
# FUNCTIONS
#
########################################################################
class TestFilter(unittest.TestCase):
    def setUp(self):
        # Data
        self.data = [
            {"age": 8, "treatment": "placebo", "group":["A", "B"]},
            {"age": 9, "treatment": "20ng ip", "group":["C"]},
            {"age": 11, "treatment": "placebo", "group":["A", "D"]},
            {"age": 25, "treatment": "20ng ip", "group":["B"]},
            {"age": 30, "treatment": "20ng pc", "group":["C", "B"]},
            {"age": 33, "treatment": "20ng pc", "group":["A", "D"]},
            {"age": 30, "treatment": "placebo", "group":["A", "B"]},
        ]


    def testFromDict(self):
        # Only required attributes are provided
        filter_obj = Filter.fromDict({
            "class": "Filter",
            "operator": "<",
            "values": 25
        })
        expected = Filter("<", 25)
        self.assertEqual(expected, filter_obj)
        # All attributes are provided
        filter_obj = Filter.fromDict({
            "action": "exclude",
            "aggregator": "nb:1",
            "class": "Filter",
            "description": "Select patients with age lower than 25 years",
            "getter": "age",
            "name": "young",
            "operator": "<",
            "values": 25
        })
        expected = Filter("<", 25, "age", "nb:1", "young", "Select patients with age lower than 25 years", "exclude")
        self.assertEqual(expected, filter_obj)


    def testToDict(self):
        # Only required attributes are provided
        filter_obj = Filter("<", 25)
        expected = {
            "action": "select",
            "aggregator": "ratio:1",
            "class": "Filter",
            "description": None,
            "getter": None,
            "name": None,
            "operator": "<",
            "values": 25
        }
        self.assertEqual(expected, filter_obj.toDict())
        # All attributes are provided
        filter_obj = Filter("<", 25, "age", "nb:1", "young", "Select patients with age lower than 25 years", "exclude")
        expected = {
            "action": "exclude",
            "aggregator": "nb:1",
            "class": "Filter",
            "description": "Select patients with age lower than 25 years",
            "getter": "age",
            "name": "young",
            "operator": "<",
            "values": 25
        }
        self.assertEqual(expected, filter_obj.toDict())


    def testSetFct(self):
        # Test operator contains with list threshold
        with self.assertRaises(AttributeError) as context:
            Filter("contains", ["ip", "pc"], "treatment")  # Select patients with ip or pc in treatment text

        # Test unknown operator
        with self.assertRaises(AttributeError) as context:
            Filter("unknown", ["ip", "pc"], "treatment")


    def testEval(self):
        # Test action: select AND numeric threshold
        filter_obj = Filter("<", 13, "age")  # Select patients with age not < 13
        expected = [True, True, True, False, False, False, False]
        observed = [filter_obj.eval(curr) for curr in self.data]
        self.assertEqual(expected, observed)

        # Test action: exclude AND numeric threshold
        filter_obj = Filter("<", 25, "age", action="exclude")  # Select patients with age not < 25 (>=25)
        expected = [False, False, False, True, True, True, True]
        observed = [filter_obj.eval(curr) for curr in self.data]
        self.assertEqual(expected, observed)

        # Test action: select AND string threshold
        filter_obj = Filter("!=", "placebo", "treatment")  # Select patients with treatment different of placebo
        expected = [False, True, False, True, True, True, False]
        observed = [filter_obj.eval(curr) for curr in self.data]
        self.assertEqual(expected, observed)

        # Test action: exclude AND string threshold
        filter_obj = Filter("=", "20ng pc", "treatment", action="exclude")  # Select patients with treatment different of 20ng pc
        expected = [True, True, True, True, False, False, True]
        observed = [filter_obj.eval(curr) for curr in self.data]
        self.assertEqual(expected, observed)

        # Test operator: contains AND string threshold
        filter_obj = Filter("contains", "20ng", "treatment")  # Select patients with 20ng in treatment text
        expected = [False, True, False, True, True, True, False]
        observed = [filter_obj.eval(curr) for curr in self.data]
        self.assertEqual(expected, observed)

        # Test operator: in AND list threshold
        filter_obj = Filter("in", ["20ng pc", "placebo"], "treatment")  # Select patients with "20ng pc" or "placebo" as treatment
        expected = [True, False, True, False, True, True, True]
        observed = [filter_obj.eval(curr) for curr in self.data]
        self.assertEqual(expected, observed)

        # Test aggregator: nb:1
        filter_obj = Filter("=", "A", "group", aggregator="nb:1")  # Select patients with at least the group A
        expected = [True, False, True, False, False, True, True]
        observed = [filter_obj.eval(curr) for curr in self.data]
        self.assertEqual(expected, observed)

        # Test aggregator: ratio:1
        filter_obj = Filter("=", "C", "group", aggregator="ratio:1")  # Select patients with 100% of grooups are equals to C
        expected = [False, True, False, False, False, False, False]
        observed = [filter_obj.eval(curr) for curr in self.data]
        self.assertEqual(expected, observed)

        # Test aggregator: ratio:0.5
        filter_obj = Filter("=", "C", "group", aggregator="ratio:0.5")  # Select patients with 50% of grooups are equals to C
        expected = [False, True, False, False, True, False, False]
        observed = [filter_obj.eval(curr) for curr in self.data]
        self.assertEqual(expected, observed)

        # Test list threshold AND operator: in AND aggregator: nb:1
        filter_obj = Filter("in", ["A", "C"], "group", aggregator="nb:1")  # Select patients with at least A or C in groups
        expected = [True, True, True, False, True, True, True]
        observed = [filter_obj.eval(curr) for curr in self.data]
        self.assertEqual(expected, observed)

        # Test list threshold AND operator: not in AND aggregator: nb:1
        filter_obj = Filter("not in", ["B", "C"], "group", aggregator="nb:1")  # Select patients with at least one group not in B and C
        expected = [True, False, True, False, False, True, True]
        observed = [filter_obj.eval(curr) for curr in self.data]
        self.assertEqual(expected, observed)

        # Test list threshold AND aggregator: ratio:1
        filter_obj = Filter("in", ["B", "C"], "group", aggregator="ratio:1")  # Select patients where groups are B and/or C
        expected = [False, True, False, True, True, False, False]
        observed = [filter_obj.eval(curr) for curr in self.data]
        self.assertEqual(expected, observed)

        filter_obj = Filter("not in", ["B", "C"], "group", aggregator="ratio:1")  # Select patients where neither group is B or C
        expected = [False, False, True, False, False, True, False]
        observed = [filter_obj.eval(curr) for curr in self.data]
        self.assertEqual(expected, observed)

        # Test None behavior in simple filter
        data_age = [{"age": None}, {"age": 10}, {"age": 1}]
        filter_obj = Filter("<", 5, "age")
        expected = [False, False, True]
        observed = [filter_obj.eval(curr) for curr in data_age]
        self.assertEqual(expected, observed)
        filter_obj = Filter(">", 5, "age")
        expected = [False, True, False]
        observed = [filter_obj.eval(curr) for curr in data_age]
        self.assertEqual(expected, observed)
        filter_obj = Filter("==", None, "age")
        expected = [True, False, False]
        observed = [filter_obj.eval(curr) for curr in data_age]
        self.assertEqual(expected, observed)
        filter_obj = Filter("<>", None, "age")
        expected = [False, True, True]
        observed = [filter_obj.eval(curr) for curr in data_age]
        self.assertEqual(expected, observed)
        filter_obj = Filter("in", [10, 1], "age")
        expected = [False, True, True]
        observed = [filter_obj.eval(curr) for curr in data_age]
        self.assertEqual(expected, observed)
        filter_obj = Filter("not in", [10, 1], "age")
        expected = [True, False, False]
        observed = [filter_obj.eval(curr) for curr in data_age]
        self.assertEqual(expected, observed)

        # Test None behavior in aggregated
        data_class = [{"class": None}, {"class": [None, "A"]}, {"class": ["B"]}]
        filter_obj = Filter("==", "A", "class", "nb:1")
        expected = [False, True, False]
        observed = [filter_obj.eval(curr) for curr in data_class]
        self.assertEqual(expected, observed)
        filter_obj = Filter("<>", "A", "class", "ratio:1")
        expected = [True, False, True]
        observed = [filter_obj.eval(curr) for curr in data_class]
        self.assertEqual(expected, observed)
        # filter_obj = Filter("==", None, "class", "nb:1")
        # expected = [False, True, False]
        # observed = [filter_obj.eval(curr) for curr in data_class]
        # self.assertEqual(expected, observed)
        filter_obj = Filter("not in", ["A", "B"], "class", "nb:1")
        expected = [True, True, False]
        observed = [filter_obj.eval(curr) for curr in data_class]
        self.assertEqual(expected, observed)
        filter_obj = Filter("not in", ["A", "B"], "class", "ratio:1")
        expected = [True, False, False]
        observed = [filter_obj.eval(curr) for curr in data_class]
        self.assertEqual(expected, observed)
        filter_obj = Filter("in", [None, "B"], "class", "ratio:1")
        expected = [True, False, True]
        observed = [filter_obj.eval(curr) for curr in data_class]
        self.assertEqual(expected, observed)
        filter_obj = Filter("in", [None, "B"], "class", "nb:1")
        expected = [True, True, True]
        observed = [filter_obj.eval(curr) for curr in data_class]
        self.assertEqual(expected, observed)


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    warnings.warn("TODO: TestFiltersCombiner")
    warnings.warn("TODO: testFiltersFromDict")
    warnings.warn("TODO: TestFilter.testGetRecordValue")
    unittest.main()
