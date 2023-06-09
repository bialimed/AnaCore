#!/usr/bin/env python3

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2017 CHU Toulouse'
__license__ = 'GNU General Public License'
__version__ = '1.2.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import os
import sys
import unittest

TEST_DIR = os.path.dirname(os.path.abspath(__file__))
PACKAGE_DIR = os.path.dirname(TEST_DIR)
sys.path.append(PACKAGE_DIR)

from anacore.filters import EmptyIterFilter, Filter, FiltersCombiner, filtersFromDict


########################################################################
#
# FUNCTIONS
#
########################################################################
class TestEmptyIterFilter(unittest.TestCase):
    def setUp(self):
        # Data
        self.data = [
            {"age": 8, "group": ["A", "B"], "treatment_duration": {"tab": 8}},
            {"age": 18, "group": ["C"], "treatment_duration": {"tab": 8, "afa": 6}},
            {"age": 28, "group": [], "treatment_duration": {}},
            {"age": 38, "group": None, "treatment_duration": None},
            {"age": 48}
        ]

    def testEvalDict(self):
        # Select
        filter_obj = EmptyIterFilter("treatment_duration")
        expected = [False, False, True, True, True]
        observed = [filter_obj.eval(curr) for curr in self.data]
        self.assertEqual(observed, expected)
        # Exclude
        filter_obj = EmptyIterFilter("treatment_duration", action="exclude")
        expected = [True, True, False, False, False]
        observed = [filter_obj.eval(curr) for curr in self.data]
        self.assertEqual(observed, expected)

    def testEvalList(self):
        # Select
        filter_obj = EmptyIterFilter("group")
        expected = [False, False, True, True, True]
        observed = [filter_obj.eval(curr) for curr in self.data]
        self.assertEqual(observed, expected)
        # Exclude
        filter_obj = EmptyIterFilter("group", action="exclude")
        expected = [True, True, False, False, False]
        observed = [filter_obj.eval(curr) for curr in self.data]
        self.assertEqual(observed, expected)

    def testFromDict(self):
        # Only required attributes are provided
        filter_obj = EmptyIterFilter.fromDict({})
        expected = EmptyIterFilter()
        self.assertEqual(expected, filter_obj)
        # All attributes are provided
        filter_obj = EmptyIterFilter.fromDict({
            "action": "exclude",
            "class": "EmptyIterFilter",
            "description": "Exclude patients in no groups",
            "getter": "group",
            "name": "empty_gp_exclusion"
        })
        expected = EmptyIterFilter("group", "empty_gp_exclusion", "Exclude patients in no groups", "exclude")
        self.assertEqual(expected, filter_obj)

    def testToDict(self):
        # Only required attributes are provided
        filter_obj = EmptyIterFilter()
        expected = {
            "action": "select",
            "class": "EmptyIterFilter",
            "description": None,
            "getter": None,
            "name": None
        }
        self.assertEqual(expected, filter_obj.toDict())
        # All attributes are provided
        filter_obj = EmptyIterFilter("group", "empty_gp_exclusion", "Exclude patients in no groups", "exclude")
        expected = {
            "action": "exclude",
            "class": "EmptyIterFilter",
            "description": "Exclude patients in no groups",
            "getter": "group",
            "name": "empty_gp_exclusion"
        }
        self.assertEqual(expected, filter_obj.toDict())


class TestFilter(unittest.TestCase):
    def setUp(self):
        # Data
        self.data = [
            {"age": 8, "treatment": "placebo", "group": ["A", "B"]},
            {"age": 9, "treatment": "20ng ip", "group": ["C"]},
            {"age": 11, "treatment": "placebo", "group": ["A", "D"]},
            {"age": 25, "treatment": "20ng ip", "group": ["B"]},
            {"age": 30, "treatment": "20ng pc", "group": ["C", "B"]},
            {"age": 33, "treatment": "20ng pc", "group": ["A", "D"]},
            {"age": 30, "treatment": "placebo", "group": ["A", "B"]},
        ]

    def testAggregator(self):
        # nb:1
        filter = Filter("in", ["A", "B"], "group", aggregator="nb:1")
        expected = [True, False, True, True, True, True, True]
        observed = [filter.eval(elt) for elt in self.data]
        self.assertEqual(observed, expected)
        # nb:2
        filter = Filter("in", ["A", "B"], "group", aggregator="nb:2")
        expected = [True, False, False, False, False, False, True]
        observed = [filter.eval(elt) for elt in self.data]
        self.assertEqual(observed, expected)
        # nb:0
        filter = Filter("in", ["A", "B"], "group", aggregator="nb:0")
        expected = [True, True, True, True, True, True, True]
        observed = [filter.eval(elt) for elt in self.data]
        self.assertEqual(observed, expected)
        # nb:1 in None
        filter = Filter("==", "pathogen", "status", aggregator="nb:1")
        expected = [False, False, False, False, False, False, False]
        observed = [filter.eval(elt) for elt in self.data]
        self.assertEqual(observed, expected)
        # ratio:1
        filter = Filter("in", ["A", "B"], "group", aggregator="ratio:1")
        expected = [True, False, False, True, False, False, True]
        observed = [filter.eval(elt) for elt in self.data]
        self.assertEqual(observed, expected)
        # ratio:0.5
        filter = Filter("in", ["A", "B"], "group", aggregator="ratio:0.5")
        expected = [True, False, True, True, True, True, True]
        observed = [filter.eval(elt) for elt in self.data]
        self.assertEqual(observed, expected)
        # ratio:None
        filter = Filter("==", "pathogen", "group", aggregator="ratio:1")
        expected = [False, False, False, False, False, False, False]
        observed = [filter.eval(elt) for elt in self.data]
        self.assertEqual(observed, expected)

    def testAggregatorFormat(self):
        filter_obj = Filter("==", "A")
        # Aggregator nb
        filter_obj.aggregator = "nb:1"
        with self.assertRaises(ValueError):
            filter_obj.aggregator = "nb:-1"
        with self.assertRaises(ValueError):
            filter_obj.aggregator = "nb:0.1"
        # Aggregator ratio
        filter_obj.aggregator = "ratio:1"
        filter_obj.aggregator = "ratio:1.0"
        with self.assertRaises(ValueError):
            filter_obj.aggregator = "ratio:-0.5"
        with self.assertRaises(ValueError):
            filter_obj.aggregator = "ratio:2"

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
            "aggregator": "nb:1",  # No sense in this context it is only for test
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
            "aggregator": None,
            "class": "Filter",
            "description": None,
            "getter": None,
            "name": None,
            "operator": "<",
            "values": 25
        }
        self.assertEqual(expected, filter_obj.toDict())
        # All attributes are provided
        filter_obj = Filter("<", 25, "age", "ratio:1", "young", "Select patients with age lower than 25 years", "exclude")
        expected = {
            "action": "exclude",
            "aggregator": "ratio:1",
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
        self.assertEqual(observed, expected)

        # Test action: exclude AND numeric threshold
        filter_obj = Filter("<", 25, "age", action="exclude")  # Select patients with age not < 25 (>=25)
        expected = [False, False, False, True, True, True, True]
        observed = [filter_obj.eval(curr) for curr in self.data]
        self.assertEqual(observed, expected)

        # Test action: select AND string threshold
        filter_obj = Filter("!=", "placebo", "treatment")  # Select patients with treatment different of placebo
        expected = [False, True, False, True, True, True, False]
        observed = [filter_obj.eval(curr) for curr in self.data]
        self.assertEqual(observed, expected)

        # Test action: exclude AND string threshold
        filter_obj = Filter("=", "20ng pc", "treatment", action="exclude")  # Select patients with treatment different of 20ng pc
        expected = [True, True, True, True, False, False, True]
        observed = [filter_obj.eval(curr) for curr in self.data]
        self.assertEqual(observed, expected)

        # Test operator: contains AND string threshold
        filter_obj = Filter("contains", "20ng", "treatment")  # Select patients with 20ng in treatment text
        expected = [False, True, False, True, True, True, False]
        observed = [filter_obj.eval(curr) for curr in self.data]
        self.assertEqual(observed, expected)

        # Test operator: contains AND string threshold for int value
        filter_obj = Filter("contains", "1", "age")  # Select patients with 1 in age
        expected = [False, False, True, False, False, False, False]
        observed = [filter_obj.eval(curr) for curr in self.data]
        self.assertEqual(observed, expected)

        # Test operator: substring of AND string threshold
        filter_obj = Filter("substring of", "large placebo", "treatment")  # Select patients where treatment text is substring of "large placebo"
        expected = [True, False, True, False, False, False, True]
        observed = [filter_obj.eval(curr) for curr in self.data]
        self.assertEqual(observed, expected)

        # Test operator: substring of AND string threshold
        filter_obj = Filter("substring of", "3087", "age")  # Select patients where age is substring of "3087"
        expected = [True, False, False, False, True, False, True]
        observed = [filter_obj.eval(curr) for curr in self.data]
        self.assertEqual(observed, expected)

        # Test operator: in AND list threshold
        filter_obj = Filter("in", ["20ng pc", "placebo"], "treatment")  # Select patients with "20ng pc" or "placebo" as treatment
        expected = [True, False, True, False, True, True, True]
        observed = [filter_obj.eval(curr) for curr in self.data]
        self.assertEqual(observed, expected)

        # Test aggregator: nb:1
        filter_obj = Filter("=", "A", "group", aggregator="nb:1")  # Select patients with at least the group A
        expected = [True, False, True, False, False, True, True]
        observed = [filter_obj.eval(curr) for curr in self.data]
        self.assertEqual(observed, expected)

        # Test aggregator: ratio:1
        filter_obj = Filter("=", "C", "group", aggregator="ratio:1")  # Select patients with 100% of grooups are equals to C
        expected = [False, True, False, False, False, False, False]
        observed = [filter_obj.eval(curr) for curr in self.data]
        self.assertEqual(observed, expected)

        # Test aggregator: ratio:0.5
        filter_obj = Filter("=", "C", "group", aggregator="ratio:0.5")  # Select patients with 50% of grooups are equals to C
        expected = [False, True, False, False, True, False, False]
        observed = [filter_obj.eval(curr) for curr in self.data]
        self.assertEqual(observed, expected)

        # Test list threshold AND operator: in AND aggregator: nb:1
        filter_obj = Filter("in", ["A", "C"], "group", aggregator="nb:1")  # Select patients with at least A or C in groups
        expected = [True, True, True, False, True, True, True]
        observed = [filter_obj.eval(curr) for curr in self.data]
        self.assertEqual(observed, expected)

        # Test list threshold AND operator: not in AND aggregator: nb:1
        filter_obj = Filter("not in", ["B", "C"], "group", aggregator="nb:1")  # Select patients with at least one group not in B and C
        expected = [True, False, True, False, False, True, True]
        observed = [filter_obj.eval(curr) for curr in self.data]
        self.assertEqual(observed, expected)

        # Test list threshold AND aggregator: ratio:1
        filter_obj = Filter("in", ["B", "C"], "group", aggregator="ratio:1")  # Select patients where groups are B and/or C
        expected = [False, True, False, True, True, False, False]
        observed = [filter_obj.eval(curr) for curr in self.data]
        self.assertEqual(observed, expected)

        filter_obj = Filter("not in", ["B", "C"], "group", aggregator="ratio:1")  # Select patients where neither group is B or C
        expected = [False, False, True, False, False, True, False]
        observed = [filter_obj.eval(curr) for curr in self.data]
        self.assertEqual(observed, expected)

        # Test None behavior in simple filter
        data_age = [{"age": None}, {"age": 10}, {"age": 1}]
        filter_obj = Filter("<", 5, "age")
        expected = [False, False, True]
        observed = [filter_obj.eval(curr) for curr in data_age]
        self.assertEqual(observed, expected)
        filter_obj = Filter(">", 5, "age")
        expected = [False, True, False]
        observed = [filter_obj.eval(curr) for curr in data_age]
        self.assertEqual(observed, expected)
        filter_obj = Filter("==", None, "age")
        expected = [True, False, False]
        observed = [filter_obj.eval(curr) for curr in data_age]
        self.assertEqual(observed, expected)
        filter_obj = Filter("<>", None, "age")
        expected = [False, True, True]
        observed = [filter_obj.eval(curr) for curr in data_age]
        self.assertEqual(observed, expected)
        filter_obj = Filter("in", [10, 1], "age")
        expected = [False, True, True]
        observed = [filter_obj.eval(curr) for curr in data_age]
        self.assertEqual(observed, expected)
        filter_obj = Filter("not in", [10, 1], "age")
        expected = [True, False, False]
        observed = [filter_obj.eval(curr) for curr in data_age]
        self.assertEqual(observed, expected)

        # Test None behavior in aggregated
        data_class = [{"class": None}, {"class": []}, {"class": [None, "A"]}, {"class": ["B"]}]
        filter_obj = Filter("==", "A", "class", "nb:1")
        expected = [False, False, True, False]
        observed = [filter_obj.eval(curr) for curr in data_class]
        self.assertEqual(observed, expected)
        filter_obj = Filter("<>", "A", "class", "nb:1")
        expected = [False, False, True, True]
        observed = [filter_obj.eval(curr) for curr in data_class]
        self.assertEqual(observed, expected)
        filter_obj = Filter("==", "A", "class", "ratio:1")
        expected = [True, True, False, False]
        observed = [filter_obj.eval(curr) for curr in data_class]
        self.assertEqual(observed, expected)
        filter_obj = Filter("<>", "A", "class", "ratio:1")
        expected = [True, True, False, True]
        observed = [filter_obj.eval(curr) for curr in data_class]
        self.assertEqual(observed, expected)
        filter_obj = Filter("==", None, "class", "nb:1")
        expected = [False, False, True, False]
        observed = [filter_obj.eval(curr) for curr in data_class]
        self.assertEqual(observed, expected)
        filter_obj = Filter("not in", ["A", "B"], "class", "nb:1")
        expected = [False, False, True, False]
        observed = [filter_obj.eval(curr) for curr in data_class]
        self.assertEqual(observed, expected)
        filter_obj = Filter("not in", ["A", "B"], "class", "ratio:1")
        expected = [True, True, False, False]
        observed = [filter_obj.eval(curr) for curr in data_class]
        self.assertEqual(observed, expected)
        filter_obj = Filter("in", [None, "B"], "class", "ratio:1")
        expected = [True, True, False, True]
        observed = [filter_obj.eval(curr) for curr in data_class]
        self.assertEqual(observed, expected)
        filter_obj = Filter("in", [None, "B"], "class", "nb:1")
        expected = [False, False, True, True]
        observed = [filter_obj.eval(curr) for curr in data_class]
        self.assertEqual(observed, expected)


class TestFiltersCombiner(unittest.TestCase):
    def setUp(self):
        self.data = [
            {"age": 8, "group": ["A", "B"]},
            {"age": 18, "group": ["C"]},
            {"age": 28, "group": []},
            {"age": 38, "group": None},
            {"age": 48}
        ]

    def testEqual(self):
        filter_def = {
            "class": "FiltersCombiner",
            "description": "Age < 10 and in group.",
            "name": "young_and_in_group",
            "operator": "and",
            "filters": [
                {
                    "class": "Filter",
                    "getter": "age",
                    "operator": "<",
                    "values": 10
                },
                {
                    "action": "exclude",
                    "class": "EmptyIterFilter",
                    "getter": "group"
                }
            ]
        }
        filters = filtersFromDict(filter_def)
        # Equal
        filters_2 = filtersFromDict(filter_def)
        self.assertEqual(filters, filters_2)
        # Equal after remove empty FiltersCombiner
        filters_2 = filtersFromDict(filter_def)
        filters_2.filters.append(FiltersCombiner([]))
        self.assertEqual(filters, filters_2)
        # Change sub
        filters_2 = filtersFromDict(filter_def)
        filters_2.filters.append(FiltersCombiner([Filter("==", "A")]))
        self.assertNotEqual(filters, filters_2)
        # Change sub filter
        filters_2 = filtersFromDict(filter_def)
        filters_2.filters[1].values = 15
        self.assertNotEqual(filters, filters_2)
        # Change name
        filters_2 = filtersFromDict(filter_def)
        filters_2.name = "other filter"
        self.assertNotEqual(filters, filters_2)

    def testOperatorAnd(self):
        filters = filtersFromDict({
            "class": "FiltersCombiner",
            "operator": "and",
            "filters": [
                {
                    "class": "Filter",
                    "getter": "age",
                    "operator": "<",
                    "values": 10
                },
                {
                    "action": "exclude",
                    "class": "EmptyIterFilter",
                    "getter": "group"
                }
            ]
        })
        observed = [filters.eval(elt) for elt in self.data]
        self.assertEqual(observed, [True, False, False, False, False])

    def testOperatorOr(self):
        filters = filtersFromDict({
            "class": "FiltersCombiner",
            "operator": "or",
            "filters": [
                {
                    "class": "Filter",
                    "getter": "age",
                    "operator": ">",
                    "values": 40
                },
                {
                    "action": "exclude",
                    "class": "EmptyIterFilter",
                    "getter": "group"
                }
            ]
        })
        observed = [filters.eval(elt) for elt in self.data]
        self.assertEqual(observed, [True, True, False, False, True])

    def testOperatorXor(self):
        filters = filtersFromDict({
            "class": "FiltersCombiner",
            "operator": "xor",
            "filters": [
                {
                    "class": "Filter",
                    "getter": "age",
                    "operator": ">",
                    "values": 30
                },
                {
                    "class": "EmptyIterFilter",
                    "getter": "group"
                }
            ]
        })
        observed = [filters.eval(elt) for elt in self.data]
        self.assertEqual(observed, [False, False, True, False, False])


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    unittest.main()
