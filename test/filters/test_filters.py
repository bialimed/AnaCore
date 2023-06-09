#!/usr/bin/env python3

__author__ = 'Rémi THEVENOUX'
__copyright__ = 'Copyright (C) 2023 CHU Toulouse'
__license__ = 'GNU General Public License'
__version__ = '1.0.0'
__status__ = 'prod'

import os
import sys
import unittest

FILTERS_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(FILTERS_DIR)
TEST_DIR = os.path.dirname(FILTERS_DIR)
PACKAGE_DIR = os.path.dirname(TEST_DIR)
sys.path.append(PACKAGE_DIR)

from anacore.filters import EmptyIterFilter, Filter, FiltersCombiner, filtersFromDict
from dict_util import flattenDict


########################################################################
#
# Test cases
#
########################################################################
class TestFilter(unittest.TestCase):
    def setUp(self):
        self.data1 = [
            {'age': 8, 'treatment': "placebo", 'checked': False, 'group': ["A", "B"]},
            {'age': 9, 'treatment': "20ng ip", 'checked': True, 'group': ["C"]},
            {'age': 11, 'treatment': "placebo", 'checked': False, 'group': ["A", "D"]},
            {'age': 25, 'treatment': "20ng ip", 'checked': True, 'group': ["B"]},
            {'age': 30, 'treatment': "20ng pc", 'checked': False, 'group': ["C", "B"]},
            {'age': 33, 'treatment': "20ng pc", 'checked': True, 'group': ["D", "A"]},
            {'age': 30, 'treatment': "placebo", 'checked': False, 'group': ["A", "B"]}
        ]
        self.data2 = [
            {},
            {'age': None, 'class': None, 'kind': None, 'checked': None},
            {'age': 0, 'class': [], 'kind': "", 'checked': False},
            {'age': 10, 'class': [None, "A"], 'kind': " ", 'checked': False},
            {'age': 1, 'class': ["B"], 'kind': "\t", 'checked': True},
            {'age': 1, 'class': ["", "C"], 'kind': "foobar", 'checked': True}
        ]

    def testFromDict(self):
        case_dict = {
            "invalid-class": {
                "None": [{'class': None, 'operator': "=", 'values': "foo", 'getter': "myKey"}, Exception],
                "invalid": [{'class': "invalid", 'operator': "=", 'values': "foo", 'getter': "myKey"}, Exception],
            },
            "simple-filter": {
                "class undefined": [
                    {'operator': "=", 'values': "foo", 'getter': "myKey"},
                    Filter("=", "foo", "myKey", action="select")
                ],
                "explicit class": [
                    {'class': "Filter", 'operator': "=", 'values': "foo", 'getter': "myKey", 'name': "aName"},
                    Filter("=", "foo", "myKey", name="aName", action="select")
                ],
                "no getter": [
                    {'class': "Filter", 'operator': "=", 'values': "foo", 'description': "aDescription"},
                    Filter("=", "foo", description="aDescription", action="select")
                ],
                "values None": [
                    {'operator': "=", 'values': None, 'getter': "myKey"},
                    Filter("=", None, "myKey", action="select")
                ],
                "values undefined": [{'operator': ">=", 'getter': "myKey"}, Exception],
                "aggregator": [
                    {'operator': "=", 'values': "foo", 'getter': "myKey", 'aggregator': "nb:2"},
                    Filter("=", "foo", "myKey", "nb:2", action="select")
                ],
                "action": [
                    {'operator': "=", 'values': "foo", 'getter': "myKey", 'action': "exclude"},
                    Filter("=", "foo", "myKey", action="exclude")
                ],
            },
            "empty-filter": {
                "minimal": [
                    {'class': "EmptyIterFilter"},
                    EmptyIterFilter(action="select")
                ],
                "case 1": [
                    {'class': "EmptyIterFilter", 'getter': "myKey", 'name': "aName"},
                    EmptyIterFilter("myKey", "aName", action="select")
                ],
                "case 2": [
                    {'class': "EmptyIterFilter", 'getter': "myKey", 'action': "exclude", 'description': "aDescription"},
                    EmptyIterFilter("myKey", None, "aDescription", "exclude")
                ],
            },
            "combiner": {
                "without operator": [
                    {
                        'class': "FiltersCombiner",
                        'filters': [{'operator': "=", 'values': "foo", 'getter': "myKey"}],
                        'name': "aName"
                    },
                    FiltersCombiner([Filter("=", "foo", "myKey")], "and", "aName")
                ],
                "or-operator": [
                    {
                        'class': "FiltersCombiner",
                        'operator': "or",
                        'filters': [
                            {'operator': "=", 'values': "foo", 'getter': "myKey"},
                            {'class': "EmptyIterFilter", 'getter': "myKey2"}
                        ],
                        'description': "aDescription",
                    },
                    FiltersCombiner(
                        [
                            Filter("=", "foo", "myKey"),
                            EmptyIterFilter("myKey2", action="select")
                        ],
                        "or", description="aDescription")
                ],
            },
        }
        case_list = flattenDict(case_dict)
        for name, (filter_hash, expected) in case_list.items():
            with self.subTest(name):
                if expected is Exception:
                    try:
                        filtersFromDict(filter_hash)
                    except Exception:
                        continue  # fail as expected
                    raise AssertionError("filtersFromDict does not fail")
                else:
                    result = filtersFromDict(filter_hash)
                    self.assertEqual(result, expected)

    # def testFromDict_unsafeMode(self):  # CHANGE: No check for safety
    #     case_dict = {
    #             "accept simple": [{'class': "Filter", 'operator': "=", 'values': "foo", 'getter': "myKey"}, True],
    #             "reject bar()": [{'class': "Filter", 'operator': "=", 'values': "foo", 'getter': "bar()"}, False],
    #             "reject function getter": [
    #                 {'class': "Filter", 'operator': "=", 'values': "foo", 'getter': lambda x: x.bar}, False
    #             ],
    #             "reject nested bar() ": [
    #                 {
    #                     'class': "FiltersCombiner",
    #                     'filters': [{'class': "Filter", 'operator': "=", 'values': "foo", 'getter': "bar()"}]
    #                 }, False
    #             ]
    #         }
    #     case_list = flattenDict(case_dict)
    #     safe_mode = False
    #     for name, (filter_hash, should_succeed) in case_list.items():
    #         with self.subTest(name):
    #             if should_succeed:
    #                 filtersFromDict(filter_hash)#, safe_mode)
    #             else:
    #                 try:
    #                     filtersFromDict(filter_hash)#, safe_mode)
    #                 except Exception:
    #                     continue  # fail as expected
    #                 raise AssertionError("filtersFromDict (unsafe mode) does not fail")

    def testBuildEmptyIterFilter(self):
        case_dict = {
            "simple": [(), True],
            "getter": [("foo"), True],
            "select": [("foo", None, None, "select"), True],
            "exclude": [("foo", None, None, "exclude"), True],
            "invalid": {
                "action": [("foo", None, None, "invalid"), False]
            }
        }
        case_list = flattenDict(case_dict)
        for name, (args, should_succeed) in case_list.items():
            with self.subTest(name):
                if should_succeed:
                    EmptyIterFilter(*args)
                else:
                    try:
                        EmptyIterFilter(*args)
                    except Exception:
                        continue  # fail as expected
                    raise AssertionError("EmptyIterFilter does not fail")

    def testEmptyIterFilter(self):
        case_dict = {
                "None": [None, True],
                "empty-array": [[], True],
                "non-empty-array": [['foo'], False],
                "string": ['foo', Exception],  # CHANGE: No check for type
                "number": [42, Exception],  # CHANGE: No check for type
                "boolean": [True, Exception],  # CHANGE: No check for type
                "object": [{}, True],  # CHANGE: dict is ok
            }
        filter = EmptyIterFilter()
        case_list = flattenDict(case_dict)
        for name, (data, expected) in case_list.items():
            with self.subTest(name):
                if expected is Exception:
                    try:
                        filter.eval(data)
                    except Exception:
                        continue  # fail as expected
                    raise AssertionError("EmptyIterFilter does not fail")
                else:
                    result = filter.eval(data)
                    self.assertEqual(result, expected)

    def testEmptyIterFilterWithGetter(self):
        filter = EmptyIterFilter("class")
        expected = [True, True, True, False, False, False]
        observed = [filter.eval(curr) for curr in self.data2]
        self.assertEqual(expected, observed)

    def testSimpleFilter(self):
        case_dict = {
            "agg": {
                "eq": {
                    "nb-1": [
                        # Select patients with at least the group A
                        Filter("=", "A", "group", "nb:1"),
                        self.data1,
                        [True, False, True, False, False, True, True]
                    ],
                    "nb-1 + None": [
                        Filter("==", "A", "class", "nb:1"),
                        self.data2,
                        [False, False, False, True, False, False]
                    ],
                    "nb-1 ref=None": [
                        Filter("==", None, "class", "nb:1"),
                        self.data2,
                        [False, False, False, True, False, False]
                    ],
                    "ratio-100": [
                        # Select patients with 100% of groups are equals to C
                        Filter("=", "C", "group", "ratio:1"),
                        self.data1,
                        [False, True, False, False, False, False, False]
                    ],
                    "ratio-50": [
                        # Select patients with 50% of groups are equals to C
                        Filter("=", "C", "group", "ratio:0.5"),
                        self.data1,
                        [False, True, False, False, True, False, False]
                    ]
                },
                "ne": {
                    "nb-1 + None": [
                        Filter("<>", "B", "class", "nb:1"),
                        self.data2,
                        [False, False, False, True, False, True]
                    ],
                    "ratio-100 + Noneish": [
                        Filter("<>", "A", "class", "ratio:1"),
                        self.data2,
                        [True, True, True, False, True, True]
                    ]
                },
                "in": {
                    "nb-1": [
                        Filter("in", ["A", "C"], "group", "nb:1"),  # Select patients with at least A or C in groups
                        self.data1,
                        [True, True, True, False, True, True, True]
                    ],
                    "nb-1 ref=[None,string]": [
                        Filter("in", [None, "B"], "class", "nb:1"),
                        self.data2,
                        [False, False, False, True, True, False]
                    ],
                    "ratio-100": [
                        Filter("in", ["B", "C"], "group", "ratio:1"), # Select patients where groups are B and/or C
                        self.data1,
                        [False, True, False, True, True, False, False]
                    ],
                    "ratio-100 ref=[None,string]": [
                        Filter("in", [None, "A"], "class", "ratio:1"),
                        self.data2,
                        [True, True, True, True, False, False]
                    ]
                },
                "not-in": {
                    "nb-1": [
                        # Select patients with at least one group not in B and C
                        Filter("not in", ["B", "C"], "group", "nb:1"),
                        self.data1,
                        [True, False, True, False, False, True, True]
                    ],
                    "nb-1 + Noneish": [
                        Filter("not in", ["A", "B"], "class", "nb:1"),
                        self.data2,
                        [False, False, False, True, False, True]
                    ],
                    "ratio-100": [
                        # Select patients where neither group is B or C
                        Filter("not in", ["B", "C"], "group", "ratio:1"),
                        self.data1,
                        [False, False, True, False, False, True, False]
                    ],
                    "ratio-100 + None": [
                        Filter("not in", ["A", None], "class", "ratio:1"),
                        self.data2,
                        [True, True, True, False, True, True]
                    ]
                },
                "empty": {
                    "nb-1": [
                        Filter("empty", None, "class", "nb:1"),
                        self.data2,
                        [False, False, False, True, False, True]
                    ],
                    "ratio-50": [
                        Filter("empty", None, "class", "ratio:.5"),
                        self.data2,
                        [True, True, True, True, False, True]
                    ]
                }
            },
            "exclude": {
                "lt": [
                    Filter("<", 25, "age", action="exclude"),
                    self.data1,
                    [False, False, False, True, True, True, True]
                ],
                "eq": [
                    Filter("=", "20ng pc", "treatment", action="exclude"),
                    self.data1,
                    [True, True, True, True, False, False, True]
                ],
            }
        }
        case_list = flattenDict(case_dict)
        for name, (filter, data, expected) in case_list.items():
            with self.subTest(name):
                observed = [filter.eval(curr) for curr in data]
                self.assertEqual(expected, observed)

    def testFiltersCombiner(self):
        case_dict = {
            "and": {
                "0 filter": [
                    FiltersCombiner([], "and"),
                    self.data1,
                    [True, True, True, True, True, True, True]
                ],
                "1 filter": [
                    FiltersCombiner([
                        Filter("<=", 25, "age"),
                    ], "and"),
                    self.data1,
                    [True, True, True, True, False, False, False]
                ],
                "2 filter": [
                    FiltersCombiner([
                        Filter(">", 10, "age"),
                        Filter("<=", 25, "age"),
                    ], "and"),
                    self.data1,
                    [False, False, True, True, False, False, False]
                ],
                "4 filter": [
                    FiltersCombiner([
                        Filter(">=", 9, "age"),
                        Filter("<=", 30, "age"),
                        Filter("!=", "20ng pc", "treatment"),
                        Filter("==", True, "checked"),
                    ], "and"),
                    self.data1,
                    [False, True, False, True, False, False, False]
                ]
            },
            "or": {
                "0 filter": [
                    FiltersCombiner([], "or"),
                    self.data1,
                    [True, True, True, True, True, True, True]
                ],
                "1 filter": [
                    FiltersCombiner([
                        Filter("<=", 25, "age"),
                    ], "or"),
                    self.data1,
                    [True, True, True, True, False, False, False]
                ],
                "2 filter": [
                    FiltersCombiner([
                        Filter("<=", 25, "age"),
                        Filter(">", 31, "age"),
                    ], "or"),
                    self.data1,
                    [True, True, True, True, False, True, False]
                ]
            },
            "xor": {
                "0 filter": [
                    # !! TODO : xor() = True ?
                    FiltersCombiner([], "xor"),
                    self.data1,
                    [True, True, True, True, True, True, True]
                ],
                "1 filter": [
                    FiltersCombiner([
                        Filter("<=", 25, "age"),
                    ], "xor"),
                    self.data1,
                    [True, True, True, True, False, False, False]
                ],
                "2 filter": [
                    FiltersCombiner([
                        Filter("<=", 25, "age"),
                        Filter("==", "placebo", "treatment"),
                    ], "xor"),
                    self.data1,
                    [False, True, False, True, False, False, True]
                ],
                "4 filter": [
                    FiltersCombiner([
                        Filter("<=", 25, "age"),
                        Filter("==", "placebo", "treatment"),
                        Filter("==", False, "checked"),
                        Filter("=", "A", "group", "nb:1")
                    ], "xor"),
                    self.data1,
                    [False, True, False, True, True, True, False]
                ]
            }
        }
        case_list = flattenDict(case_dict)
        for caseName, (filter, data, expected) in case_list.items():
            with self.subTest(caseName):
                observed = [filter.eval(curr) for curr in data]
                self.assertEqual(expected, observed)

    def testBuildFiltersCombiner(self):
        # All cases are invalid FiltersCombiner args
        case_dict = {
            "operator": ([], "foo"),
            "filters": {
                "undefined": (),
                "None": (None),
                "not array": {
                    "string": ("foo"),
                    "int": (42),
                },
                # "not Filter or FilterCombiner": (["foo"]),  # CHANGE: No check for type
            }
        }
        case_list = flattenDict(case_dict)
        for name, args in case_list.items():
            with self.subTest(name):
                try:
                    FiltersCombiner(*args)
                except Exception:
                    continue  # fail as expected
                raise AssertionError("FiltersCombiner build does not fail")


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    unittest.main()
