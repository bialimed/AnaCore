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
from anacore.filters import Filter, FiltersCombiner

TEST_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PACKAGE_DIR = os.path.dirname(TEST_DIR)
sys.path.append(PACKAGE_DIR)


#######
# must be implemented
#
class EmptyArrayFilter():
    def __init__(self, getter=None, name=None, description=None, action="select", safe=True):
        pass

    def eval(self, data):
        pass


def filtersFromDict(dict, safeMode=True):
    # Add 'safeMode' to existing filtersFromDict method
    pass

#######


data1 = [
  {'age': 8, 'treatment': "placebo", 'checked': False, 'group': ["A", "B"]},
  {'age': 9, 'treatment': "20ng ip", 'checked': True, 'group': ["C"]},
  {'age': 11, 'treatment': "placebo", 'checked': False, 'group': ["A", "D"]},
  {'age': 25, 'treatment': "20ng ip", 'checked': True, 'group': ["B"]},
  {'age': 30, 'treatment': "20ng pc", 'checked': False, 'group': ["C", "B"]},
  {'age': 33, 'treatment': "20ng pc", 'checked': True, 'group': ["D", "A"]},
  {'age': 30, 'treatment': "placebo", 'checked': False, 'group': ["A", "B"]},
]
data2 = [
  {},
  {'age': None, 'class': None, 'kind': None, 'checked': None},
  {'age': 0, 'class': [], 'kind': "", 'checked': False},
  {'age': 10, 'class': [None, "A"], 'kind': " ", 'checked': False},
  {'age': 1, 'class': ["B"], 'kind': "\t", 'checked': True},
  {'age': 1, 'class': ["", "C"], 'kind': "foobar", 'checked': True}
]


########################################################################
#
# Test cases
#
########################################################################
class TestFilter(unittest.TestCase):

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
                    {'class': "EmptyArrayFilter"},
                    EmptyArrayFilter(action="select")
                ],
                "case 1": [
                    {'class': "EmptyArrayFilter", 'getter': "myKey", 'name': "aName"},
                    EmptyArrayFilter("myKey", "aName", action="select")
                ],
                "case 2": [
                    {'class': "EmptyArrayFilter", 'getter': "myKey", 'action': "exclude", 'description': "aDescription"},
                    EmptyArrayFilter("myKey", None, "aDescription", "exclude")
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
                            {'class': "EmptyArrayFilter", 'getter': "myKey2"}
                        ],
                        'description': "aDescription",
                    },
                    FiltersCombiner(
                        [
                            Filter("=", "foo", "myKey"),
                            EmptyArrayFilter("myKey2", action="select")
                        ],
                        "or", description="aDescription")
                ],
            },
        }
        case_list = flattenDict(case_dict)
        for caseName, caseValues in case_list.items():
            with self.subTest(caseName):
                expected = caseValues[1]
                if expected is Exception:
                    try:
                        filtersFromDict(caseValues[0])
                    except Exception:
                        continue  # fail as expected
                    raise AssertionError("filtersFromDict does not fail")
                else:
                    result = filtersFromDict(caseValues[0])
                    self.assertEquals(result, expected)

    def testFromDict_unsafeMode(self):
        case_dict = {
                "accept simple": [{'class': "Filter", 'operator': "=", 'values': "foo", 'getter': "myKey"}, True],
                "reject bar()": [{'class': "Filter", 'operator': "=", 'values': "foo", 'getter': "bar()"}, False],
                "reject function getter": [
                    {'class': "Filter", 'operator': "=", 'values': "foo", 'getter': lambda x: x.bar}, False],
                "reject nested bar() ": [
                    {
                        'class': "FiltersCombiner",
                        'filters': [{'class': "Filter", 'operator': "=", 'values': "foo", 'getter': "bar()"}]
                    }, False
                ]
            }
        case_list = flattenDict(case_dict)
        safeMode = False
        for caseName, caseValues in case_list.items():
            with self.subTest(caseName):
                shouldSucceed = caseValues[1]
                if shouldSucceed:
                    filtersFromDict(caseValues[0], safeMode)
                else:
                    try:
                        filtersFromDict(caseValues[0], safeMode)
                    except Exception:
                        continue  # fail as expected
                    raise AssertionError("filtersFromDict (unsafe mode) does not fail")

    def testBuildEmptyArrayFilter(self):
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
        for caseName, caseValues in case_list.items():
            with self.subTest(caseName):
                args = caseValues[0]
                shouldSucceed = caseValues[1]
                if shouldSucceed:
                    EmptyArrayFilter(*args)
                else:
                    try:
                        EmptyArrayFilter(*args)
                    except Exception:
                        continue  # fail as expected
                    raise AssertionError("EmptyArrayFilter does not fail")

    def testEmptyArrayFilter(self):
        case_dict = {
                "None": [None, True],
                "empty-array": [[], True],
                "non-empty-array": [['foo'], False],
                "string": ['foo', Exception],
                "number": [42, Exception],
                "boolean": [True, Exception],
                "object": [{}, Exception],
            }

        filter = EmptyArrayFilter()
        case_list = flattenDict(case_dict)
        for caseName, caseValues in case_list.items():
            with self.subTest(caseName):
                data = caseValues[0]
                expected = caseValues[1]
                if expected is Exception:
                    try:
                        filter.eval(data)
                    except Exception:
                        continue  # fail as expected
                    raise AssertionError("EmptyArrayFilter does not fail")
                else:
                    result = filter.eval(data)
                    self.assertEqual(result, expected)

    def testEmptyArrayFilterWithGetter(self):
        filter = EmptyArrayFilter("class")
        expected = [True, True, True, False, False, False]
        observed = [filter.eval(curr) for curr in data2]
        self.assertEqual(expected, observed)

    def testSimpleFilter(self):
        case_dict = {
            "agg": {
                "eq": {
                    "nb-1": [
                        # Select patients with at least the group A
                        Filter("=", "A", "group", "nb:1"),
                        data1,
                        [True, False, True, False, False, True, True]
                    ],
                    "nb-1 + None": [
                        Filter("==", "A", "class", "nb:1"),
                        data2,
                        [False, False, False, True, False, False]
                    ],
                    "nb-1 ref=None": [
                        Filter("==", None, "class", "nb:1"),
                        data2,
                        [False, False, False, True, False, False]
                    ],
                    "ratio-100": [
                        # Select patients with 100% of groups are equals to C
                        Filter("=", "C", "group", "ratio:1"),
                        data1,
                        [False, True, False, False, False, False, False]
                    ],
                    "ratio-50": [
                        # Select patients with 50% of groups are equals to C
                        Filter("=", "C", "group", "ratio:0.5"),
                        data1,
                        [False, True, False, False, True, False, False]
                    ]
                },
                "ne": {
                    "nb-1 + None": [
                        Filter("<>", "B", "class", "nb:1"),
                        data2,
                        [False, False, False, True, False, True]
                    ],
                    "ratio-100 + Noneish": [
                        Filter("<>", "A", "class", "ratio:1"),
                        data2,
                        [True, True, True, False, True, True]
                    ]
                },
                "in": {
                    "nb-1": [
                        Filter("in", ["A", "C"], "group", "nb:1"),  # Select patients with at least A or C in groups
                        data1,
                        [True, True, True, False, True, True, True]
                    ],
                    "nb-1 ref=[None,string]": [
                        Filter("in", [None, "B"], "class", "nb:1"),
                        data2,
                        [False, False, False, True, True, False]
                    ],
                    "ratio-100": [
                        Filter("in", ["B", "C"], "group", "ratio:1"), # Select patients where groups are B and/or C
                        data1,
                        [False, True, False, True, True, False, False]
                    ],
                    "ratio-100 ref=[None,string]": [
                        Filter("in", [None, "A"], "class", "ratio:1"),
                        data2,
                        [True, True, True, True, False, False]
                    ]
                },
                "not-in": {
                    "nb-1": [
                        # Select patients with at least one group not in B and C
                        Filter("not in", ["B", "C"], "group", "nb:1"),
                        data1,
                        [True, False, True, False, False, True, True]
                    ],
                    "nb-1 + Noneish": [
                        Filter("not in", ["A", "B"], "class", "nb:1"),
                        data2,
                        [False, False, False, True, False, True]
                    ],
                    "ratio-100": [
                        # Select patients where neither group is B or C
                        Filter("not in", ["B", "C"], "group", "ratio:1"),
                        data1,
                        [False, False, True, False, False, True, False]
                    ],
                    "ratio-100 + None": [
                        Filter("not in", ["A", None], "class", "ratio:1"),
                        data2,
                        [True, True, True, False, True, True]
                    ]
                },
                "empty": {
                    "nb-1": [
                        Filter("empty", None, "class", "nb:1"),
                        data2,
                        [False, False, False, True, False, True]
                    ],
                    "ratio-50": [
                        Filter("empty", None, "class", "ratio:.5"),
                        data2,
                        [True, True, True, True, False, True]
                    ]
                }
            },
            "exclude": {
                "lt": [
                    Filter("<", 25, "age", action="exclude"),
                    data1, 
                    [False, False, False, True, True, True, True]
                ],
                "eq": [
                    Filter("=", "20ng pc", "treatment", action="exclude"),
                    data1, 
                    [True, True, True, True, False, False, True]
                ],
            }
        }
        case_list = flattenDict(case_dict)
        for caseName, caseValues in case_list.items():
            with self.subTest(caseName):
                filter = caseValues[0]
                expected = caseValues[2]
                observed = [filter.eval(curr) for curr in caseValues[1]]
                self.assertEqual(expected, observed)

    def testFiltersCombiner(self):
        case_dict = {
            "and": {
                "0 filter": [
                    FiltersCombiner([], "and"),
                    data1,
                    [True, True, True, True, True, True, True]
                ],
                "1 filter": [
                    FiltersCombiner([
                        Filter("<=", 25, "age"),
                    ], "and"),
                    data1,
                    [True, True, True, True, False, False, False]
                ],
                "2 filter": [
                    FiltersCombiner([
                        Filter(">", 10, "age"),
                        Filter("<=", 25, "age"),
                    ], "and"),
                    data1,
                    [False, False, True, True, False, False, False]
                ],
                "4 filter": [
                    FiltersCombiner([
                        Filter(">=", 9, "age"),
                        Filter("<=", 30, "age"),
                        Filter("!=", "20ng pc", "treatment"),
                        Filter("==", True, "checked"),
                    ], "and"),
                    data1,
                    [False, True, False, True, False, False, False]
                ]
            },
            "or": {
                "0 filter": [
                    FiltersCombiner([], "or"),
                    data1,
                    [True, True, True, True, True, True, True]
                ],
                "1 filter": [
                    FiltersCombiner([
                        Filter("<=", 25, "age"),
                    ], "or"),
                    data1,
                    [True, True, True, True, False, False, False]
                ],
                "2 filter": [
                    FiltersCombiner([
                        Filter("<=", 25, "age"),
                        Filter(">", 31, "age"),
                    ], "or"),
                    data1,
                    [True, True, True, True, False, True, False]
                ]
            },
            "xor": {
                "0 filter": [
                    # !! TODO : xor() = True ?
                    FiltersCombiner([], "xor"),
                    data1,
                    [True, True, True, True, True, True, True]
                ],
                "1 filter": [
                    FiltersCombiner([
                        Filter("<=", 25, "age"),
                    ], "xor"),
                    data1,
                    [True, True, True, True, False, False, False]
                ],
                "2 filter": [
                    FiltersCombiner([
                        Filter("<=", 25, "age"),
                        Filter("==", "placebo", "treatment"),
                    ], "xor"),
                    data1,
                    [False, True, False, True, False, False, True]
                ],
                "4 filter": [
                    FiltersCombiner([
                        Filter("<=", 25, "age"),
                        Filter("==", "placebo", "treatment"),
                        Filter("==", False, "checked"),
                        Filter("=", "A", "group", "nb:1")
                    ], "xor"),
                    data1,
                    [False, True, False, True, True, True, False]
                ]
            }
        }
        case_list = flattenDict(case_dict)
        for caseName, caseValues in case_list.items():
            with self.subTest(caseName):
                filter = caseValues[0]
                expected = caseValues[2]
                observed = [filter.eval(curr) for curr in caseValues[1]]
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
                "not Filter or FilterCombiner": (["foo"]),
            }
        }

        case_list = flattenDict(case_dict)
        for caseName, args in case_list.items():
            with self.subTest(caseName):
                try:
                    FiltersCombiner(*args)
                except Exception:
                    continue  # fail as expected
                raise AssertionError("FiltersCombiner build does not fail")
