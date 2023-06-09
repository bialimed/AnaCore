#!/usr/bin/env python3

__author__ = 'Rémi THEVENOUX'
__copyright__ = 'Copyright (C) 2023 CHU Toulouse'
__license__ = 'GNU General Public License'
__version__ = '1.0.0'
__status__ = 'prod'

import datetime
import os
import sys
import unittest

FILTERS_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(FILTERS_DIR)
TEST_DIR = os.path.dirname(FILTERS_DIR)
PACKAGE_DIR = os.path.dirname(TEST_DIR)
sys.path.append(PACKAGE_DIR)

from anacore.filters import Filter
from dict_util import flattenDict


########################################################################
#
# FUNCTIONS
#
########################################################################
class EvalFunction:
    @staticmethod
    def build(operator, operand):
        f = Filter(operator, operand)
        return lambda val, ref: f._evalFct(val, ref)


########################################################################
#
# Test cases
#
########################################################################
class TestEvalFunction(unittest.TestCase):

    def applyEvalTestSuite(self, case_dict):
        case_list = flattenDict(case_dict)
        for name, (val, operator, ref, expected) in case_list.items():
            with self.subTest(name):
                eval_fct = EvalFunction.build(operator, ref)
                if expected is Exception:
                    try:
                        eval_fct(val, ref)
                    except Exception:
                        return  # fail as expected
                    raise AssertionError("EvalFunction does not fail")
                else:
                    observed = eval_fct(val, ref)
                    self.assertEqual(observed, expected)

    def applyBuildTestSuite(self, case_dict):
        case_list = flattenDict(case_dict)
        for name, (operator, ref, expected) in case_list.items():
            with self.subTest(name):
                if expected:  # should succeed
                    EvalFunction.build(operator, ref)
                else:  # should failed
                    try:
                        EvalFunction.build(operator, ref)
                    except Exception:
                        return  # fail as expected
                    raise AssertionError("EvalFunction build does not fail")

    def orderOperatorBuildTests(self, operator):
        """
        Return a tests suite for order operators ( '<', '<=', '>=' and '>' )
        """
        return {
                # Accept : String, Number
                "string": [operator, "hello", True],
                "number": [operator, 42, True],
                # Reject : Boolean and None, Array, Object, Date
                "boolean": [operator, False, False],
                "None": [operator, None, False],
                "array": [operator, ["a", "b"], False],
                "object": [operator, {'a': "b"}, False],
                "date": [operator, datetime.datetime(2020, 2, 2), False],
            }

    def orderOperatorEvalTests(self, operator, less, equal, greater):
        """
        Return a tests suite for order operators ( '<', '<=', '>=' and '>' )
        """
        return {
            "number": {
                "less": [19, operator, 42, less],
                "equal": [42, operator, 42, equal],
                "greater": [50, operator, 42, greater],
                "less float": [2.14, operator, 3.14, less],
                "equal float": [3.14, operator, 3.14, equal],
                "greater float": [4.14, operator, 3.14, greater],
                "None": [None, operator, 42, False],
                # fail if first operand is not a number
                "invalid-first-operand": {
                    "array": [[4, 2], operator, 42, Exception],
                    "boolean": [True, operator, 42, Exception],
                    "object": [{}, operator, 42, Exception],
                    "date": [datetime.datetime(2020, 2, 2), operator, 42, Exception],
                    "string": ["foo", operator, 42, Exception]  # no string (left) with number (rigth)
                }
            },
            "string": {
                "less": ['aaa', operator, 'abc', less],
                "equal": ['abc', operator, 'abc', equal],
                "greater": ['bbb', operator, 'abc', greater],
                "case sensitive": ['ABC', operator, 'aaa', less],  # case sensitive
                "None": [None, operator, 'aa', False],
                # fail if first operand is not a string
                "invalid-first-operand": {
                    "array": [[4, 2], operator, "foo", Exception],
                    "boolean": [True, operator, "foo", Exception],
                    "object": [{}, operator, "foo", Exception],
                    "date": [datetime.datetime(2020, 2, 2), operator, "foo", Exception],
                    "number": [42, operator, "foo", Exception],  # no number (left) with string (rigth)
                },
            }
        }

    def equalityOperatorBuildTests(self, operator):
        """
        Return a test suite , for equality operators ( '==' and '!=' )
        """
        return {
            # Accept : String, Number, Boolean and None
            "string": [operator, "hello", True],
            "number": [operator, 42, True],
            "boolean": [operator, False, True],
            "none": [operator, None, True],
            # Reject : Array, Object, Date
            # "array": [operator, ["a", "b"], False],
            # "object": [operator, {'a': "b"}, False],
            # "date": [operator, datetime.date(2020, 2, 2), False],
        }

    def equalityOperatorEvalTests(self, operator, equal):
        """
        Return a test suite , for equality operators ( '==' and '!=' )
        """
        return {
            "string": {
                "equal": ['foo', operator, 'foo', equal],
                "not equal": ['bar', operator, 'foo', not equal],
                "case-sensitive": ['FOO', operator, 'foo', not equal],  # case sensitive
                "empty": {
                    "empty": ['', operator, '', equal],
                    "None": [None, operator, '', not equal],
                },
                "None": [None, operator, 'foo', not equal],
                "array": [['a', 'b'], operator, 'ab', not equal],
            },
            "number": {
                "integer": {
                    "equal": [42, operator, 42, equal],
                    "not equal": [17, operator, 42, not equal],
                    "None": [None, operator, 42, not equal],
                    "string": ['42', operator, 42, not equal],
                },
                "float": {
                    "equal": [3.14, operator, 3.14, equal],
                    "not equal": [1.41, operator, 3.14, not equal],
                    "None": [None, operator, 3.14, not equal],
                    "string": ['3.14', operator, 3.14, not equal],
                },
            },
            "boolean": {
                "t/t": [True, operator, True, equal],
                "t/f": [True, operator, False, not equal],
                "f/t": [False, operator, True, not equal],
                "f/f": [False, operator, False, equal],
                "'True'": ['True', operator, True, not equal],
                "'False'": ['False', operator, False, not equal],
                "1": [1, operator, True, not equal],
                "0": [0, operator, False, not equal],
                "None/t": [None, operator, True, not equal],
                "None/f": [None, operator, False, not equal],
            },
            "None": {
                "None": [None, operator, None, equal]
            }
        }

    def inOperatorBuildTests(self, operator):
        """
        Return a tests suite for 'in' and 'not in' operators
        """
        return {
            # Accept : Array of String or Array of Number
            "array-string": [operator, ["a", "b"], True],
            "array-number": [operator, [4.2, 8], True],
            "array-with-None": [operator, ["a", "b", None], True],
            "array-mixed": [operator, [42, "foo", None], True],
            # Reject :
            "empty-array": [operator, [], False],
            "string": [operator, "hello", False],
            "number": [operator, 42, False],
            "boolean": [operator, False, False],
            "None": [operator, None, False],
            "object": [operator, {'a': "b"}, False],
            "date": [operator, datetime.datetime(2020, 2, 2), False],
        }

    def inOperatorEvalTests(self, operator, _in):
        """
        Return a tests suite for 'in' and 'not in' operators
        """
        return {
            "array": {
                "string-present": ["a", operator, ["a", "b"], _in],
                "string-absent-1": ["a", operator, ["ab", "cd"], not _in],
                "string-absent-2": ["ab", operator, ["a", "b"], not _in],
                "number-present": [4, operator, [2, 4, 6, 8], _in],
                "number-absent": [5, operator, [2, 4, 6, 8], not _in],
                "float-number-present": [4.2, operator, [4.2, 8.6], _in],
                "float-number-absent": [4.2, operator, [3.14, 8.6], not _in],
                "number-as-string": [4, operator, ["4", "a"], not _in],
                "None-present": [None, operator, ["a", None], _in],
                "None-absent": [None, operator, ["a", "b"], not _in],
            }
        }

    def testBuildEq(self):
        suite = self.equalityOperatorBuildTests('==')
        self.applyBuildTestSuite(suite)

    def testEvalEq(self):
        suite = self.equalityOperatorEvalTests('==', True)
        self.applyEvalTestSuite(suite)

    def testBuildNe(self):
        suite = self.equalityOperatorBuildTests('!=')
        self.applyBuildTestSuite(suite)

    def testEvalNe(self):
        suite = self.equalityOperatorEvalTests('!=', False)
        self.applyEvalTestSuite(suite)

    def testBuildLt(self):
        suite = self.orderOperatorBuildTests('<')
        self.applyBuildTestSuite(suite)

    def testEvalLt(self):
        suite = self.orderOperatorEvalTests('<', True, False, False)
        self.applyEvalTestSuite(suite)

    def testBuildLe(self):
        suite = self.orderOperatorBuildTests('<=')
        self.applyBuildTestSuite(suite)

    def testEvalLe(self):
        suite = self.orderOperatorEvalTests('<=', True, True, False)
        self.applyEvalTestSuite(suite)

    def testBuildGe(self):
        suite = self.orderOperatorBuildTests('>=')
        self.applyBuildTestSuite(suite)

    def testEvalGe(self):
        suite = self.orderOperatorEvalTests('>=', False, True, True)
        self.applyEvalTestSuite(suite)

    def testBuildGt(self):
        suite = self.orderOperatorBuildTests('>')
        self.applyBuildTestSuite(suite)

    def testEvalGt(self):
        suite = self.orderOperatorEvalTests('>', False, False, True)
        self.applyEvalTestSuite(suite)

    def testBuildIn(self):
        suite = self.inOperatorBuildTests('in')
        self.applyBuildTestSuite(suite)

    def testEvalIn(self):
        suite = self.inOperatorEvalTests('in', True)
        self.applyEvalTestSuite(suite)

    def testBuildNotIn(self):
        suite = self.inOperatorBuildTests('not in')
        self.applyBuildTestSuite(suite)

    def testEvalNotIn(self):
        suite = self.inOperatorEvalTests('not in', False)
        self.applyEvalTestSuite(suite)

    def testBuildContain(self):
        suite = {
            # Accept : String only
            "string": ["contains", "hello", True],
            "numeric-string": ["contains", "42", True],
            # Reject :  Number, Boolean, Null, Array, Object, Date
            "number": ["contains", 42, False],
            "boolean": ["contains", False, False],
            "None": ["contains", None, False],
            "array": ["contains", ["a", "b"], False],
            "object": ["contains", {'a': "b"}, False],
            "date": ["contains", datetime.datetime(2020, 2, 2), False],
        }
        self.applyBuildTestSuite(suite)

    def testEvalContain(self):
        suite = {
            "full match": ["Hello World", "contains", "Hello World", True],
            "case-sensitive": ["hello world", "contains", "Hello World", False],
            "partial match": ["abc def", "contains", "c d", True],
            "None": [None, "contains", "foo", False],
            "number": [4862, "contains", "86", True],  # number are convert to string
            # Test invalid first operand
            # "boolean": [True, "contains", "foo", Exception],  # CHANGE: No check for type
            # "array": [['a'], "contains", "foo", Exception],  # CHANGE: No check for type
            # "date": [datetime.datetime(2020, 2, 2), "contains", "foo", Exception],  # CHANGE: No check for type
        }
        self.applyEvalTestSuite(suite)

    def testBuildEmpty(self):
        suite = {
            # Accept : undefined and None
            "None": ["empty", None, True],
            # Reject : Any not-None value
            "string": ["empty", "hello", False],
            "number": ["empty", 42, False],
            "boolean": ["empty", False, False],
            "array": ["empty", ["a", "b"], False],
            "object": ["empty", {'a': "b"}, False],
            "date": ["empty", datetime.datetime(2020, 2, 2), False],
        }
        self.applyBuildTestSuite(suite)

    def testEvalEmpty(self):
        suite = {
            "None": [None, "empty", None, True],
            "string": {
                "empty": ["", "empty", None, True],
                "space": [" ", "empty", None, False],
                "tab": ["\t", "empty", None, False],
                "new line": ["\n", "empty", None, False],
                "no sapce": ["foo", "empty", None, False],
                "not only space": ["  a  ", "empty", None, False],
            },
            "number": {
                "0": [0, "empty", None, False],
                "42": [42, "empty", None, False],
            },
            "boolean": {
                "True": [True, "empty", None, False],
                "False": [False, "empty", None, False],
            },
            "array": {
                "empty": [[], "empty", None, False],
                "not empty": [['foo'], "empty", None, False],
            },
            "object": {
                "empty": [{}, "empty", None, False],
                "not empty": [{'foo': 'bar'}, "empty", None, False],
            },
            "date": {
                "empty": [datetime.datetime(2020, 2, 2), "empty", None, False],
            }
        }
        self.applyEvalTestSuite(suite)


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    unittest.main()
