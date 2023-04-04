#!/usr/bin/env python3

__author__ = 'Rémi THEVENOUX'
__copyright__ = 'Copyright (C) 2023 CHU Toulouse'
__license__ = 'GNU General Public License'
__version__ = '1.0.0'
__status__ = 'prod'

import os
import sys
import unittest

TEST_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PACKAGE_DIR = os.path.dirname(TEST_DIR)
sys.path.append(PACKAGE_DIR)


class GetterPath:
    @staticmethod
    def parse(input, safe=True):
        """
        Return a 'getter chain', an array of simple getter description.

        Parse input if it is a string, call validate is it an array.
        If 'safe'=True, method-getter are not allowed
        Throws an exception if input is invalid.
        """
        pass

    @staticmethod
    def validate(getter_chain, safe=True):
        """
        Throw an error if getter_chain is invalid
        (e.g.: not an array of dict, missing/invalid key in dict)

        If 'safe'=True, method-getter are not allowed
        """
        pass

    @staticmethod
    def apply(data, getter_chain):
        """
        Apply the getter_chain on the given data
        Return the value or throw an error if getter_chain is invalid
        """
        pass

    @staticmethod
    def equals(getter1, getter2):
        """
        Call parse() on both values and return True if they match
        """


########################################################################
#
# Test cases
#
########################################################################
class TestFilter(unittest.TestCase):
    def setUp(self):
        pass

    def assertParseEqual(self, sequence, expected, safe):
        if safe is None:
            result = GetterPath.parse(sequence)
            self.assertEquals(expected, result)
        else:
            result = GetterPath.parse(sequence, safe)
            self.assertEquals(expected, result)

    def assertParseFail(self, sequence, safe):
        try:
            if safe is None:
                GetterPath.parse(sequence)
            else:
                GetterPath.parse(sequence, safe)
        except Exception:
            return  # fail as expected
        raise AssertionError("Parse does not fail")

    def assertValidationEquals(self, array, valid, safe):
        if valid:
            if safe is None:
                GetterPath.validate(array)
            else:
                GetterPath.validate(array, safe)
        else:
            try:
                if safe is None:
                    GetterPath.validate(array)
                else:
                    GetterPath.validate(array, safe)
            except Exception:
                return  # fail as expected
            raise AssertionError("unexpected validation")

    def assertApplyEqual(self, input, data, expected):
        getterChains = GetterPath.parse(input)
        result = GetterPath.apply(data, getterChains)
        self.assertEquals(expected, result)

    def assertApplyFail(self, input, data):
        getterChains = GetterPath.parse(input)
        try:
            GetterPath.apply(data, getterChains)
        except Exception:
            return  # fail as expected
        raise AssertionError("Apply does not fail")

    def assertPathEquals(self, path1, path2, expected):
        result = GetterPath.equals(path1, path2)
        self.assertEquals(expected, result)

    def flattenDict(self, dict, result=None, prefix=None):
        if result is None:
            result = {}

        for key, value in dict.items():
            long_key = key if prefix is None else f'{prefix}.{key}'
            if type(value) == list:
                result[long_key] = value
            else:
                self.flattenDict(value, result, long_key)

        return result

    def testParse(self):
        case_dict = {
            "nullish": {
              "none": [None, []],
              "empty string": ["", []],
              "empty array": [[], []],
            },
            "key": {
              "simple": ["foo", [{'kind': 'key', 'key': 'foo'}]],
              "space": [" ", [{'kind': 'key', 'key': ' '}]],
              "unicode": ["\u2F54", [{'kind': 'key', 'key': '\u2F54'}]],
              "*-starting-key": {# check no conflict with iterator getter '*'
                "single": ["*a", [{"kind": 'key', "key": '*a'}]],
                "double": ["**.foo", [{"kind": 'key', "key": '**'}, {"kind": 'key', "key": 'foo'}]],
              }
            },
            "iterable": {
              "simple": ["*", [{"kind": 'iterable'}]],
            },
            "method": {
              "simple name": ["foo()", [{"kind": 'method', "name": 'foo'}]],
              "snake case": ["foo_bar_v2()", [{"kind": 'method', "name": 'foo_bar_v2'}]],
              "camel case": ["fooBarV2()", [{"kind": 'method', "name": 'fooBarV2'}]],
              "exotic": {
                "special character": ["%$\";@!: =/*()", [{"kind": 'method', "name": '%$\";@!: =/*'}]],
                "space": [" foo(bar) ()", [{"kind": 'method', "name": ' foo(bar) '}]],
                "parenthesis only": ["()()", [{"kind": 'method', "name": '()'}]],
                "parenthesis in the middle": ["foo()bar()", [{"kind": 'method', "name": 'foo()bar'}]],
                "multiple parenthesis": ["foo(())()", [{"kind": 'method', "name": 'foo(())'}]],
              },
              "error": {
                "no name": ["()", Exception],
              }
            },
            "chaining": {
              "keys": ["foo.bar.baz", [
                {"kind": 'key', "key": 'foo'},
                {"kind": 'key', "key": 'bar'},
                {"kind": 'key', "key": 'baz'},
              ]],
              "iterables": ["*.*", [
                {"kind": 'iterable'},
                {"kind": 'iterable'},
              ]],
              "methods": ["foo().bar()", [
                {"kind": 'method', "name": 'foo'},
                {"kind": 'method', "name": 'bar'}
              ]],
              "mix": {
                "m.k.i.k.m": ["foo().bar.*.biz.fuz()", [
                  {"kind": 'method', "name": 'foo'},
                  {"kind": 'key', "key": 'bar'},
                  {"kind": 'iterable'},
                  {"kind": 'key', "key": 'biz'},
                  {"kind": 'method', "name": 'fuz'},
                ]],
                "i.m.i": ["*.foo().*", [
                  {"kind": 'iterable'},
                  {"kind": 'method', "name": 'foo'},
                  {"kind": 'iterable'},
                ]],
                "k.k.m.k": ["bar.biz.foo().baz", [
                  {"kind": 'key', "key": 'bar'},
                  {"kind": 'key', "key": 'biz'},
                  {"kind": 'method', "name": 'foo'},
                  {"kind": 'key', "key": 'baz'},
                ]],
              }
            },
            "unsafe mode": {
              "accept-key-getter": ["foo", [{"kind": 'key', "key": 'foo'}], False],
              "reject-method-geter": ["foo()", Exception, False],
            },
            "error": {
              "dot-alone": [".", Exception],
              "initial-dot": [".foo", Exception],
              "double-dot": ["foo..bar", Exception],
              "final-dot": {
                "after key": ["foo.", Exception],
                "after '*'": ["*.", Exception],
                "after 'foo()'": ["foo().", Exception],
              },
            }
        }

        case_list = self.flattenDict(case_dict)
        for caseName, caseValues in case_list.items():
            with self.subTest(caseName):
                safe = caseValues[2] if len(caseValues) > 2 else None
                if caseValues[1] is Exception:
                    self.assertParseFail(caseValues[0], safe)
                else:
                    self.assertParseEqual(caseValues[0], caseValues[1], safe)

    def testValidate(self):
        case_dict = {
          "not-an-array": ["foo", False],
          "token-without-kind": [[{"a": 42}], False],
          "invalid-kind": [[{"kind": "invalid"}], False],
          "key-getter": {
            "valid": [[{"kind": "key", "key": "aKey"}], True],
            "missing key": [[{"kind": "key"}], False],
            "non-string key": [[{"kind": "key", "key": []}], False],
          },
          "iterable-getter": {
            "valid": [[{"kind": "iterable"}], True],
          },
          "method-getter": {
            "valid": {
              "without args": [[{"kind": "method", "name": "aName"}], True],  # arguments is optionnal
              "with args": [[{"kind": "method", "name": "aName", "arguments": [True, 'foo']}], True],
            },
            "missing name": [[{"kind": "method"}], False],
            "non-string name": [[{"kind": "method", "name": []}], False],
            "empty name": [[{"kind": "method", "name": ""}], False],
            "non-array args": [[{"kind": "method", "name": "aName", "arguments": True}], False],
          },
          "unsafe mode": {
            "accept-without-method": [[{"kind": "key", "key": "keyA"}, {"kind": "iterable", "key": "keyB"}], True, False],
            "reject-method": [[{"kind": "key", "key": "keyA"}, {"kind": "method", "name": "aName"}], False, False]
          }
        }

        case_list = self.flattenDict(case_dict)
        for caseName, caseValues in case_list.items():
            with self.subTest(caseName):
                safe = caseValues[2] if len(caseValues) > 2 else None
                self.assertValidationEquals(caseValues[0], caseValues[1], safe)

    def testEquals(self):
        case_dict = {
            "single-key-equals": ["foo", "foo", True],
            "single-key-differs-1": ["foo", "bar", False],
            "single-key-differs-2": ["foo", "*", False],
            "sub-key-equals": ["foo.bar", "foo.bar", True],
            "sub-key-differs": ["foo.bar", "fuz", False],
            "string-vs-array": ["foo.bar", [{'kind': 'key', 'key': 'foo'}, {'kind': 'key', 'key': 'bar'}], True],
            "array-vs-array": [
                [{'kind': 'key', 'key': 'foo'}, {'kind': 'key', 'key': 'bar'}],
                [{'kind': 'key', 'key': 'foo'}, {'kind': 'key', 'key': 'bar'}],
                True],
            "iterable": [[{'kind': 'key', 'key': 'foo'}, {'kind': 'iterable'}], "foo.*", True],
            "method": {
                "simple": [[{'kind': 'method', 'name': 'foo'}], "foo()", True],
                "args-equals": [
                    [{'kind': 'method', 'name': 'foo', 'arguments': ['a', 1]}],
                    [{'kind': 'method', 'name': 'foo', 'arguments': ['a', 1]}],
                    True],
                "args-differs-1": [
                    [{'kind': 'method', 'name': 'foo', 'arguments': ['a', 1]}],
                    [{'kind': 'method', 'name': 'foo'}],
                    False],
                "args-differs-2": [
                    [{'kind': 'method', 'name': 'foo', 'arguments': ['a', 1]}],
                    [{'kind': 'method', 'name': 'foo', 'arguments': ['a']}],
                    False],
                "args-differs-3": [
                    [{'kind': 'method', 'name': 'foo', 'arguments': ['a', 1]}],
                    [{'kind': 'method', 'name': 'foo', 'arguments': ['a', 2]}],
                    False],
                "comlpex-args-equals": [
                    [{'kind': 'method', 'name': 'foo', 'arguments': ['a', {'b': 'b'}, [1, 2]]}],
                    [{'kind': 'method', 'name': 'foo', 'arguments': ['a', {'b': 'b'}, [1, 2]]}],
                    True],
            }
        }

        case_list = self.flattenDict(case_dict)
        for caseName, caseValues in case_list.items():
            with self.subTest(caseName):
                self.assertPathEquals(caseValues[0], caseValues[1], caseValues[2])

    def testApply(self):
        myObject = {
          'aInt': 2,
          'aString': "Hello",
          'method1': lambda self: self.aInt,
          'method2': lambda self, value: value * self.aInt
        }

        case_dict = {
          "empty": {
            "none": ['', None, None],
            "string": ['', 'value', 'value'],
            "number": ['', 42, 42],
            "boolean": {
              "True": ['', True, True],
              "False": ['', False, False],
            },
            "array": ['', ['ax', 'bx'], ['ax', 'bx']],
            "object": ['', {'k': 'v'}, {'k': 'v'}],
          },
          "key": {
            "none": {
              "string-key": ["foo", None, None],
              "numeric-key": ["0", None, None],
            },
            "number": {
              "integer": ["foo", 42, Exception],
              "float": ["foo", 3.14, Exception],
            },
            "string": {
              "simple": ["foo", "hello", Exception],
              "numeric-key": ["0", "hello", Exception],  # Can not access n-th character
            },
            "boolean": {
              "True": ["foo", True, Exception],
              "False": ["foo", False, Exception],
            },
            "object": {
              "none-value": ["k", {'x': 'x', 'k': None}, None],
              "string-value": ["k", {'k': 'foo', 'x': 'x'}, "foo"],
              "number-value": ["k", {'k': 42, 'x': 'x'}, 42],
              "True-value": ["k", {'k': True, 'x': 'x'}, True],
              "False-value": ["k", {'k': False, 'x': 'x'}, False],
              "array-value": ["k", {'k': ["a", "b"], 'x': 'x'}, ["a", "b"]],
              "object-value": ["k", {'k': {'sub': 'v'}, 'x': 'x'}, {'sub': 'v'}],
              "numeric-key 1": ["1", {'1': 'a', 2: 'b'}, 'a'],
              "numeric-key 2": ["2", {'1': 'a', 2: 'b'}, 'b'],
            },
            "array": {
              "index '0'": ["arr.0", {'arr': ['a', 'b', 'c', 'd']}, 'a'],
              "index '2'": ["arr.2", {'arr': ['a', 'b', 'c', 'd']}, 'c'],
              "negative index": ["arr.-1", {'arr': ['a', 'b', 'c', 'd']}, Exception],
              "out of bound index": ["arr.20", {'arr': ['a', 'b', 'c', 'd']}, Exception],
              "non-numeric index": ["arr.foo", {'arr': ['a', 'b', 'c', 'd']}, Exception],
            },
            "chaining": {
              "chaining": ["obj.sub", {'obj': {'sub': 'val'}}, 'val'],
              "invalid key 1": ["foo.bar", {'foo': {'biz': 42}}, None],
              "invalid key 2": ["foo.bar", {'foo': None}, None],
              "invalid key 3": ["foo.bar", {}, None],
            }
          },
          "iterable": {
            "none": ["*", None, []],
            "string": ["*", "foo", Exception],
            "number": {
              "integer": ["*", 42, Exception],
              "float": ["*", 3.14, Exception],
            },
            "boolean": {
              "True": ["*", True, Exception],
              "False": ["*", False, Exception],
            },
            "array": {
              "'*'": ["*", ['ax', 'bx'], ['ax', 'bx']],
              "'key.*'": ["arr.*", {'arr': ['ax', 'bx']}, ['ax', 'bx']],
              "'*.subkey'": ["*.age", [{'age': 10}, {'age': 15}], [10, 15]],
              "'key.*.subkey'": ["key.*.age", {'key': [{'age': 10}, {'age': 15}]}, [10, 15]],
            },
            "object": {
              "'*'": ["*", {'k1': 'ax', 'k2': 'bx'}, ['ax', 'bx']],
              "'key.*'": ["key.*", {'key': {'k1': 'ax', 'k2': 'bx'}}, ['ax', 'bx']],
              "'*.subkey'": ["*.age", {'k1': {'age': 10}, 'k2': {'age': 15}}, [10, 15]],
              "'key.*.subkey'": ["key.*.age", {'key': {'k1': {'age': 10}, 'k2': {'age': 15}}}, [10, 15]]
            },
            "multiple": {
              "array+array": ["aa.*.bb.*.cc", {
                'aa': [
                  {'bb': [{'cc': "a1"}, {'cc': "a2"}]},
                  {'bb': [{'cc': "b1"}, {'cc': "b2"}]}
                ]
              }, ["a1", "a2", "b1", "b2"]],

              "object+array": ["*.*.bb", {
                'a1': [{'bb': "a1"}, {'bb': "a2"}],
                'a2': [{'bb': "b1"}, {'bb': "b2"}]
              }, ["a1", "a2", "b1", "b2"]],

              "object+object": ["*.*", {
                'aa': {'a1': "a1", 'a2': "a2"},
                'bb': {'b1': "b1", 'b2': "b2"}
              }, ["a1", "a2", "b1", "b2"]],

              "array+object": ["*.*", [
                {'a1': "a1", 'a2': "a2"},
                {'b1': "b1", 'b2': "b2"}
              ], ["a1", "a2", "b1", "b2"]],
            }
          },
          "method": {
            "none": ['method()', None, None],
            "string": ['toUpperCase()', 'hello', 'HELLO'],
            "boolean": ['toString()', True, 'True'],
            "number": {
              "without arg": ['toString()', 233, '233'],
              "with arg": [[{'kind': 'method', 'name': 'toString', 'arguments': [16]}], 233, 'e9'],
            },
            "array": ['toString()', ['a', 'b'], "a,b"],
            "object": {
              "no arg array": ['method1()', myObject, 2],
              "empty arg array": [[{'kind': 'method', 'name': 'method1', 'arguments': []}], myObject, 2],
              "1 param": [[{'kind': 'method', 'name': 'method2', 'arguments': [3]}], myObject, 6],
            },
            "invalid": {
              "method-does-not-exist": ["i_am_not_a_method()", myObject, Exception],
              "not-a-method": ["aInt()", myObject, Exception],
            }
          }
        }

        case_list = self.flattenDict(case_dict)
        for caseName, caseValues in case_list.items():
            with self.subTest(caseName):
                if caseValues[2] is Exception:
                    self.assertApplyFail(caseValues[0], caseValues[1])
                else:
                    self.assertApplyEqual(caseValues[0], caseValues[1], caseValues[2])


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    unittest.main()
