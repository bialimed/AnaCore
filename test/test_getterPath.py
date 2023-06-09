#!/usr/bin/env python3

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2023 CHU Toulouse'
__license__ = 'GNU General Public License'
__version__ = '1.0.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import os
import sys
import unittest

TEST_DIR = os.path.dirname(os.path.abspath(__file__))
PACKAGE_DIR = os.path.dirname(TEST_DIR)
sys.path.append(PACKAGE_DIR)

from anacore.getterPath import GetterPath, IterableToken, KeyToken, MethodToken


########################################################################
#
# FUNCTIONS
#
########################################################################
class TestGetterPath(unittest.TestCase):
    def testDict(self):
        test_cases = [
            {  # key => return record[key]
                "getter": GetterPath([KeyToken("id")]),
                "expected": [{"kind": "key", "key": "id"}]
            },
            {  # key.sub_key => return record[key][sub_key]
                "getter": GetterPath([
                    KeyToken("status"),
                    KeyToken("classifier1")
                ]),
                "expected": [{"kind": "key", "key": "status"}, {"kind": "key", "key": "classifier1"}]
            },
            {  # None => return record
                "getter": GetterPath(),
                "expected": []
            },
            {  # key.sub_key().sub_sub_key => return record[key].sub_key()[sub_sub_key]
                "getter": GetterPath([
                    KeyToken("obj"),
                    MethodToken("getInfo"),
                    KeyToken("B")
                ]),
                "expected": [
                    {"kind": "key", "key": "obj"},
                    {"kind": "method", "name": "getInfo"},
                    {"kind": "key", "key": "B"}
                ]
            },
            {  # key.*.sub_key => return record[key].map(lambda elt: elt[sub_key])
                "getter": GetterPath([
                    KeyToken("class"),
                    IterableToken(
                        GetterPath([
                            KeyToken("status")
                        ])
                    )
                ]),
                "expected": [
                    {"kind": "key", "key": "class"},
                    {"kind": "iterable"},
                    {"kind": "key", "key": "status"}
                ]
            },
            {  # key.1.sub_key => return record[key][1][sub_key]
                "getter": GetterPath([
                    KeyToken("class"),
                    KeyToken(1),
                    KeyToken("status")
                ]),
                "expected": [
                    {"kind": "key", "key": "class"},
                    {"kind": "key", "key": 1},
                    {"kind": "key", "key": "status"}
                ]
            },
            {  # key.method() => return record[key].method()
                "getter": GetterPath([
                    KeyToken("obj"),
                    MethodToken("getAD")
                ]),
                "expected": [
                    {"kind": "key", "key": "obj"},
                    {"kind": "method", "name": "getAD"}
                ]
            }
        ]
        for curr_test in test_cases:
            self.assertEqual(
                curr_test["getter"].toDict(), curr_test["expected"]
            )

    def testFromDict(self):
        test_cases = [
            {  # key => return record[key]
                "path": [{"kind": "key", "key": "id"}],
                "expected": GetterPath([KeyToken("id")])
            },
            {  # key.sub_key => return record[key][sub_key]
                "path": [
                    {"kind": "key", "key": "status"},
                    {"kind": "key", "key": "classifier1"}
                ],
                "expected": GetterPath([
                    KeyToken("status"),
                    KeyToken("classifier1")
                ])
            },
            {  # None => return record
                "path": None,
                "expected": GetterPath([])
            },
            {  # None => return record
                "path": [],
                "expected": GetterPath([])
            },
            {  # key.sub_key().sub_sub_key => return record[key].sub_key()[sub_sub_key]
                "path": [
                    {"kind": "key", "key": "obj"},
                    {"kind": "method", "name": "getInfo", "arguments": []},
                    {"kind": "key", "key": "B"}
                ],
                "expected": GetterPath([
                    KeyToken("obj"),
                    MethodToken("getInfo"),
                    KeyToken("B")
                ])
            },
            {  # key.*.*.sub_key()
                "path": [
                    {"kind": "key", "key": "shapes"},
                    {"kind": "iterable"},
                    {"kind": "iterable"},
                    {"kind": "method", "name": "scaleUp", "arguments": [2]}
                ],
                "expected": GetterPath([
                    KeyToken("shapes"),
                    IterableToken(
                        GetterPath([
                            IterableToken(
                                GetterPath([
                                    MethodToken("scaleUp", [2])
                                ])
                            )
                        ])
                    )
                ])
            },
            {  # key.*.sub_key => return record[key].map(lambda elt: elt[sub_key])
                "path": [
                    {"kind": "key", "key": "class"},
                    {"kind": "iterable"},
                    {"kind": "key", "key": "status"}
                ],
                "expected": GetterPath([
                    KeyToken("class"),
                    IterableToken(
                        GetterPath([
                            KeyToken("status")
                        ])
                    )
                ])
            },
            {  # key.1.sub_key => return record[key][1][sub_key]
                "path": [
                    {"kind": "key", "key": "class"},
                    {"kind": "key", "key": 1},
                    {"kind": "key", "key": "status"}
                ],
                "expected": GetterPath([
                    KeyToken("class"),
                    KeyToken(1),
                    KeyToken("status")
                ])
            },
            {  # key.method() => return record[key].method()
                "path": [
                    {"kind": "key", "key": "class"},
                    {"kind": "method", "name": "getAD", "arguments": []}
                ],
                "expected": GetterPath([
                    KeyToken("class"),
                    MethodToken("getAD")
                ])
            }
        ]
        for curr_test in test_cases:
            self.assertEqual(
                GetterPath.fromDict(curr_test["path"]),
                curr_test["expected"]
            )

    def testFromStr(self):
        test_cases = [
            {  # key => return record[key]
                "path": "id",
                "expected": [KeyToken("id")]},
            {  # key.sub_key => return record[key][sub_key]
                "path": "status.classifier1",
                "expected": [KeyToken("status"), KeyToken("classifier1")]
            },
            {  # None => return record
                "path": None,
                "expected": []
            },
            {  # None => return record
                "path": "",
                "expected": []
            },
            {  # key.sub_key().sub_sub_key => return record[key].sub_key()[sub_sub_key]
                "path": "obj.getInfo().B",
                "expected": [
                    KeyToken("obj"),
                    MethodToken("getInfo"),
                    KeyToken("B")
                ]
            },
            {  # key.*.*.sub_key()
                "path": "shapes.*.*.getArea()",
                "expected": [
                    KeyToken("shapes"),
                    IterableToken(
                        GetterPath([
                            IterableToken(
                                GetterPath([
                                    MethodToken("getArea", [])
                                ])
                            )
                        ])
                    )
                ]
            },
            {  # key.*.sub_key => return record[key].map(lambda elt: elt[sub_key])
                "path": "class.*.status",
                "expected": [
                    KeyToken("class"),
                    IterableToken(
                        GetterPath([
                            KeyToken("status")
                        ])
                    )
                ]
            },
            {  # key.1.sub_key => return record[key][1][sub_key]
                "path": "class.1.status",
                "expected": [
                    KeyToken("class"),
                    KeyToken("1"),
                    KeyToken("status")
                ]
            },
            {  # key.method() => return record[key].method()
                "path": "obj.getAD()",
                "expected": [
                    KeyToken("obj"),
                    MethodToken("getAD")
                ]
            }
        ]
        for curr_test in test_cases:
            self.assertEqual(
                GetterPath.fromStr(curr_test["path"]).tokens_chain,
                curr_test["expected"]
            )

    def testGet(self):
        class FakeObject:
            def __init__(self, data):
                self.data = data

            def getInfo(self, mult=None, sub=None):
                data = {k: v for k, v in self.data.items()}
                if mult:
                    data = {k: v * mult for k, v in data.items()}
                if sub:
                    data = {k: v - sub for k, v in data.items()}
                return data
        data = [
            {
                "id": 1,
                "status": {
                    "classifier1": "risk",
                    "classifier2": "undetermined"
                },
                "class": [
                    {"classifier": "1", "status": "risk", "score": 0.9},
                    {"classifier": "2", "status": "undetermined", "score": 0.2}
                ],
                "obj": FakeObject({"A": 10, "B": 22})
            },
            {
                "id": 2,
                "status": {
                    "classifier1": "risk",
                    "classifier2": "risk"
                },
                "class": [
                    {"classifier": "1", "status": "risk", "score": 1.0},
                    {"classifier": "2", "status": "risk", "score": 0.97}
                ],
                "obj": FakeObject({"B": 9, "C": 30})
            }
        ]
        test_cases = [
            {  # Unknown key
                "path": "nothing",
                "expected": [None, None]
            },
            {  # key => return record[key]
                "path": "id",
                "expected": [1, 2]
            },
            {  # Unknown key and continue token
                "path": "nothing.classifier1",
                "expected": [None, None]
            },
            {  # Unknown sub key
                "path": "status.classifier3",
                "expected": [None, None]
            },
            {  # key.sub_key => return record[key][sub_key]
                "path": "status.classifier1",
                "expected": ["risk", "risk"]
            },
            {  # None => return record
                "path": None,
                "expected": [data[0], data[1]]
            },
            {  # None => return record
                "path": "",
                "expected": [data[0], data[1]]
            },
            {  # key on object
                "path": "obj.data",
                "expected": [{"A": 10, "B": 22}, {"B": 9, "C": 30}]
            },
            {  # Unknown method and continue token
                "path": "obj.getUnknown().B",
                "expected": [None, None]
            },
            {  # key.sub_key(x, y) => return record[key].sub_key(x, y)
                "path": [{"kind": "key", "key": "obj"}, {"kind": "method", "name": "getInfo", "arguments": [2, 1]}],
                "expected": [{"A": 19, "B": 43}, {"B": 17, "C": 59}]
            },
            {  # key.sub_key().sub_sub_key => return record[key].sub_key()[sub_sub_key]
                "path": "obj.getInfo().B",
                "expected": [22, 9]
            },
            {  # Unknown after iter
                "path": "class.*.nothing",
                "expected": [[None, None], [None, None]]
            },
            {  # Unknown before iter
                "path": "nothing.*.nothing",
                "expected": [None, None]
            },
            {  # key.* in dict => return record[key].map(lambda elt: elt)
                "path": "status.*",
                "expected": [["risk", "undetermined"], ["risk", "risk"]]
            },
            {  # key.*.sub_key => return record[key].map(lambda elt: elt[sub_key])
                "path": "class.*.status",
                "expected": [["risk", "undetermined"], ["risk", "risk"]]
            },
            {  # key.1.sub_key => return record[key][1][sub_key]
                "path": "class.1.status",
                "expected": ["undetermined", "risk"]
            }
        ]
        for curr_test in test_cases:
            obs = list()
            for curr in data:
                if isinstance(curr_test["path"], str):
                    getter = GetterPath.fromStr(curr_test["path"])
                else:
                    getter = GetterPath.fromDict(curr_test["path"])
                obs.append(getter.get(curr))
            self.assertEqual(obs, curr_test["expected"])

    def testStr(self):
        test_cases = [
            {  # key => return record[key]
                "getter": GetterPath([KeyToken("id")]),
                "expected": "id"
            },
            {  # key.sub_key => return record[key][sub_key]
                "getter": GetterPath([
                    KeyToken("status"),
                    KeyToken("classifier1")
                ]),
                "expected": "status.classifier1"
            },
            {  # None => return record
                "getter": GetterPath(),
                "expected": ""
            },
            {  # key.sub_key().sub_sub_key => return record[key].sub_key()[sub_sub_key]
                "getter": GetterPath([
                    KeyToken("obj"),
                    MethodToken("getInfo"),
                    KeyToken("B")
                ]),
                "expected": "obj.getInfo().B"
            },
            {  # key.* => return record[key].map(lambda elt: elt)
                "getter": GetterPath([
                    KeyToken("class"),
                    IterableToken(
                        GetterPath([])
                    )
                ]),
                "expected": "class.*"
            },
            {  # key.*.sub_key => return record[key].map(lambda elt: elt[sub_key])
                "getter": GetterPath([
                    KeyToken("class"),
                    IterableToken(
                        GetterPath([
                            KeyToken("status")
                        ])
                    )
                ]),
                "expected": "class.*.status"
            },
            {  # key.1.sub_key => return record[key][1][sub_key]
                "getter": GetterPath([
                    KeyToken("class"),
                    KeyToken("1"),
                    KeyToken("status")
                ]),
                "expected": "class.1.status"
            },
            {  # key.method() => return record[key].method()
                "getter": GetterPath([
                    KeyToken("obj"),
                    MethodToken("getAD")
                ]),
                "expected": "obj.getAD()"
            }
        ]
        for curr_test in test_cases:
            self.assertEqual(
                curr_test["getter"].toStr(), curr_test["expected"]
            )


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    unittest.main()
