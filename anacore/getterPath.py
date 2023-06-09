# -*- coding: utf-8 -*-
"""Classes to get specific deep property value from complex objects based on
path (as dict or string).

**GetterPath as string:**

  It is composed of accessors and level separators. Levels are separated by ``.``. Accessors can be:

    * ``attr_name`` for property in object or key in dict or index in list.
    * ``method_name()`` for method in object. If you need arguments in method use GetterPath as dict.
    * ``*`` for iteration in each item of an iterable.

  Examples:
    Case 1:
      .. highlight:: python
      .. code-block:: python

        object = {
            "id": 8,
            "class": {"status": "pathogen", "score": 0.9}
        }
        path = "class.status"
        getter = GetterPath.fromStr(path)
        getter.get(object)  # Return: "pathogen"

    Case 2:
      .. highlight:: python
      .. code-block:: python

        object = {
            "layer": 2,
            "shapes": [<Square (l:2, w:2)>, <Rectange (l:2, w:5)>]
        }
        path = "shapes.*.getArea()"
        getter = GetterPath.fromStr(path)
        getter.get(object)  # Return: [4, 10]

    Case 3:
      .. highlight:: python
      .. code-block:: python

        object = [
            {"class": "stable", "scores": [0.9, 0.05, 0.05]},
            {"class": "stable", "scores": [0.7, 0.2, 0.1]}
        ]
        path = "*.scores.0"
        getter = GetterPath.fromStr(path)
        getter.get(object)  # Return: [0.9, 0.7]


**GetterPath as dict:**

  Is a list of accessors ordered by increasing depth. Accessors can be:

    * `{"kind": "key", "key": attr_name}` for property in object or key in dict or index in list.
    * `{"kind": "method", "name": method_name, "arguments": args_list}` for method in object.
    * `{"kind": "iterable"}` for iteration in each item of an iterable.

  Examples:
    Case 1:
      .. highlight:: python
      .. code-block:: python

        object = {
            "id": 8,
            "class": {"status": "pathogen", "score": 0.9}
        }
        path = [
            {"kind": "key", "key": "class"},
            {"kind": "key", "key": "status"}
        ]
        getter = GetterPath.fromDict(path)
        getter.get(object)  # Return: "pathogen"

    Case 2:
      .. highlight:: python
      .. code-block:: python

        object = {
            "layer": 2,
            "shapes": [<Square (l:2, w:2)>, <Rectange (l:2, w:5)>]
        }
        path = [
            {"kind": "key", "key": "shapes"},
            {"kind": "iterate"},
            {"kind": "method", "name": "scaleUp", "arguments": [2]},
        ]
        getter = GetterPath.fromDict(path)
        getter.get(object)  # Return: [<Square (l:4, w:4)>, <Rectange (l:4, w:10)>]

    Case 3:
      .. highlight:: python
      .. code-block:: python

        object = [
            {"class": "stable", "scores": [0.9, 0.05, 0.05]},
            {"class": "stable", "scores": [0.7, 0.2, 0.1]}
        ]
        path = [
            {"kind": "iterate"},
            {"kind": "key", "key": "scores"},
            {"kind": "key", "key": 0},
        ]
        getter = GetterPath.fromDict(path)
        getter.get(object)  # Return: [0.9, 0.7]
"""

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2023 CHU Toulouse'
__license__ = 'GNU General Public License'
__version__ = '1.0.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'


class GetterPath:
    """Class to manage an accessor on deep property value from complex objects based on a path (as dict or string)."""

    def __init__(self, tokens_chain=None):
        """
        Build and return an instance of GetterPath.

        :param tokens_chain: Accessors tokens for sub-levels.
        :type tokens_chain: list
        :return: The new instance.
        :rtype: getterPath.GetterPath
        """
        self.tokens_chain = tokens_chain if tokens_chain else []

    def __eq__(self, other):
        """
        Return True if self and other have same path.

        :param other: GetterPath compared to self.
        :type other: GetterPath
        :return: True if self and other have same path.
        :rtype: bool
        """
        return self.toDict() == other.toDict()

    @staticmethod
    def fromDict(array):
        """
        Return GetterPath instance corresponding to path.

        :param array: List to define accessor on deep property value from complex objects.
        :type array: list
        :return: GetterPath instance corresponding to path.
        :rtype: GetterPath
        """
        if array is None or len(array) == 0:
            return GetterPath()
        # Create tokens chain
        getter = GetterPath()
        tokens_chain = getter.tokens_chain
        for token_dict in array:
            if token_dict["kind"] == "iterable":
                sub_getter = GetterPath()  # Management of list and list iterate in path is different
                tokens_chain.append(IterableToken(sub_getter))
                tokens_chain = sub_getter.tokens_chain
            elif token_dict["kind"] == "method":
                args = token_dict["arguments"] if "arguments" in token_dict else list()
                tokens_chain.append(
                    MethodToken(token_dict["name"], args)
                )
            elif token_dict["kind"] == "key":
                tokens_chain.append(KeyToken(token_dict["key"]))
            else:
                raise ValueError("Kind is not one of ['iterable', 'key', 'method'] in token {}.".format(token_dict))
        return getter

    @staticmethod
    def fromStr(path):
        """
        Return GetterPath instance corresponding to path.

        :param path: String to define accessor on deep property value from complex objects.
        :type path: str
        :return: GetterPath instance corresponding to path.
        :rtype: GetterPath
        """
        if path is None or path == "":
            return GetterPath()
        # Create tokens chain
        getter = GetterPath()
        tokens_chain = getter.tokens_chain
        for token_str in path.split("."):
            if token_str == "*":  # *.x, x.*, x.*.y
                sub_getter = GetterPath()  # Management of list and list iterate in path is different
                tokens_chain.append(IterableToken(sub_getter))
                tokens_chain = sub_getter.tokens_chain
            elif token_str.endswith("()"):
                tokens_chain.append(MethodToken(token_str[:-2]))
            else:
                tokens_chain.append(KeyToken(token_str))
        return getter

    def get(self, item):
        """
        Return value corresponding to getter from item.

        :param item: Object on wich getter is apply.
        :type item: *
        :return: Value corresponding to getter from item or None if element is undefined.
        :rtype: * | None
        """
        val = item
        for curr_token in self.tokens_chain:
            val = curr_token.get(val)
            if val is None:
                return None
        return val

    def toDict(self):
        """
        Return getter path as dict path.

        :return: Getter path as dict path.
        :rtype: str
        """
        hash = list()
        for curr_token in self.tokens_chain:
            if isinstance(curr_token, IterableToken):
                for sub_token_dict in curr_token.toDict():
                    hash.append(sub_token_dict)
            else:
                hash.append(curr_token.toDict())
        return hash

    def toStr(self):
        """
        Return getter path as str path.

        :return: Getter path as str path.
        :rtype: str
        """
        return ".".join([curr_token.toStr() for curr_token in self.tokens_chain])


class IterableToken:
    """Getter to manage sub-level getters on each item of an iterable."""

    def __init__(self, sub_getter=None):
        """
        Build and return an instance of IterableToken.

        :param sub_getter: Sub-levels accessors from iteation point to deeper accessor.
        :type sub_getter: GetterPath
        :return: The new instance.
        :rtype: getterPath.IterableToken
        """
        self.sub_getter = sub_getter if sub_getter else GetterPath()

    def __eq__(self, other):
        """
        Return True if self and other have same path.

        :param other: IterableToken compared to self.
        :type other: IterableToken
        :return: True if self and other have same path.
        :rtype: bool
        """
        return self.toDict() == other.toDict()

    def get(self, item):
        """
        Return value corresponding to getter from item.

        :param item: Object on wich getter is apply.
        :type item: *
        :return: Values corresponding to getter from item.
        :rtype: list
        """
        if item is None:
            return list()
        if issubclass(item.__class__, dict):
            item = [elt_value for elt_key, elt_value in item.items()]
        elif issubclass(item.__class__, list):
            pass
        else:
            raise Exception("Iteration cannot be used on {}.".format(item))
        values = list()
        for elt in item:
            sub_value = self.sub_getter.get(elt)
            if isinstance(sub_value, list):
                # Flatten list and sub list
                for sub_elt in sub_value:
                    values.append(sub_elt)
            else:
                values.append(sub_value)
        return values

    def toDict(self):
        """
        Return getter token and sub-token as dict path.

        :return: Getter token and sub-token as dict path.
        :rtype: str
        """
        hash = [{"kind": "iterable"}]
        for curr_token_dict in self.sub_getter.toDict():
            hash.append(curr_token_dict)
        return hash

    def toStr(self):
        """
        Return getter token and sub-token as str path.

        :return: Getter token and sub-token as str path.
        :rtype: str
        """
        if len(self.sub_getter.tokens_chain) > 0:
            return "*.{}".format(self.sub_getter.toStr())
        return "*"


class KeyToken:
    """Getter to manage get by property in object or key in dict or index in list."""

    def __init__(self, key):
        """
        Build and return an instance of KeyToken.

        :param key: Property name or key name or index.
        :type key: str
        :return: The new instance.
        :rtype: getterPath.KeyToken
        """
        self.key = key

    def __eq__(self, other):
        """
        Return True if self and other have same path.

        :param other: KeyToken compared to self.
        :type other: KeyToken
        :return: True if self and other have same path.
        :rtype: bool
        """
        return self.toDict() == other.toDict()

    def get(self, item):
        """
        Return value corresponding to getter from item.

        :param item: Object on wich getter is apply.
        :type item: *
        :return: Value corresponding to getter from item or None if element is undefined.
        :rtype: * | None
        """
        if item is None:
            return None
        if issubclass(item.__class__, dict):
            return item.get(self.key, None)
        if issubclass(item.__class__, list):
            return item[int(self.key)]
        if issubclass(item.__class__, object):
            return getattr(item, self.key, None)
        raise Exception("The key {} cannot be used on item {}.".format(self.key, item))

    def toDict(self):
        """
        Return getter token as dict.

        :return: Getter token as dict.
        :rtype: str
        """
        return {"kind": "key", "key": self.key}

    def toStr(self):
        """
        Return getter token as str.

        :return: Getter token as str.
        :rtype: str
        """
        return str(self.key)


class MethodToken:
    """Getter to manage get by method in object."""

    def __init__(self, name, args=None):
        """
        Build and return an instance of MethodToken.

        :param name: method name.
        :type name: str
        :param args: List of values for method arguments.
        :type args: list
        :return: The new instance.
        :rtype: getterPath.MethodToken
        """
        self.name = name
        self.args = args if args else list()

    def __eq__(self, other):
        """
        Return True if self and other have same path.

        :param other: MethodToken compared to self.
        :type other: MethodToken
        :return: True if self and other have same path.
        :rtype: bool
        """
        return self.toDict() == other.toDict()

    def get(self, item):
        """
        Return value corresponding to getter from item.

        :param item: Object on wich getter is apply.
        :type item: *
        :return: Value corresponding to getter from item or None if element is undefined.
        :rtype: * | None
        """
        method = getattr(item, self.name, None)
        if method is None:
            return None
        if len(self.args) == 0:
            return method()
        return method(*self.args)

    def toDict(self):
        """
        Return getter token as dict.

        :return: Getter token as dict.
        :rtype: str
        """
        hash = {"kind": "method", "name": self.name}
        if self.args:
            hash["arguments"] = self.args
        return hash

    def toStr(self):
        """
        Return getter token as str.

        :return: Getter token as str.
        :rtype: str
        """
        if len(self.args) != 0:
            raise Exception("Method with parameters is not implemented for GetterPath as str. Please dict representation instead.")
        return "{}()".format(self.name)
