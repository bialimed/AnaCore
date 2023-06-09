# -*- coding: utf-8 -*-
"""Classes and functions for managing and processing complex filters like included
filters and combined filters with differents operators.

Filters can be:
  * ``EmptyIterfilter``: Filter class to select or exclude empty iterable value
    (see example above: "Filter patients with empty list of group").
  * ``Filter``: Filter class to select or exclude based on comparison operator
    and reference/threshold value (see example above: "Filter on age"). Operator
    attribute defines the comparison to process bewteen a reference/threshold
    value and the item value (see Filters.setFct).
  * ``FiltersCombiner``: Filter class to manage a combination of filters
    (EmptyIterfilter, Filter and FiltersCombiner). These filters can be combined
    by "and", "or" or "xor" combination (see example above: ilter patients with
    (group with at least C or age lower than 20) and age != None).

Get value from object when eval filter:
  Evaluated value from evaluated object is accessed with a GetterPath as string
  (see `anacore.getterPath.GetterPath`). Briefely, GetterPath describes the
  attributes/key/method to apply on object to get value. This GetterPath is
  provide in ``getter`` attribute in ``EmptyIterfilter`` and ``Filter``.
  Otherwise all the object is evaluated.
  
  Example:
    .. highlight:: python
    .. code-block:: python

      patient = {
          "id": 1,
          "classification": {
              "clinvar": "pathogenic",
              "cosmic": "unknown_significance"
          }
      }
      
      filter = Filter("=", "pathogenic", "classification.clinvar")
      filter.eval(patient)  # Return True

Examples:
  .. highlight:: python
  .. code-block:: python

    patients = [
        {"id": 1, "age": 12, "treatment": "drug A", "group": ["A", "B"]},
        {"id": 2, "age": 32, "treatment": None, "group": ["C"]},
        {"id": 3, "age": 64, "treatment": "drug B", "group": []},
        {"id": 4, "age": None, "treatment": "placebo", "group": None}
    ]

  :Filter on simple value:

  * Simple filter on age:

    .. highlight:: python
    .. code-block:: python

      filter = Filter("<", 18, "age")
      [elt["id"] for elt in patients if filter.eval(elt)]  # Return: [1]

  * Simple filter select patients where treatment contains drug:

    .. highlight:: python
    .. code-block:: python

      filter = Filter("contains", "drug", "treatment")
      [elt["id"] for elt in patients if filter.eval(elt)]  # Return: [1, 3]

  * Simple filter exclude patients where treatment contains drug:

    .. highlight:: python
    .. code-block:: python

      filter = Filter("contains", "drug", "treatment", action="exclude")
      [elt["id"] for elt in patients if filter.eval(elt)]  # Return: [2, 4]

  :Filter on list:

  * Filter patients with empty list of group:

    .. highlight:: python
    .. code-block:: python

      filter = EmptyIterFilter("group")
      [elt["id"] for elt in patients if filter.eval(elt)]  # Return: [3, 4]


  * Filter patients with at least one group equal to B:

    .. highlight:: python
    .. code-block:: python

      filter = Filter("=", "B", "group", aggregator="nb:1")
      [elt["id"] for elt in patients if filter.eval(elt)]  # Return: [1]

  * Filter patients with at least one group different to B:

    .. highlight:: python
    .. code-block:: python

      filter = Filter("<>", "B", "group", aggregator="nb:1")
      [elt["id"] for elt in patients if filter.eval(elt)]  # Return: [1, 2]

  * Filter patients with at all groups different to B:

    .. highlight:: python
    .. code-block:: python

      filter = Filter("<>", "B", "group", aggregator="ratio:1")
      [elt["id"] for elt in patients if filter.eval(elt)]  # Return: [2, 3, 4]

  :Combine filters:

  * Filter patients with (group with at least C or age lower than 20) and age != None:

    .. highlight:: python
    .. code-block:: python

      filter = FiltersCombiner(
        [
            FiltersCombiner(
                [
                    Filter("==", "C", "group", aggregator="nb:1"),
                    Filter("<", 20, "age")
                ],
                "or"
            ),
            Filter("<>", None, "age")
        ],
        "and"
      )
      [elt["id"] for elt in patients if filter.eval(elt)]  # Return: [1, 2]

Define filter in JSON file:
  Filters can be defined in JSON file and load from it.

  Example of JSON:
    .. highlight:: json
    .. code-block:: json

      {
          "class": "FiltersCombiner",
          "filters": [
              {
                  "class": "FiltersCombiner",
                  "filters": [
                      {"aggregator": "nb:1", "class": "Filter", "getter": "group", "operator": "=", "values": "C"},
                      {"class": "Filter", "getter": "age", "operator": "<", "values": 20}
                  ],
                  "operator": "and"
              },
              {
                  "class": "Filter",
                  "getter": "age",
                  "operator": "<>",
                  "values": null
              },
          ],
          "operator": "or"
      }

  Load filters:
    .. highlight:: python
    .. code-block:: python

      import json
      from anacore.filters import filtersFromDict
    
      with open("filters.json") as reader:
          rules = json.load(reader)
          filter = filtersFromDict(rules)
"""

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2017 CHU Toulouse'
__license__ = 'GNU General Public License'
__version__ = '2.0.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

from anacore.getterPath import GetterPath
from copy import deepcopy


class AbstractAtomicFilter:
    def __init__(self, getter=None, name=None, description=None, action="select"):
        """
        Build and return an instance of AbstractAtomicFilter.

        :param getter: The property or the function applied on evaluated item to get the evaluated value. See Filter.setGetter().
        :type getter: str or callable
        :param name: The name of the filter.
        :type name: str
        :param description: The description of the filter.
        :type description: str
        :param action: Determines if the item fitting the filter are kept ("select") or rejected ("exclude") by filter.
        :type action: str
        :return: The new instance.
        :rtype: filters.AbstractAtomicFilter
        """
        self.action = action
        if action not in ["select", "exclude"]:
            raise Exception('The action "{}" is invalid. An action must be "select" or "exclude".'.format(action))
        self.description = description
        self.name = name
        self.setGetter(getter)

    def __eq__(self, other):
        """
        Return True if self and other have same filtering behaviour.

        :param other: Filter compared to self.
        :type other: AbstractAtomicFilter
        :return: True if self and other have same filtering behaviour.
        :rtype: bool
        """
        self_dict = {key: val for key, val in self.__dict__.items() if not (key.startswith("_") and key.endswith("Fct"))}
        other_dict = {key: val for key, val in other.__dict__.items() if not (key.startswith("_") and key.endswith("Fct"))}
        return self_dict == other_dict

    def eval(self, item):
        """
        Return True if the item fits the filter and action is select or if the item does not fit the filter and action is exclude.

        :param item: The evaluated element.
        :type item: *
        :return: True if the item fits the filter and action is select or if the item does not fit the filter and action is exclude.
        :rtype: bool
        """
        raise NotImplementedError

    def setGetter(self, new):
        """
        Change the getter used to retrieve evaluated value from evaluated item.

        :param new: [str|callable] The new value of the getter. This value can take 3 forms:
          (1) If the getter is None: the evaluated item himself is used in evaluation.
          (2) If the getter is callable: this is the result of this call applied on item that is evaluated. Use a callable is not compatible with Filter.toDict().
          (3) If the getter is a string: An instance of  GetterPath is created from getter an it is used to get a specific property.
        :type new: str
        """
        self.getter = new
        self._getterFct = None
        if callable(new):
            self._getterFct = new
        else:
            self._getterFct = GetterPath.fromStr(new).get

    def toDict(self):
        """
        Return the dictionary corresponding to the instance of Filter.

        :return: The dictionary corresponding to the instance of Filter.
        :rtype: dict
        """
        obj_dict = {"class": self.__class__.__name__}
        for attr_name in dir(self):
            if not attr_name.startswith("_"):
                attr_value = self.__getattribute__(attr_name)
                if not callable(attr_value):
                    obj_dict[attr_name] = attr_value
        return obj_dict


class EmptyIterFilter(AbstractAtomicFilter):
    """
    Special filter on empty iterable.

    :Code example:
        .. highlight:: python
        .. code-block:: python

            patients = [
                {"age": 8, "treatment": "placebo", "group": ["A", "B"]},
                {"age": 9, "treatment": "20ng ip", "group": ["C"]},
                {"age": 11, "treatment": "placebo", "group": []},
                {"age": 18, "treatment": "placebo", "group": None},
                {"age": 25, "treatment": "20ng ip"}
            ]
            # Select patient in no group
            wout_gp_filter = EmptyIterFilter("group")
            for curr in patients:
                if wout_gp_filter.eval(curr):
                    print(curr)

            # Result>
            # {"age": 11, "treatment": "placebo", "group":[]}
            # {"age": 18, "treatment": "placebo", "group": None}
            # {"age": 25, "treatment": "20ng ip"}
    """

    def eval(self, item):
        """
        Return True if the item fits the filter and action is select or if the item does not fit the filter and action is exclude.

        :param item: The evaluated element.
        :type item: *
        :return: True if the item fits the filter and action is select or if the item does not fit the filter and action is exclude.
        :rtype: bool
        """
        is_valid = False
        # Get value
        value = self._getterFct(item)
        # Eval item
        if value is None:  # Missing item is equivalent to empty
            is_valid = True
        else:
            try:
                iter(value)
            except Exception:
                raise ValueError("Value inspected by filter is not iterable in item {}.".format(item))
            if isinstance(value, str):
                raise ValueError("Value inspected by filter is string in item {}.".format(item))
            is_valid = len(value) == 0
        # Apply action
        if self.action == "exclude":
            is_valid = not is_valid
        # Return
        return is_valid

    @staticmethod
    def fromDict(filter_desc):
        """
        Return an instance of EmptyIterFilter corresponding to the parameters in dictionary.

        :param filter_desc: The parameters.
        :type: dict
        :return: An instance of Filter corresponding to the parameters in dictionary.
        :rtype: filters.EmptyIterFilter
        """
        cleaned_desc = deepcopy(filter_desc)
        if "class" in cleaned_desc:
            cleaned_desc.pop('class', None)
        return EmptyIterFilter(**cleaned_desc)


class Filter(AbstractAtomicFilter):
    """
    Manage one filter.

    :Code example:
        .. highlight:: python
        .. code-block:: python

            patients = [
                {"age": 8, "treatment": "placebo", "group":["A", "B"]},
                {"age": 9, "treatment": "20ng ip", "group":["C"]},
                {"age": 11, "treatment": "placebo", "group":["A", "D"]},
                {"age": 25, "treatment": "20ng ip", "group":["B"]},
                {"age": 30, "treatment": "20ng pc", "group":["C", "B"]},
                {"age": 33, "treatment": "20ng pc", "group":["A", "D"]},
                {"age": 30, "treatment": "placebo", "group":["A", "B"]},
            ]
            # Select children
            children_filter = Filter("<", 13, "age")
            for curr in patients:
                if children_filter.eval(curr):
                    print(curr)
            # Select not placebo
            treatment_filter = Filter("=", "placebo", "treatment", action="exclude")
            for curr in patients:
                if treatment_filter.eval(curr):
                    print(curr)
            # Select patients with at least group A or group C
            group_filter = Filter("in", ["A", "C"], "group", aggregator="nb:1")
            for curr in patients:
                if group_filter.eval(curr):
                    print(curr)
    """

    def __init__(self, operator, values, getter=None, aggregator=None, name=None, description=None, action="select"):
        """
        Build and return an instance of FiltersCombiner.

        :param values: The value(s) used as threshold/reference in item evaluation.
        :type values: *
        :param operator: The operator used in evaluation. For authorized values see doctring of Filters.setFct.
        :type operator: str
        :param getter: The property or the function applied on evaluated item to get the evaluated value. See Filter.setGetter().
        :type getter: str or callable
        :param aggregator: The rule used when evaluated value is a list. It determines the number of element in list that must fit the filter.
        :type aggregator: str
        :param name: The name of the filter.
        :type name: str
        :param description: The description of the filter.
        :type description: str
        :param action: Determines if the item fitting the filter are kept ("select") or rejected ("exclude") by filter.
        :type action: str
        :return: The new instance.
        :rtype: filters.FiltersCombiner
        """
        super().__init__(getter, name, description, action)
        self.name = name
        self._operator = None
        self.values = values
        self.aggregator = aggregator
        self.operator = operator

    @property
    def aggregator(self):
        """
        Return rule used when evaluated value is a list. It determines the number of element in list that must fit the filter.

        :return: Return rule used when evaluated value is a list.
        :rtype: str
        """
        return self._aggregator

    @aggregator.setter
    def aggregator(self, new):
        """
        Change the aggregator used to evaluate a list of values.

        :param new: [str] The new value of the aggregator. The aggregator must have one of the following form: "nb:INT" or "ratio:FLOAT". The filter is ok if the evaluated list returned by getter has:
          (1) a number of values fitting the Filter superior than INT if aggregator is nb:INT.
          (2) a ratio of values fitting the Filter superior than FLOAT if aggregator is ratio:FLOAT.
        :type new: str
        """
        self._aggregator = new  # nb:X, ratio:X.X
        self._aggregatorEvalFct = None
        if new is not None:
            agg_type, agg_threshold = new.split(":")
            if agg_type == "nb":
                if not agg_threshold.isdigit():
                    raise ValueError("The value {} in aggregator {} is not a positive integer.".format(agg_threshold, new))
                int_agg_threshold = int(agg_threshold)
                self._aggregatorEvalFct = lambda nb, length: nb >= int_agg_threshold
            elif agg_type == "ratio":
                float_agg_threshold = None
                try:
                    float_agg_threshold = float(agg_threshold)
                except Exception:
                    raise ValueError("The value {} in aggregator {} is not a float.".format(agg_threshold, new))
                if float_agg_threshold < 0 or float_agg_threshold > 1:
                    raise ValueError("The value {} in aggregator {} is not between 0.0 and 1.0.".format(agg_threshold, new))
                self._aggregatorEvalFct = lambda nb, length: nb >= float_agg_threshold * length
            else:
                raise AttributeError('The aggregator "{}" is invalid. An aggregator must start with "nb" or "ratio".'.format(new))

    def eval(self, item):
        """
        Return True if the item fits the filter and action is select or if the item does not fit the filter and action is exclude.

        :param item: The evaluated element.
        :type item: *
        :return: True if the item fits the filter and action is select or if the item does not fit the filter and action is exclude.
        :rtype: bool
        """
        is_valid = False
        # Eval item
        value = self._getterFct(item)
        if self.aggregator is not None:
            if value is None:  # None is like empty
                value = list()
            is_valid_list = [self._evalFct(curr_val, self.values) for curr_val in value]
            if self._aggregatorEvalFct(is_valid_list.count(True), len(value)):
                is_valid = True
        else:
            is_valid = self._evalFct(value, self.values)
        # Apply action
        if self.action == "exclude":
            is_valid = not is_valid
        # Return
        return is_valid

    @staticmethod
    def fromDict(filter_desc):
        """
        Return an instance of Filter corresponding to the parameters in dictionary.

        :param filter_desc: The parameters.
        :type: dict
        :return: An instance of Filter corresponding to the parameters in dictionary.
        :rtype: filters.Filter
        """
        cleaned_desc = deepcopy(filter_desc)
        if "class" in cleaned_desc:
            cleaned_desc.pop('class', None)
        return Filter(**cleaned_desc)

    @property
    def operator(self):
        """
        Return operator used in evaluation.

        :return: Operator used in evaluation.
        :type: str
        """
        return self._operator

    @operator.setter
    def operator(self, new):
        """
        Change the operator used in evaluation.

        :param new: The new value of the operator.
        :type new: str
        """
        self._operator = new
        self.setFct()

    def setFct(self):
        """
        Set the evaluation function according to the operator.

        Authorized operators:
            * ``=``, ``==``, ``eq``: Equality on reference/threshold and value return by evaluated item.
            * ``!=``, ``<>``, ``ne``: Inequality on reference/threshold and value return by evaluated item.
            * ``<=``, ``le``: Value is lower or equal than reference/threshold. Return False if value is None.
            * ``>=``, ``ge``: Value is greater or equal than reference/threshold. Return False if value is None.
            * ``<``, ``lt``:  Value is lower than reference/threshold. Return False if value is None.
            * ``>``, ``gt``:  Value is greater than reference/threshold. Return False if value is None.
            * ``in``: Value is in list of authorized values.
            * ``not in``: Value is not in list of excluded values.
            * ``contains``: Value to string contains the reference.
            * ``substring of``: Value to string is a substring of the reference.
            * ``empty``: Value is None or is an empty string.
        """
        if self.operator in ["=", "==", "eq"]:
            if type(self.values) == bool:
                self._evalFct = lambda val, ref: val == ref and type(val) == bool
            else:
                self._evalFct = lambda val, ref: val == ref and type(val) != bool
        elif self.operator in ["!=", "<>", "ne"]:
            if type(self.values) == bool:
                self._evalFct = lambda val, ref: val != ref or type(val) != bool
            else:
                self._evalFct = lambda val, ref: val != ref or type(val) == bool
        elif self.operator in ["<=", "le"]:
            if self.values is None or isinstance(self.values, bool):
                raise AttributeError("Reference value in filter cannot be '{}' in '{}' comparison.".format(self.values, self.operator))
            self._evalFct = lambda val, ref: False if val is None else val <= ref
        elif self.operator in [">=", "ge"]:
            if self.values is None or isinstance(self.values, bool):
                raise AttributeError("Reference value in filter cannot be '{}' in '{}' comparison.".format(self.values, self.operator))
            self._evalFct = lambda val, ref: False if val is None else val >= ref
        elif self.operator in ["<", "lt"]:
            if self.values is None or isinstance(self.values, bool):
                raise AttributeError("Reference value in filter cannot be '{}' in '{}' comparison.".format(self.values, self.operator))
            self._evalFct = lambda val, ref: False if val is None else val < ref
        elif self.operator in [">", "gt"]:
            if self.values is None or isinstance(self.values, bool):
                raise AttributeError("Reference value in filter cannot be '{}' in '{}' comparison.".format(self.values, self.operator))
            self._evalFct = lambda val, ref: False if val is None else val > ref
        elif self.operator == "in":  # Compare to list of accepted values. Filter group in ["A", "B"] with data={"group": "A"}
            if self.values is None or isinstance(self.values, str) or not (hasattr(self.values, "__iter__") or hasattr(self.values, "__contains__")):
                raise AttributeError("Reference value in filter cannot be '{}' in '{}' comparison.".format(self.values, self.operator))
            if len(self.values) == 0:
                raise AttributeError("Reference value in filter cannot be empty in '{}' comparison.".format(self.operator))
            self._evalFct = lambda val, ref: val in ref
        elif self.operator == "contains":  # Value as string contains reference
            if not issubclass(self.values.__class__, str):
                raise AttributeError('Reference value in filter must be a string for "contains" operator. If you want to apply "contains" with several possibility you must declare n filters and combine them with FiltersCombiner.')
            self._evalFct = lambda val, ref: False if val is None else ref in str(val)
        elif self.operator == "substring of":  # Value as string is substring of reference
            if not issubclass(self.values.__class__, str):
                raise AttributeError('Reference value in filter must be a string for "substring of" operator.')
            self._evalFct = lambda val, ref: False if val is None else str(val) in ref
        elif self.operator == "not in":  # Compare to list of rejected values. Filter group not in ["A", "B"] with data={"group": "C"}
            if self.values is None or isinstance(self.values, str) or not (hasattr(self.values, "__iter__") or hasattr(self.values, "__contains__")):
                raise AttributeError("Reference value in filter cannot be '{}' in '{}' comparison.".format(self.values, self.operator))
            if len(self.values) == 0:
                raise AttributeError("Reference value in filter cannot be empty in '{}' comparison.".format(self.operator))
            self._evalFct = lambda val, ref: val not in ref
        elif self.operator == "empty":  # Values is missing or is an empty string
            if self.values is not None:
                raise AttributeError("Reference value in filter with operator 'empty' must be None.")
            self._evalFct = lambda val, ref: val is None or (issubclass(val.__class__, str) and val == "")
        else:
            raise AttributeError('The operator "{}" is not implemented.'.format(self.operator))

    @property
    def values(self):
        """
        Return value(s) used as threshold/reference in item evaluation.

        :return: Value(s) used as threshold/reference in item evaluation.
        :rtype: *
        """
        return self._values

    @values.setter
    def values(self, new):
        """
        Change reference values used in evaluation.

        :param new: The new value for values.
        :type new: *
        """
        self._values = new
        if self._operator:
            self.setFct()


class FiltersCombiner:
    """
    Wrap as an unique filter a combination of several Filter and FiltersCombiner.

    :Code example:
        .. highlight:: python
        .. code-block:: python

            ...
            age_filter = Filter("<", 10, "age")
            treatment_filter = filter("in", ["20ng", "20ng + pc"], "treatment")
            filters = FiltersCombiner(
                [age_filter, treatment_filter],
                "and",
                "young_with_treatment"
            )
            for curr in patients:
                if filters.eval(curr):
                    print(curr)
    """

    def __init__(self, filters, operator="and", name=None, description=None):
        """
        Build and return an instance of FiltersCombiner.

        :param filters: List of Filter and FiltersCombiner evaluated by the instance.
        :type filters: list
        :param operator: The type of evaluation processed between filters. With "and" all the filters must have an valid evaluation, with "or" at least one filter must have an valid evaluation, with "xor" only one filter must have an valid evaluation.
        :type operator: str
        :param name: The name of the filter.
        :type name: str
        :param description: The description of the filter.
        :type description: str
        :return: The new instance.
        :rtype: filters.FiltersCombiner
        """
        self.description = description
        self.filters = filters
        self.name = name
        self.operator = operator
        if operator not in {"and", "or", "xor"}:
            raise Exception(
                'The operator "{}" is not valid for {}.'.format(
                    operator, self.__class__.__name__
                )
            )

    def __eq__(self, other):
        """
        Return True if self and other have same filtering behaviour and name and description.

        :param other: Filter compared to self.
        :type other: ComplexFilters
        :return: True if self and other have same filtering behaviour and name and description.
        :rtype: bool
        """
        self_dict = {
            "filters": [elt for elt in self.filters if not issubclass(elt.__class__, FiltersCombiner) or len(elt.filters) > 0],
            "operator": self.operator,
            "name": self.name,
            "description": self.description
        }
        other_dict = {
            "filters": [elt for elt in other.filters if not issubclass(elt.__class__, FiltersCombiner) or len(elt.filters) > 0],
            "operator": other.operator,
            "name": other.name,
            "description": other.description
        }
        return self_dict == other_dict

    def eval(self, item):
        """
        Return True if the item fits the filters regarding to the operator.

        :param item: The evaluated element.
        :type: *
        :return: True if the item fits the filters regarding to the operator.
        :rtype: bool
        :todo: xor.
        """
        is_valid = True
        nb_filters = len(self.filters)
        if nb_filters > 0:
            is_valid = None
            nb_true = 0
            idx = 0
            while(idx < nb_filters and is_valid is None):
                if self.filters[idx].eval(item):  # Is valid
                    if self.operator == "or":
                        is_valid = True
                    nb_true += 1
                else:  # Is invalid
                    if self.operator == "and":
                        is_valid = False
                idx += 1
            if self.operator == "xor":
                if nb_true == 1:
                    is_valid = True
                else:
                    is_valid = False
            if is_valid is None:  # Without value to test
                if self.operator == "and":  # The operator is AND and no false has been encountered
                    is_valid = True
                else:  # The operator is OR or XOR and no true has been encountered
                    is_valid = False
        return is_valid


def filtersFromDict(filters):
    """
    Return filter from filters descriptor.

    :param filters: Filters parameters.
    :type filters: dict
    :return: The filter corresponding to descriptor.
    :rtype: filters.EmptyIterFilter or filters.Filter or filters.FiltersCombiner
    """
    if "class" not in filters or filters["class"] == "Filter":
        return Filter.fromDict(filters)
    elif filters["class"] == "FiltersCombiner":
        cleaned_desc = deepcopy(filters)
        cleaned_desc.pop('class', None)
        cleaned_desc["filters"] = [filtersFromDict(sub_filter) for sub_filter in filters["filters"]]
        return FiltersCombiner(**cleaned_desc)
    elif filters["class"] == "EmptyIterFilter":
        return EmptyIterFilter.fromDict(filters)
    else:
        raise Exception('The class "{}" is invalid for filters dictionary. The value must be "EmptyIterFiler", "Filter" or "FiltersCombiner".'.format(filters["class"]))
