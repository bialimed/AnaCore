#
# Copyright (C) 2017 IUCT-O
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2017 IUCT-O'
__license__ = 'GNU General Public License'
__version__ = '1.2.2'
__email__ = 'frederic.escudie@iuct-oncopole.fr'
__status__ = 'prod'

from copy import deepcopy



class FiltersCombiner:
    """
    @summary: Wrap as an unique filter a combination of several Filter and
    FiltersCombiner.
    @synopsis:
        ...
        age_filter = Filter(10, "<", "age")
        treatment_filter = filter(["20ng", "20ng + pc"], "in", "treatment")
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
        @param filters: [list] List of Filter and FiltersCombiner evaluated by
        the instance.
        @param operator: [str] The type of evaluation processed between filters.
        With "and" all the filters must have an valid evaluation, with "or" at
        least one filter must have an valid evaluation, with "xor" only one
        filter must have an valid evaluation.
        @param name: [str] The name of the filter.
        @param description: [str] The description of the filter.
        @todo: xor.
        """
        self.description = description
        self.filters = filters
        self.name = name
        self.operator = operator
        if operator not in ["and", "or"]:
            raise Exception(
                'The operator "{}" is not valid for {}.'.format(
                    operator, self.__class__.__name__
                )
            )

    def eval(self, item):
        """
        @summary: Returns True if the item fits the filters regarding to the
        operator.
        @param item: [*] The evaluated element.
        @return: [bool] True if the item fits the filters regarding to the
        operator.
        @todo: xor.
        """
        is_valid = True
        nb_filters = len(self.filters)
        if nb_filters > 0:
            is_valid = None
            idx = 0
            while(idx < nb_filters and is_valid is None):
                if self.filters[idx].eval(item):
                    if self.operator == "or":
                        is_valid = True
                else:
                    if self.operator == "and":
                        is_valid = False
                idx += 1
            if is_valid is None:
                if self.operator == "and":  # The operator is AND and no false has been encountered
                    is_valid = True
                else:  # The operator is OR and no true has been encountered
                    is_valid = False
        return is_valid


class Filter:
    """
    @summary: Manages a filter.
    @synopsis:
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
        children_filter = Filter(13, "<", "age")
        for curr in patients:
            if children_filter.eval(curr):
                print(curr)
        # Select not placebo
        treatment_filter = Filter("placebo", "=", "treatment", action="exclude")
        for curr in patients:
            if treatment_filter.eval(curr):
                print(curr)
        # Select patients with at least group A or group C
        group_filter = Filter(["A", "C"], "in", "group", aggregator="nb:1")
        for curr in patients:
            if group_filter.eval(curr):
                print(curr)
    """
    def __init__(self, values, operator, getter=None, aggregator="ratio:1", name=None, description=None, action="select"):
        """
        @param values: [*] The value(s) used as threshold/reference in item
        evaluation.
        @param operator: [str] The operator used in evaluation.
        @param getter: [str|callable] The property or the function applied on
        evaluated item to get the evaluated value. See Filter.setGetter().
        @param aggregator: [str] The rule used when evaluated value is a list.
        It determines the number of element in list that must fit the filter.
        @param name: [str] The name of the filter.
        @param description: [str] The description of the filter.
        @param action: [str] Determines if the item fitting the filter are kept
        ("select") or rejected ("exclude") by filter.
        """
        self.action = action
        if action not in ["select", "exclude"]:
            raise Exception('The action "{}" is invalid. An action must be "select" or "exclude".'.format(action))
        self.description = description
        self.name = name
        self.values = values
        self.setAggregator(aggregator)
        self.setGetter(getter)
        self.setOperator(operator)

    @staticmethod
    def fromDict(filter_desc):
        """
        @summary: Returns an instance of Filter corresponding to the parameters
        in dictionary.
        @param filter_desc: [dict] The parameters.
        @return: [Filter] An instance of Filter corresponding to the parameters
        in dictionary.
        """
        cleaned_desc = deepcopy(filter_desc)
        if "class" in cleaned_desc:
            cleaned_desc.pop('class', None)
        return Filter(**cleaned_desc)

    def setOperator(self, new):
        """
        @summary: Changes the operator used in evaluation.
        @param new: [str] The new value of the operator.
        """
        self.operator = new
        self.setFct()

    def setAggregator(self, new):
        """
        @summary: Changes the aggregator used to evaluate a list of values.
        @param new: [str] The new value of the aggregator. The aggregator must
        have one of the following form: "nb:INT" or "ratio:FLOAT". The filter
        is ok if the evaluated list returned by getter has:
          - a number of values fitting the Filter superior than INT if
          aggregator is nb:INT.
          - a ratio of values fitting the Filter superior than FLOAT if
          aggregator is ratio:FLOAT.
        """
        self.aggregator = new  # nb:X, ratio:X.X
        agg_type, agg_threshold = new.split(":")
        if agg_type == "nb":
            self._aggregatorEval = lambda nb, length: nb >= int(agg_threshold)
        elif agg_type == "ratio":
            self._aggregatorEval = lambda nb, length: nb >= float(agg_threshold) * length
        else:
            raise Exception('The aggregator "{}" is invalid. An aggregator must start with "nb" or "ratio".'.format(new))

    def setGetter(self, new):
        """
        @summary: Changes the getter used to retrieve evaluated value from evaluated item.
        @param new: [str|callable] The new value of the getter.
        This value can take 3 forms:
            - If the getter is None: the evaluated item himself is used in
            evaluation.
            - If the getter is callable: this is the result of this call applied
            on item that is evaluated. Use a callable is not compatible with
            Filter.toDict().
            - If the getter is a string: the getter is parsed by getRecordValue
            to retrieve a specific property.
        """
        self.getter = new
        self._getterFct = None
        if new is None:
            def returnItem(item):
                return item
            self._getterFct = returnItem
        else:
            if callable(new):
                self._getterFct = new
            else:
                self._getterFct = self.getRecordValue

    def getRecordValue(self, record, _getter=None):
        key = _getter
        if _getter is None:
            key = self.getter
        sub_key = None
        if "." in key:
            key, sub_key = key.split(".", 1)
        # Get value from current key
        value = record
        if key.startswith("m:"):
            method_name = key[2:]
            method_param = None
            if "(" in method_name:  # Parameters must be numeric or string
                if method_name.endswith("()"):
                    method_name = method_name[:-2]
                else:
                    method_name, method_param = method_name.split("(", 1)
                    method_param = [elt.strip() for elt in method_param[:-1].split(",")]
            method = value.__getattribute__(method_name)
            if method_param is None:
                value = method()
            else:
                value = method(*method_param)
        elif key.startswith("i:"):
            if key != "i:":
                value = value.__getattribute__(key[2:])
                if issubclass(value.__class__, dict):
                    new_value = list()
                    for elt_key, elt_value in value.items():
                        new_value.append(elt_value)
                    value = new_value
                elif issubclass(value.__class__, list):
                    new_value = list()
                    for curr_elt in value:
                        new_value.append(curr_elt)
                    value = new_value
                else:
                    raise Exception("The iteration has been ask on non iterable object {}.".format(key[2:]))
        elif issubclass(value.__class__, dict):
            value = value[key]
        elif issubclass(value.__class__, list):
            value = value[int(key)]
        elif issubclass(value.__class__, object):
            value = value.__getattribute__(key)
        # Get value from current sub-key
        if sub_key is not None:
            if key.startswith("i:"):
                new_value = list()
                for curr_elt in value:
                    new_value.append(self.getRecordValue(curr_elt, sub_key))
                value = new_value
            else:
                value = self.getRecordValue(value, sub_key)
        # return
        return value

    def eval(self, item):
        """
        @summary: Returns True if the item fits the filter and action is select
        or if the item does not fit the filter and action is exclude.
        @param item: [*] The evaluated element.
        @return: [bool] True if the item fits the filter and action is select or
        if the item does not fit the filter and action is exclude.
        """
        is_valid = False
        # Eval item
        value = self._getterFct(item)
        if issubclass(value.__class__, list):
            is_valid_list = [self._evalFct(curr_val, self.values) for curr_val in value]
            if self._aggregatorEval(is_valid_list.count(True), len(value)):
                is_valid = True
        else:
            is_valid = self._evalFct(value, self.values)
        # Apply action
        if self.action == "exclude":
            is_valid = not is_valid
        # Return
        return is_valid

    def setFct(self):
        """
        @summary: Sets the evaluation function according to the operator.
        """
        if self.operator in ["=", "==", "eq"]:
            self._evalFct = lambda val, ref: val == ref
        elif self.operator in ["!=", "<>", "ne"]:
            self._evalFct = lambda val, ref: val != ref
        elif self.operator in ["<=", "le"]:
            self._evalFct = lambda val, ref: val <= ref
        elif self.operator in [">=", "ge"]:
            self._evalFct = lambda val, ref: val >= ref
        elif self.operator in ["<", "lt"]:
            self._evalFct = lambda val, ref: val < ref
        elif self.operator in [">", "gt"]:
            self._evalFct = lambda val, ref: val > ref
        elif self.operator == "in":
            self._evalFct = lambda val, ref: val in ref
        elif self.operator == "not in":
            self._evalFct = lambda val, ref: val not in ref
        else:
            raise Exception('The operator "{}" is not implemented.'.format(self.operator))


def filtersFromDict(filters):
    """
    @summary: Returns filter from filters descriptor.
    @param filters: [dict] Filters parameters.
    @return: [Filter|FiltersCombiner] The filter corresponding to descriptor.
    """
    if "class" not in filters or filters["class"] == "Filter":
        return Filter.fromDict(filters)
    elif filters["class"] == "FiltersCombiner":
        cleaned_desc = deepcopy(filters)
        cleaned_desc.pop('class', None)
        cleaned_desc["filters"] = [filtersFromDict(sub_filter) for sub_filter in filters["filters"]]
        return FiltersCombiner(**cleaned_desc)
    else:
        raise Exception('The class "{}" is invalid for filters dictionary. The value must be "Filter" or "FiltersCombiner".'.format(filters["class"]))
