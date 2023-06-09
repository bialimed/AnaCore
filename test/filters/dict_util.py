# -*- coding: utf-8 -*-

def flattenDict(data, result=None, prefix=None):
    if result is None:
        result = {}
    for key, value in data.items():
        long_key = key if prefix is None else "{}.{}".format(prefix, key)
        if type(value) != dict:
            result[long_key] = value
        else:
            flattenDict(value, result, long_key)
    return result
