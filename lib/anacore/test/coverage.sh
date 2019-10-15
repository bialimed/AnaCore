#!/bin/bash
TEST_DIR=`dirname $0`
cd ${TEST_DIR}

if [ -e .coverage ]; then
    coverage erase
fi
coverage run -m unittest discover . > /dev/null 2> /dev/null && \
coverage report
