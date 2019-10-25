#!/bin/bash
script_path=`realpath $0`
coverage_dir=`dirname ${script_path}`
test_dir=`dirname ${coverage_dir}`
cd ${test_dir}

# Get code coverage
if [ -e .coverage ]; then
    coverage erase
fi
coverage run --source=anacore -m unittest discover . > /dev/null 2> /dev/null && \
coverage report && \
coverage erase
