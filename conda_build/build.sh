#!/bin/bash
set -ex

# Install AnaCore
cd lib
${PYTHON} -m pip install -vv .

# Install AnaCore-bin
cd -
${PYTHON} -m pip install -vv .
