# Package on pypi

## Build dependencies

    python3 -m venv venv_build
    source venv_build/bin/activate
    python3 -m pip install --upgrade build
    deactivate

## Create package

    source venv_build/bin/activate
    python -m build
    deactivate

## Test

    python3 -m venv venv_test
    source venv_test/bin/activate
    python3 -m pip install dist/anacore-*.tar.gz
    ##### Test
    deactivate
    rm -r venv_test

### Push package

    twine upload dist/*
    rm dist/*
