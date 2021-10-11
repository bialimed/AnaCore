# Create and upload PyPI package

    python setup.py sdist  # Create package
    twine upload dist/*  # Upload to PyPI
    rm dist/*
