#!/usr/bin/env python
import os
import re
from setuptools import find_packages, setup


def get_long_description():
    content = ""
    with open("README.md", encoding='utf-8') as FH:
        content = FH.read()
    return content


def get_version():
    version = None
    notes_filepath = "RELEASES_NOTES.md"
    if os.path.exists(notes_filepath):
        with open(notes_filepath) as FH:
            first_line = FH.readline()
            version = re.search(r"^\#\s+.+\s+(.+)\s+\[", first_line).groups()[0]  # Example: "# v2.5.0 [DEV]"
    return version


def load_requirements(path):
    requirements = []
    with open(path) as FH:
        requirements = [elt.strip().replace(" ", "") for elt in FH]
    return requirements


setup(
    name='anacore',
    version=get_version(),
    description='Libraries for managing standard file formats and objects from NGS.',
    author='Frederic Escudie',
    author_email='escudie.frederic@iuct-oncopole.fr',
    license='GNU GPL v3',
    packages=find_packages(),
    install_requires=load_requirements("requirements.txt"),
    url='https://github.com/bialimed/anacore',
    python_requires='>=3.5',
    keywords='bio NGS',
    #long_description_content_type='text/markdown',
    #long_description=get_long_description()
)
