#!/usr/bin/env python
import os
import re
from distutils.core import setup


def get_version():
    version = None
    notes_filepath = "RELEASES_NOTES.md"
    if os.path.exists(notes_filepath):
        with open(notes_filepath) as FH:
            first_line = FH.readline()
            version = re.search("^\#\s+.+\s+(.+)\s+\[", first_line).groups()[0]  # Example: "# v2.5.0 [DEV]"
    return version


def load_requirements(path):
    requirements = []
    with open(path) as FH:
        requirements = [elt.strip().replace(" ", "") for elt in FH]
    return requirements


setup(
    name='anacore',
    version=get_version(),
    description='Libraries to read/write/process standard files in NGS.',
    long_description='Anapath Core is package containing libraries to read/write/process standard files in NGS. It has been developped for Anatomo-Cytopathologie department of IUCT Oncopole.',
    author='Frederic Escudie',
    author_email='escudie.frederic@iuct-oncopole.fr',
    license='GNU GPL v3',
    packages=["anacore"],
    install_requires=load_requirements("requirements.txt")
)
