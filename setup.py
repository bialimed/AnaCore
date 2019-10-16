#!/usr/bin/env python
from distutils.core import setup


def load_requirements(path):
    requirements = []
    with open(path) as FH:
        requirements = [elt.strip().replace(" ", "") for elt in FH]
    return requirements


setup(
    name='anacore',
    version='2.5.0',
    description='Libraries to read/write/process standard files in NGS.',
    long_description='Anapath Core is package containing libraries to read/write/process standard files in NGS. It has been developped for Anatomo-Cytopathologie department of IUCT Oncopole.',
    author='Frederic Escudie',
    author_email='escudie.frederic@iuct-oncopole.fr',
    license='GNU GPL v3',
    packages=["anacore"],
    install_requires=load_requirements("requirements.txt")
)
