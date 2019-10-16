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


def load_scripts(path):
    scripts = []
    for filename in os.listdir(path):
        filepath = os.path.join(path, filename)
        if os.path.isdir(filepath):
            if filename != "test":
                load_scripts(filepath)
        else:
            if filename.endswith(".py") and not filename.startswith("__"):
                scripts.append(filepath)
    return scripts


def load_requirements(path):
    requirements = []
    with open(path) as FH:
        requirements = [elt.strip().replace(" ", "") for elt in FH]
    return requirements


setup(
    name='anacore-bin',
    version=get_version(),
    description='Scripts to process NGS analyses.',
    author='Frederic Escudie',
    author_email='escudie.frederic@iuct-oncopole.fr',
    license='GNU GPL v3',
    packages=["bin"],
    install_requires=load_requirements("requirements.txt"),
    scripts=load_scripts("bin")
)
