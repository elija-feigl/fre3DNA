from fre3Dna.version import get_version
from setuptools import setup, find_packages

"""Lattice free 3D DNA Origami."""

with open("README.md", "r") as fh:
    long_description = fh.read()

with open("LICENSE", "r") as fh:
    license = fh.read()

setup(
    name="fre3Dna",
    version=get_version(),
    author="Elija Feigl",
    author_email="elija.feigl@tum.de",
    description=__doc__,
    license=license,
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/elija-feigl/fre3Dna",
    packages=find_packages(),
    include_package_data=True,
    install_requires=(
        'click',
    ),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: GNU General Public License Version 3",
        "Operating System :: OS Independent",
    ],
    entry_points='''
        [console_scripts]
        fre3Dna=fre3Dna.scripts.fre3Dna:cli
    ''',
)
