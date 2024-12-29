#!/usr/bin/env python

"""The setup script."""

from setuptools import setup, find_packages

with open('README.rst') as readme_file:
    readme = readme_file.read()

requirements = [
    "anytree>=2.8.0",
    "biopython>=1.79",
    "click>=7.0",
    "dendropy>=5.0.0",
    "Deprecated>=1.2.13",
    "ete3>=3.0.0",
    "iterpop>=0.4.0",
    "nanto>=0.1.1",
    "networkx>=2.5",
    "numpy>=1.21.5",
    "opytional>=0.1.0",
    "pandas>=1.1.0",
    "sortedcontainers>=2.4.0",
    "treeswift",
    "typing-extensions>=4.7.1",
    "validators>=0.20.0",
    "yarl>=1.9.3",
]

setup_requirements = ['pytest-runner', 'tox' ]

test_requirements = ['pytest>=3', 'phylotrackpy>=0.1.16' ]

setup(
    author="Matthew Andres Moreno",
    author_email='m.more500@gmail.com',
    python_requires='>=3.8',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Programming Language :: Python :: 3.12',
        'Programming Language :: Python :: 3.13',
    ],
    description="alifedata-phyloinformatics-convert helps apply traditional phyloinformatics software to alife standardized data",
    entry_points='''
    [console_scripts]
    alifedata-phyloinformatics-convert=alifedata_phyloinformatics_convert.cli:cli
    ''',
    install_requires=requirements,
    extras_require={
        'test': test_requirements,
        'docs': ['sphinx', 'alabaster', 'twine==1.14.0', 'mkdocs==1.2.3', 'jinja2==3.0.0'],
    },
    license="MIT license",
    long_description=readme,
    include_package_data=True,
    keywords='alifedata-phyloinformatics-convert',
    name='alifedata-phyloinformatics-convert',
    packages=find_packages(include=['alifedata_phyloinformatics_convert', 'alifedata_phyloinformatics_convert.*']),
    setup_requires=setup_requirements,
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/mmore500/alifedata-phyloinformatics-convert',
    version='0.18.0',
    zip_safe=False,
)
