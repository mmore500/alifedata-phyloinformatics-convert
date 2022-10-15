#!/usr/bin/env python

"""The setup script."""

from setuptools import setup, find_packages

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = [
    "biopython>=1.79",
    "click>=7.0",
    "dendropy>=4.5.2",
    "iterpop>=0.4.0",
    "lyncs-utils>=0.3.2",
    "nanto>=0.1.1",
    "numpy>=1.21.5",
    "opytional>=0.1.0",
    "pandas>=1.1.0",
    "sortedcontainers>=2.4.0",
]

setup_requirements = ['pytest-runner', ]

test_requirements = ['pytest>=3', ]

setup(
    author="Matthew Andres Moreno",
    author_email='m.more500@gmail.com',
    python_requires='>=3.6',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
    ],
    description="alifedata-phyloinformatics-convert helps apply traditional phyloinformatics software to alife standardized data",
    entry_points='''
    [console_scripts]
    alifedata-phyloinformatics-convert=alifedata_phyloinformatics_convert.cli:cli
    ''',
    install_requires=requirements,
    license="MIT license",
    long_description=readme + '\n\n' + history,
    include_package_data=True,
    keywords='alifedata-phyloinformatics-convert',
    name='alifedata-phyloinformatics-convert',
    packages=find_packages(include=['alifedata_phyloinformatics_convert', 'alifedata_phyloinformatics_convert.*']),
    setup_requires=setup_requirements,
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/mmore500/alifedata-phyloinformatics-convert',
    version='0.7.0',
    zip_safe=False,
)
