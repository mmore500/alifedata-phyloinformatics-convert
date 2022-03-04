#!/usr/bin/env python

"""The setup script."""

from setuptools import setup, find_packages

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = [
    "dendropy>=4.5.2",
    "iterpop>=0.4.0",
    "lyncs-utils>=0.3.2",
    "opytional>=0.1.0",
    "pandas>=1.1.2",
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
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
    ],
    description="alifedata-phyloinformatics-conversion helps apply traditional phyloinformatics software to alife standardized data",
    install_requires=requirements,
    license="MIT license",
    long_description=readme + '\n\n' + history,
    include_package_data=True,
    keywords='alifedata-phyloinformatics-conversion',
    name='alifedata-phyloinformatics-conversion',
    packages=find_packages(include=['alifedata-phyloinformatics-conversion', 'alifedata-phyloinformatics-conversion.*']),
    setup_requires=setup_requirements,
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/mmore500/alifedata-phyloinformatics-conversion',
    version='0.0.0',
    zip_safe=False,
)
