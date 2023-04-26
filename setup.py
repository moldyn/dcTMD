# -*- coding: utf-8 -*-
import pathlib
from collections import defaultdict
from setuptools import setup, find_packages


def get_extra_requirements(path, add_all=True):
    """Parse extra-requirements file."""
    with open(path) as depfile:
        extra_deps = defaultdict(set)
        for line in depfile:
            if not line.startswith('#'):
                if ':' not in line:
                    raise ValueError(
                        f'Dependency in {path} not correct formatted: {line}',
                    )
                dep, tags = line.split(':')
                tags = {tag.strip() for tag in tags.split(',')}
                for tag in tags:
                    extra_deps[tag].add(dep)

        # add tag `all` at the end
        if add_all:
            extra_deps['all'] = {
                tag for tags in extra_deps.values() for tag in tags
            }
    return extra_deps


# The directory containing this file
HERE = pathlib.Path(__file__).parent
# The text of the README file
README = (HERE / 'README.md').read_text()

setup(
    name='dcTMD',
    version='0.3.0',
    description='Analyse targeted molecular dynamics data with dcTMD',
    long_description=README,
    long_description_content_type='text/markdown',
    keywords=[
        'enhanced sampling',
        'friction',
        'MD analysis',
    ],
    author='taenzel, dieJaegerIn, braniii, floWneffetS',
    url='https://github.com/moldyn/dcTMD',
    license='MIT License',
    classifiers=[
        'License :: OSI Approved :: MIT License',
        'Intended Audience :: Science/Research',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Topic :: Scientific/Engineering :: Physics',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Software Development :: Libraries :: Python Modules',
    ],
    project_urls={
        'Documentation': 'https://moldyn.github.io/dcTMD',
        'Source Code': 'https://github.com/moldyn/dcTMD',
        'Bug Tracker': 'https://github.com/moldyn/dcTMD/issues',
    },
    package_dir={'': 'src'},
    packages=find_packages(where='src'),
    include_package_data=True,
    python_requires='>=3.8',
    install_requires=[
        'numpy>=1.21.0',
        'scikit-learn',
        'beartype>=0.10.4',
        'scipy',
        'tqdm',
        'click>=7.0.0',
        'typing_extensions>=3.9.0;python_version<"3.9"',
        'matplotlib>=3.7',
    ],
    extras_require=get_extra_requirements('extra-requirements.txt'),
    entry_points={
        'console_scripts': [
            'dcTMD = dcTMD.__main__:main',
        ],
    }
)
