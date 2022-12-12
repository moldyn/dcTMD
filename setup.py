from setuptools import setup, find_packages

setup(
    name='dcTMD',
    version='0.2.1',
    description='Analyse targeted molecular dynamics data with dcTMD',
    #long_description=README,
    keywords=[
        'enhance sampling',
        'friction',
        'MD analysis',
    ],
    author='taenzel, dieJaegerIn, floWneffetS',
    project_urls={
        #'Documentation': 'https://moldyn.github.io/MoSAIC',
        'Source Code': 'https://github.com/moldyn/dcTMD',
        #'Bug Tracker': 'https://github.com/moldyn/MoSAIC/issues',
    },
    package_dir={'': 'src'},
    packages=find_packages(where='src'),
    include_package_data=True,
    python_requires='>=3.8',
        install_requires=[
        'numpy>=1.21.0',
        'matplotlib',
        'beartype>=0.10.4',
        'scipy',
        'tqdm',
        'click>=7.0.0',
        'typing_extensions>=3.9.0;python_version<"3.9"',
    ],
    entry_points={
        'console_scripts': [
            'dcTMD = dcTMD.__main__:main',
        ],
    }
)