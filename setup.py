from setuptools import setup, find_packages

setup(
    name='dcTMD',
    version='0.2.0',
    description='analyse targeted MD data with dcTMD',
    #long_description=README,
    keywords=[
        'enhance ampling',
        'friction',
        'MD analysis',
    ],
    author='dieJaegerIn and taenzel',
    package_dir={'': 'src'},
    packages=find_packages(where='src'),
    include_package_data=True,
    python_requires='>=3.8',
    entry_points={
        'console_scripts': [
            'dcTMD = dcTMD.__main__:main',
        ],
    }
)