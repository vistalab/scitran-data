"""

XXX


"""

from setuptools import setup, find_packages

setup(
    name='scitran.data',
    version='0.0.1',
    packages=find_packages(),
    package_data={
        'scitran': ['data/*.json'],
    },
    namespace_packages=['scitran'],
)
