import os
from setuptools import setup, find_packages
import versioneer  # https://github.com/warner/python-versioneer

# Constants
THIS_DIR = os.path.dirname(__file__)

# Package description
with open("README.md", "r") as fh:
    long_description = fh.read()

# Package requirements helper
def read_requirements_from_file(filepath):
    '''Read a list of requirements from the given file and split into a
    list of strings. It is assumed that the file is a flat
    list with one requirement per line.
    :param filepath: Path to the file to read
    :return: A list of strings containing the requirements
    '''
    with open(filepath, 'rU') as req_file:
        return req_file.readlines()

setup_args = dict(
    install_requires=read_requirements_from_file(
        os.path.join(
            THIS_DIR,
            'requirements.txt')),
    tests_require=read_requirements_from_file(
        os.path.join(
            THIS_DIR,
            'requirements-dev.txt')))

# Author list
authors = [
    'Marshall McDonnell',
    'Peter Peterson',
    'Elliot Oram',
]

# Main setup
setup(
    name='mantid_total_scattering',
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    author=','.join(authors)
    author_email='mcondonnellmd@ornl.gov',
    url='https://github.com/neutrons/mantid_total_scattering',
    description='Mantid Total Scattering Reduction',
    long_description_content_type="text/markdown",
    license='GPL License (version 3)',
    entry_points={
      'console_scripts': [
          "mantidtotalscattering = total_scattering.cli:main"
      ]
    },
    packages=find_packages(),
    include_package_data=True,
    setup_requires=[],
    install_requires=setup_args['install_requires'],
    tests_require=setup_args['install_requires'] + setup_args['tests_require'],
    test_suite='tests',
)
