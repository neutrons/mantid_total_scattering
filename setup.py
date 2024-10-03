import os
from setuptools import setup, find_packages
import versioneer  # https://github.com/warner/python-versioneer

# Redirect stdout and stderr to null
os.close(1)
os.close(2)
os.open(os.devnull, os.O_RDWR)
os.dup(1)
os.dup(2)

# Constants
THIS_DIR = os.path.dirname(__file__)

# Package description
with open("README.md", "r") as fh:
    readme = fh.read()


# Package requirements helper
def read_requirements_from_file(filepath):
    '''Read a list of requirements from the given file and split into a
    list of strings. It is assumed that the file is a flat
    list with one requirement per line.
    :param filepath: Path to the file to read
    :return: A list of strings containing the requirements
    '''
    with open(filepath, 'r') as req_file:
        lines = req_file.readlines()
        lines_return = list()
        for line in lines:
            if "https://" in line:
                lines_return.append(f"pyoncat @ {line}")
            else:
                lines_return.append(line)
        return lines_return


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
    'Yuanpeng Zhang',
    'Marshall McDonnell',
    'Peter Peterson',
    'Coleman Kendrick',
    'Jenna DeLozier',
    'Elliot Oram',
    'Donnie Earnest',
]

# Main setup
setup(
    name='mantid_total_scattering',
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    author=','.join(authors),
    author_email='mcondonnellmd@ornl.gov',
    url='https://github.com/neutrons/mantid_total_scattering',
    description='Mantid Total Scattering Reduction',
    long_description=readme,
    long_description_content_type="text/markdown",
    license='GPL License (version 3)',
    packages=find_packages(),
    include_package_data=True,
    setup_requires=[],
    install_requires=setup_args['install_requires'],
    tests_require=setup_args['install_requires'] + setup_args['tests_require'],
    test_suite='tests',
    entry_points={
        'console_scripts': [
            "mantidtotalscattering = total_scattering.cli:main"
        ]
    },
)
