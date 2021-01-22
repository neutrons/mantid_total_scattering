import os
from setuptools import setup, find_packages
import versioneer  # https://github.com/warner/python-versioneer

# Constants
THIS_DIR = os.path.dirname(__file__)

# Package description
with open("README.md", "r") as fh:
    readme = fh.read()

requirements = []

setup_requirements = ['pytest-runner', ]

test_requirements = ['pytest>=3', ]

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
    author=','.join(authors),
    author_email='mcondonnellmd@ornl.gov',
    url='https://github.com/neutrons/mantid_total_scattering',
    description='Mantid Total Scattering Reduction',
    long_description=readme,
    long_description_content_type="text/markdown",
    license='GPL License (version 3)',
    entry_points={
      'console_scripts': [
          "mantidtotalscattering = total_scattering.cli:main"
      ]
    },
    packages=find_packages(),
    include_package_data=True,
    setup_requires=setup_requirements,
    tests_require=test_requirements,
    test_suite='tests',
)
