from setuptools import setup, find_packages
# import versioneer  # https://github.com/warner/python-versioneer

setup(name='mantid_total_scattering',
      version=1.0,  # versioneer.get_version(),
      # cmdclass=versioneer.get_cmdclass(),
      description='Need a description',
      author='Marshall, Elliot, Pete',
      author_email='mcondonnellmd@ornl.gov',
      url='https://github.com/marshallmcdonnell/mantid_total_scattering',
      long_description='''Should have a longer description''',
      license='GPL License (version 3)',
      scripts=['bin/mantidtotalscattering'],
      packages=find_packages(),
      # package_data={'': ['*.ui', '*.png', '*.qrc', '*.json']},
      include_package_data=True,
      install_requires=[],
      setup_requires=[],
      test_suite='tests'
      )
