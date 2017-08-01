from setuptools import setup, find_packages

packages = ['.']
packages.extend(filter(lambda x: x.find('athena') == 0, find_packages()))
setup(
  name = "athena",
  version = "0.1",
  packages = packages,
  entry_points = {
    'console_scripts': [ 'athena = main:main' ]
  },
  install_requires=[
    'ipython-cluster-helper>=0.5.2',
    'pysam>=0.9',
  ],
  author='Alex Bishara',
  description='athena assembler',
)
