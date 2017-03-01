import os
import sys
import shutil
from subprocess import call
from setuptools import setup
from warnings import warn

if sys.version_info.major != 3:
    raise RuntimeError('Magic requires Python 3')


# install phenograph
call(['pip3', 'install', 'git+https://github.com/jacoblevine/phenograph.git'])


setup(name='magic',
      version='0.0',
      description='MAGIC',
      author='',
      author_email='',
      package_dir={'': 'src'},
      packages=['magic'],
      install_requires=[
          'numpy>=1.10.0',
          'pandas>=0.18.0',
          'scipy>=0.14.0',
          'matplotlib',
          'seaborn',
          'sklearn',
          'networkx',
          'fcsparser',
          'statsmodels'],
      scripts=['src/magic/magic_gui.py'],
      )


# get location of setup.py
setup_dir = os.path.dirname(os.path.realpath(__file__))

# Copy test data
data_dir = os.path.expanduser('~/.magic/data')
if os.path.isdir(data_dir):
    shutil.rmtree(data_dir)
shutil.copytree(setup_dir + '/data/', data_dir)

