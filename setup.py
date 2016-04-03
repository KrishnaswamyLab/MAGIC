import os
import shutil
from subprocess import call
from setuptools import setup
from warnings import warn

# install phenograph
call(['pip3', 'install', 'git+https://github.com/jacoblevine/phenograph.git'])


setup(name='wishbone',
      version='1.0',
      description='Wishbone algorithm for identifying bifurcating trajectories from single-cell data',
      author='Manu Setty',
      author_email='manu.setty@columbia.edu',
      package_dir={'': 'src'},
      install_requires=[
          'numpy>=1.10.0',
          'pandas>=0.18.0',
          'scipy>=0.14.0',
          'tsne',
          'matplotlib',
          'seaborn',
          'sklearn',
          'networkx',
          'fcsparser',
          'statsmodels'],
      )


# get location of setup.py
setup_dir = os.path.dirname(os.path.realpath(__file__))

# install GSEA, diffusion components
tools_dir = os.path.expanduser('~/.wishbone/tools')
if os.path.isdir(tools_dir):
    shutil.rmtree(tools_dir)
shutil.copytree(setup_dir + '/tools/', tools_dir)
shutil.unpack_archive(tools_dir + '/DiffusionGeometry.zip', tools_dir +
                      '/DiffusionGeometry/')
shutil.unpack_archive(tools_dir + '/mouse_gene_sets.tar.gz', tools_dir)
shutil.unpack_archive(tools_dir + '/human_gene_sets.tar.gz', tools_dir)

# Copy test data
data_dir = os.path.expanduser('~/.wishbone/data')
if os.path.isdir(data_dir):
    shutil.rmtree(data_dir)
shutil.copytree(setup_dir + '/data/', data_dir)

# Create directory for GSEA reports
os.makedirs( os.path.expanduser('~/.wishbone/gsea/'), exist_ok=True )