import os
import sys
from setuptools import setup

install_requires = [
    'numpy>=1.14.0',
    'pandas>=0.21.0',
    'scipy>=1.1.0',
    'matplotlib',
    'scikit-learn>=0.19.1',
    'tasklogger>=0.2.1',
    'graphtools>=0.1.9',
    'scprep>=0.7.1'
]

test_requires = [
    'nose2',
]

if sys.version_info[0] == 3:
    test_requires += ['anndata']

doc_requires = [
    'sphinx',
    'sphinxcontrib-napoleon',
]

if sys.version_info[:2] < (2, 7) or (3, 0) <= sys.version_info[:2] < (3, 5):
    raise RuntimeError("Python version 2.7 or >=3.5 required.")

version_py = os.path.join(os.path.dirname(
    __file__), 'magic', 'version.py')
version = open(version_py).read().strip().split(
    '=')[-1].replace('"', '').strip()

readme = open('README.rst').read()

setup(name='magic-impute',
      version=version,
      description='MAGIC',
      author='',
      author_email='',
      packages=['magic', ],
      license='GNU General Public License Version 2',
      install_requires=install_requires,
      extras_require={'test': test_requires,
                      'doc': doc_requires},
      test_suite='nose2.collector.collector',
      long_description=readme,
      url='https://github.com/KrishnaswamyLab/MAGIC',
      download_url="https://github.com/KrishnaswamyLab/MAGIC/archive/v{}.tar.gz".format(
          version),
      keywords=['visualization',
                'big-data',
                'dimensionality-reduction',
                'embedding',
                'manifold-learning',
                'computational-biology'],
      classifiers=[
          'Development Status :: 5 - Production/Stable',
          'Environment :: Console',
          'Framework :: Jupyter',
          'Intended Audience :: Developers',
          'Intended Audience :: Science/Research',
          'Natural Language :: English',
          'Operating System :: MacOS :: MacOS X',
          'Operating System :: Microsoft :: Windows',
          'Operating System :: POSIX :: Linux',
          'Programming Language :: Python :: 2',
          'Programming Language :: Python :: 2.7',
          'Programming Language :: Python :: 3',
          'Programming Language :: Python :: 3.5',
          'Programming Language :: Python :: 3.6',
          'Topic :: Scientific/Engineering :: Bio-Informatics',
      ]
      )
