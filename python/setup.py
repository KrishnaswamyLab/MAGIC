import os
import sys
from setuptools import setup

install_requires = [
    "numpy>=1.14.0",
    "scipy>=1.1.0",
    "matplotlib",
    "scikit-learn>=0.19.1",
    "future",
    "tasklogger>=1.0.0",
    "graphtools>=1.4.0",
    "pandas>=0.25",
    "scprep>=1.0",
]

test_requires = [
    "nose2",
]

if sys.version_info[0] == 3:
    test_requires += ["anndata"]

doc_requires = [
    "sphinx",
    "sphinxcontrib-napoleon",
]

if sys.version_info[:2] < (3, 5):
    raise RuntimeError("Python version >=3.5 required.")
elif sys.version_info[:2] >= (3, 6):
    test_requires += ["black"]

version_py = os.path.join(os.path.dirname(__file__), "magic", "version.py")
version = open(version_py).read().strip().split("=")[-1].replace('"', "").strip()

readme = open("README.rst").read()

setup(
    name="magic-impute",
    version=version,
    description="MAGIC",
    author="",
    author_email="",
    packages=[
        "magic",
    ],
    license="GNU General Public License Version 2",
    install_requires=install_requires,
    extras_require={"test": test_requires, "doc": doc_requires},
    test_suite="nose2.collector.collector",
    long_description=readme,
    url="https://github.com/KrishnaswamyLab/MAGIC",
    download_url="https://github.com/KrishnaswamyLab/MAGIC/archive/v{}.tar.gz".format(
        version
    ),
    keywords=[
        "visualization",
        "big-data",
        "dimensionality-reduction",
        "embedding",
        "manifold-learning",
        "computational-biology",
    ],
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Environment :: Console",
        "Framework :: Jupyter",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "Natural Language :: English",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: Microsoft :: Windows",
        "Operating System :: POSIX :: Linux",
        "Programming Language :: Python :: 2",
        "Programming Language :: Python :: 2.7",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.5",
        "Programming Language :: Python :: 3.6",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
)
