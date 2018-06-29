#!/bin/bash

# install python
if [[ $TRAVIS_OS_NAME == "linux" ]]; then
    wget https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh -O miniconda.sh
elif [[ $TRAVIS_OS_NAME == "osx" ]]; then
    wget https://repo.continuum.io/miniconda/Miniconda2-latest-MacOSX-x86_64.sh -O miniconda.sh
fi

bash miniconda.sh -b -p $HOME/miniconda
hash -r
conda config --set always_yes yes --set changeps1 no
conda update -q conda
conda info -a
pip install --upgrade pip
