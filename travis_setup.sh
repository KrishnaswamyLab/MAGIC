#!/bin/bash

# install python
if [[ $TRAVIS_OS_NAME == "linux" ]]; then
    wget https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh -O miniconda2.sh
    wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda3.sh
elif [[ $TRAVIS_OS_NAME == "osx" ]]; then
    wget https://repo.continuum.io/miniconda/Miniconda2-latest-MacOSX-x86_64.sh -O miniconda2.sh
    wget https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh -O miniconda3.sh
fi

bash miniconda2.sh -b -p $HOME/miniconda2
hash -r
$HOME/miniconda2/bin/conda config --set always_yes yes --set changeps1 no
$HOME/miniconda2/bin/conda update -q conda
$HOME/miniconda2/bin/conda info -a
$HOME/miniconda2/bin/pip install --upgrade pip

bash miniconda3.sh -b -p $HOME/miniconda3
hash -r
$HOME/miniconda3/bin/conda config --set always_yes yes --set changeps1 no
$HOME/miniconda3/bin/conda update -q conda
$HOME/miniconda3/bin/conda info -a
$HOME/miniconda3/bin/pip install --upgrade pip
