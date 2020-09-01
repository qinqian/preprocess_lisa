#!/bin/bash -ex

#wget -c https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
#bash Miniconda3-latest-Linux-x86_64.sh

conda install anaconda-client
conda create -n lisa anaconda python=3
source activate lisa

while read requirement; do conda install --yes $requirement; done < requirements.txt
conda install -c daler bedtools

# double check
pip install -r requirements.txt

