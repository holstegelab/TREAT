#!/bin/bash

echo "** Installing TREAT and Otter"

echo "**** Installing Conda environment for TREAT"

conda env create -f treat.yml

echo "**** TREAT successfully installed"

echo "**** Activating TREAT environment"
conda init bash
source ~/.bash_profile
conda activate newtreat

chmod +x ../../bin/TREAT.py

echo "**** Now installing Otter"

cd ../../bin/
git clone --branch development --recursive https://github.com/holstegelab/otter.git
cd otter/include/WFA2-lib
make clean setup lib_wfa
cd ../../
mkdir build
make packages
make
cd ..

echo "export PATH=${PWD}/otter/build/:$PATH" >> activate_env.sh
echo "export PATH+=:${PWD}/:$PATH" >> activate_env.sh
chmod +x activate_env.sh
mkdir -p ${CONDA_PREFIX}/etc/conda/activate.d/
mv activate_env.sh $CONDA_PREFIX/etc/conda/activate.d/
