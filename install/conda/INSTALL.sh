#!/bin/bash

echo "** Installing TREAT and Otter"

echo "**** Installing Conda environment for TREAT"

conda env create -f treat.yml

echo "**** TREAT successfully installed"

chmod +x ../../bin/TREAT.py
conda init bash
conda activate newtreat

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

echo "export PATH=../../bin/otter/build/:$PATH" >> activate_env.sh
echo "export PATH+=:../../bin/:$PATH" >> activate_env.sh
chmod +x ../../bin/activate_env.sh
mv activate_env.sh $CONDA_PREFIX/etc/conda/activate.d/
