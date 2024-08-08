#!/bin/bash

echo "** Installing TREAT and Otter"

echo "**** Installing Conda environment for TREAT"

conda env create -f treat.yml

echo "**** TREAT successfully installed"

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

conda init bash
source ~/.bash_profile
conda activate newtreat

echo "export PATH=${PWD}/otter/build/:$PATH" >> activate_env.sh
echo "export PATH+=:${PWD}/:$PATH" >> activate_env.sh
chmod +x activate_env.sh
mv activate_env.sh $CONDA_PREFIX/etc/conda/activate.d/

echo "**** Installation is over. If you noticed errors in the last commands, that means that the Conda environment was not activated successfully. It's OK, but you may need to manually add TREAT and Otter to your executables."
