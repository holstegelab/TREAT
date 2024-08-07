#!/bin/bash

echo "** Installing TREAT and Otter"

echo "**** Installing Conda environment for TREAT"

conda env create -f treat.yml

echo "**** TREAT successfully installed"

chmod +x $PWD/TREAT.py
conda activate treat

echo "**** Now installing HTS lib"

wget https://github.com/samtools/htslib/releases/download/1.19.1/htslib-1.19.1.tar.bz2 .
tar -xvf htslib-1.19.1.tar.bz2
cd htslib-1.19.1
./configure --prefix=$PWD
make
make install

echo "**** Installing Otter"
cd ..
export CPATH=$CPATH:$PWD/htslib-1.19.1
export LIBRARY_PATH=$LIBRARY_PATH:$PWD/htslib-1.19.1
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PWD/htslib-1.19.1
echo export CPATH=$CPATH:$PWD/htslib-1.19.1 > $PWD/activate_env.sh
echo export LIBRARY_PATH=$LIBRARY_PATH:$PWD/htslib-1.19.1 >> $PWD/activate_env.sh
echo export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PWD/htslib-1.19.1 >> $PWD/activate_env.sh
git clone https://github.com/holstegelab/otter.git && cd otter
mkdir build
make packages
make

echo "** Exporting packages to main path"
cd ../../

echo "export PATH=$PWD/bin/otter/build/:$PATH" >> $PWD/bin/activate_env.sh
echo "export PATH+=:$PWD/bin/:$PATH" >> $PWD/bin/activate_env.sh
chmod +x $PWD/bin/activate_env.sh
mv $PWD/bin/activate_env.sh $CONDA_PREFIX/etc/conda/activate.d/
