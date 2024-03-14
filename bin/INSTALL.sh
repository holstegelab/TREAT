echo "** Installing TREAT and Otter"

echo "**** Installing Conda environment for TREAT"

conda env create -f treat.yml

echo "**** Conda installation done"

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
git clone https://github.com/holstegelab/otter.git && cd otter
mkdir build
make packages
make

echo "** Activating treat environment"
conda activate treat

echo "** Exporting packages to main path"
cd ..
export PATH=$PWD/:$PATH
export PATH=$PWD/otter/build/:$PATH
