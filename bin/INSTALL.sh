echo "** Installing TREAT and Otter"

echo "**** Installing Conda environment for TREAT"

conda env create -f newtreat.yml

echo "**** Conda installation done"

conda activate newtreat

echo "**** Now installing HTS lib"

wget https://github.com/samtools/htslib/releases/download/1.19.1/htslib-1.19.1.tar.bz2 .
tar -xf htslib-1.19.1.tar.bz2
cd htsblib-1.19.1
./configure
make
make install





