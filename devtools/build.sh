TOP=`pwd`/..

CHARM=$TOP/charm


cd $CHARM
./build charm++ multicore-linux-x86_64 --with-production
#./build charm++ mpi-linux-x86_64 --with-production

make -j

apt update && apt install -y libfftw3-dev libfftw3-mpi-dev tcl8.6 tcl8.6-dev sfftw-dev


cd $TOP
rm -fr Linux-x86_64-g++
./config Linux-x86_64-g++ --charm-arch multicore-linux-x86_64 --with-cuda
cd $TOP/Linux-x86_64-g++
make -j 
