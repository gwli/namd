TOP=`pwd`/..
BUILD="$TOP/Linux-x86_64-g++"
exec_proc=$BUILD/namd2
data=$BUILD/src/alanin
cmd="$exec_proc $data"
sp="/raid/tools/SP/NsightSystems-linux-public-2018.0.1.61-df7ba6e/Target-x86_64/x86_64/nsys"


apt update && apt install -y libfftw3-dev libfftw3-mpi-dev tcl8.6 tcl8.6-dev sfftw-dev

cd $BUILD

rm -fr ./*.qdstrm
$sp profile --delay=5 --duration=10 -o namd.qdstrm -t cuda,cublas,cudnn,osrt $cmd
mv namd.qdstrm $TOP/devtools
