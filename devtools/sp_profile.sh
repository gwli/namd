TOP=`pwd`/..

BUILD="$TOP/Linux-x86_64-g++"
exec_proc=$BUILD/namd2
data=$BUILD/src/alanin
cmd="$exec_proc $data"
sp="/raid/tools/SP/NsightSystems-linux-public-2018.0.0.170-5fe72bc/Target-x86_64/x86_64/sp"

cd $BUILD

rm -fr ./*.qdstrm
$sp profile --duration=10 -o lmmps.qdstrm -t cuda,cublas,curand,cudnn,osrt $cmd
