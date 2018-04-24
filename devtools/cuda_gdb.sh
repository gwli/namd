TOP=`pwd`/..
BUILD="$TOP/Linux-x86_64-g++"
exec_proc=$BUILD/namd2
data=$BUILD/src/alanin

cd $BUILD
cmd="$exec_proc $data"

cuda-gdb  --args $cmd
