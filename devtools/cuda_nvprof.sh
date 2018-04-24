TOP=`pwd`/..

BUILD="$TOP/Linux-x86_64-g++"
exec_proc=$BUILD/namd2
data=$BUILD/src/alanin
cmd="$exec_proc $data"

cd $BUILD
nvprof --csv $cmd
mvprof --query-events
nvprof --query-metrics
nvprof --events inst_executed $cmd
nvprof --metrics all $cmd
nvprof --metrics inst_executed $cmd
nvprof --aggregate-mode-off --source-level-analysis global_access,shared_access,branch,instruction_execution,pc_sampling  $cmd

