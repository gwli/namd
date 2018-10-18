tagname="qafarm-sh01:5000/namd"
curr=`pwd`
TOP="/raid/HPC/namd"

cd $TOP
echo $curr

#cd $TOP
echo `pwd`
nvidia-docker build -t $tagname -f $curr/Dockerfile . && \
docker push  $tagname


