#!/bin/bash

#$ -l hostname=cpu-2-0|cpu-2-1|cpu-2-2|cpu-2-3|cpu-2-4|cpu-2-5|cpu-2-6|cpu-2-7

#PBS -j oe
#PBS -q all.q
#PBS -l nodes=2:ppn=96
#PBS -l walltime=TRY2:00:00
#PBS -V

# Export these environmental variables
#$ -v CUDA_HOME=/usr/local/cuda-11.0/,AMBERHOME=/home/cuda/amber21-cuda11
#$ -v LD_LIBRARY_PATH=/usr/local/cuda-11.0/lib64:/home/cuda/amber21-cuda11/lib:/opt/gridengine/lib/linux-x64:/opt/openmpi/lib:/opt/python/lib/:/home/zhuqiang/software/gcc/9.5.0/lib64
#$ -v PATH=/usr/local/cuda-11.0/bin:/home/cuda/amber21-cuda11/bin:/bin/opt/openmpi/bin:/usr/lib64/qt-3.3/bin:/usr/local/bin:/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/sbin:/opt/bio/ncbi/bin:/opt/bio/mpiblast/bin:/opt/bio/EMBOSS/bin:/opt/bio/clustalw/bin:/opt/bio/tcoffee/bin:/opt/bio/hmmer/bin:/opt/bio/phylip/exe:/opt/bio/mrbayes:/opt/bio/fasta:/opt/bio/glimmer/bin:/opt/bio/glimmer/scripts:/opt/bio/gromacs/bin:/opt/bio/gmap/bin:/opt/bio/tigr/bin:/opt/bio/autodocksuite/bin:/opt/bio/wgs/bin:/opt/eclipse:/opt/ganglia/bin:/opt/ganglia/sbin:/usr/java/latest/bin:/opt/rocks/bin:/opt/rocks/sbin:/opt/condor/bin:/opt/condor/sbin:/opt/gridengine/bin/linux-x64
#$ -M yongxian.wu@uci.edu
#PBS -A account_no
#PBS -N pmemd.MPI

# Set the no of threads, should match the nodes value above.
#$ -pe smp 112

#########################################
# Macro setting for experiment related things
exp_name=$1
cwrk=$2
output_dir_name=$3
parmfile=$4
prepared_trajin_content=$5

echo "********************************************************************************"
echo "Current experiment setting:"
echo "exp_name=$exp_name"
echo "cwrk=$cwrk"
echo "output_dir_name=$output_dir_name"
echo "parmfile=$parmfile"
echo "trajin_content=$prepared_trajin_content"
echo "********************************************************************************"

#########################################
# Set the AMBER ENV variables

AMBERHOME=/home8/yxwu/AMBER/bin/amber23-yongxian-MPI

source $AMBERHOME/amber.sh
export PATH=/home8/yxwu/apps/bin/gcc-11.3.0/bin:$PATH
export LD_LIBRARY_PATH=/home8/yxwu/apps/bin/gcc-11.3.0/lib64:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/home8/yxwu/apps/bin/cuda_11.7.0/lib64:$LD_LIBRARY_PATH
export PATH=/home8/yxwu/apps/bin/openmpi-4.1.5/bin:$PATH

echo "AMBERHOME=$AMBERHOME"
echo "LD_LIBRARY_PATH=$LD_LIBRARY_PATH"
echo "HOSTNAME=$HOSTNAME"
echo "--------------------------------------------------------------------------------"
# Set the working directory
echo "Entering Working Directory"
#$ -cwd

# Set location of your original working directory
# current location of 01_min
echo $cwrk
cd $cwrk

#########################################
# Set the MD engine
exe=$AMBERHOME/bin/cpptraj

#########################################
# Create experiment running output dir
# mkdir $output_dir_name
cd $output_dir_name

echo "Current Working Directory"
echo $(pwd)

trajin_content="${prepared_trajin_content//##NEWLINE##/$'\n'}"

cat >mdin <<EOF
parm $parmfile

$trajin_content

radial out rdf_ghh.out 0.02 9.0 @H*
go
quit
EOF

echo "--------------------------------------------------------------------------------"
echo "cmd:"
echo "$exe -i mdin"
echo "--------------------------------------------------------------------------------"

$exe -i mdin

####################################
echo "Local job completed"
echo $(date)
exit
