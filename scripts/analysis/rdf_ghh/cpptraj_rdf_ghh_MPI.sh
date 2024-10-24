#!/bin/bash

#SBATCH --job-name=rdf_goo
#SBATCH --output=res_%j.out
#SBATCH --error=res_%j.err
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --exclusive

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

DO_PARALLEL="mpirun -np 32 --oversubscribe"

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
echo "$DO_PARALLEL $exe -i mdin"
echo "--------------------------------------------------------------------------------"

$DO_PARALLEL $exe -i mdin

####################################
echo "Local job completed"
echo $(date)
exit
