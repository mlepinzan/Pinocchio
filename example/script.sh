#!/bin/bash
#SBATCH -J HEFFTE
#SBATCH -p boost_usr_prod
##SBATCH --qos=boost_qos_dbg
##SBATCH --qos=boost_qos_bprod
#SBATCH -A INA24_C6B04
#SBATCH -N 32
#SBATCH --ntasks-per-node=4
#SBATCH -c 8
#SBATCH --mem=0
#SBATCH --exclusive
#SBATCH --gres=gpu:4
#SBATCH -o hinocchio.out
#SBATCH -e hinocchio.err
#SBATCH --time=00:30:00

source /leonardo_work/INA24_C6B04/ginocchio.sh

cd ../src
make clean && make
cp hinocchio.x ../example/hinocchio_run.x

cd ../example

export OMP_PLACES=cores
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}

mpirun -np ${SLURM_NTASKS} --map-by ppr:${SLURM_NTASKS_PER_NODE}:node:pe=${OMP_NUM_THREADS} ./hinocchio_run.x parameter_file
