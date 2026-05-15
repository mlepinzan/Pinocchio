#!/bin/bash

#SBATCH --job-name=Test_NewClusterMock_GPU_version
#SBATCH --partition=boost_usr_prod
#SBATCH --account=INA24_C7B09   
#SBATCH --nodes=16
#SBATCH --ntasks=256             # Total MPI tasks 
#SBATCH --ntasks-per-node=16     # Explicitly specify MPI tasks per node
#SBATCH --cpus-per-task=2        # Number of OpenMP threads per MPI task
#SBATCH --gres=gpu:4             # Request 4 GPUs per node
#SBATCH --gpu-bind=closest       # Bind MPI tasks to the closest GPUs
##SBATCH --gpu-bind=map_gpu:0,0,0,0,1,1,1,1,2,2,2,2,3,3,3,3
##SBATCH --gpu-bind=map_gpu:0,0,1,1,2,2,3,3
##SBATCH --qos=boost_qos_bprod
#SBATCH --exclusive
#SBATCH --output=PLarge-%j.out
#SBATCH --error=PLarge-%j.err
#SBATCH --time=5:00:00

ml gcc/12.2 binutils nvhpc/24.3 openmpi/4.1.6--nvhpc--24.3 cuda/12.3

# edit the libraries' path
LIB=/leonardo_scratch/fast/CNHPC_1498509/lib/nvhpc-23.11
FTTW_LIB=${LIB}/fftw/fftw-3.3.10/lib
GSL_LIB=${LIB}/gsl/gsl-2.7.1/lib
PFFT_LIB=${LIB}/pfft/pfft/lib
PMT_LIB=/leonardo_scratch/fast/CNHPC_1498509/lib/pmt/local/lib64

LIB_DIR=/leonardo/pub/userexternal/dgoz0000/libraries
HEFF_LIBR=${LIB_DIR}/heffte_nvidia/lib64

LIB_Mass_shell=/leonardo_scratch/large/userexternal/tbatalha/Pinocchio/dep
HEALPIX_LIB=${LIB_Mass_shell}/Healpix_3.83/lib
CFITSIO_LIB=${LIB_Mass_shell}/lib

export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${FTTW_LIB}:${GSL_LIB}:${PFFT_LIB}:${PMT_LIB}:${HEFF_LIBR}:${HEALPIX_LIB}:${CFITSIO_LIB}

export OMP_TARGET_OFFLOAD=mandatory
export OMP_PROC_BIND=true
export OMP_WAIT_POLICY=ACTIVE
export OMP_PLACES=cores

####################################

# === Set Compilation Flags === #
export COMPILER_CC=mpicc
export COMPILER_CPP=mpic++
export DEBUG=NO
export OMP=YES
export GPU=NO       
export FULL_GPU=YES
export SCOREP=NO
export NSIGHT=NO
export ENERGY_CPU=NO
export ENERGY_GPU=NO

# === Set Directories === #
BASE_DIR="/leonardo_work/CNHPC_1498509/PINOCCHIO_final_GPU"
COMPILE_DIR="${BASE_DIR}/src" # GPU

# === Define Box Sizes & Grid Sizes === #
BOX_SIZES=(3870)
N_PARTICLES=(2160)
MAX_MEM=(31000)
MAX_MEM_PER_PARTICLE=(350)


# === Determine Execution Mode (CPU or GPU) === #
if [[ "$GPU" == "YES" || "$FULL_GPU" == "YES" ]]; then
    EXECUTABLE="Pinocchio_GPU"
    RUN_TYPE="GPU"
else
    EXECUTABLE="Pinocchio_CPU"
    RUN_TYPE="CPU"
fi

# === Get Submission Directory === #
SUBMIT_DIR=$(pwd)  # Captures the directory where sbatch is run

# === Compile the Code === #
cd ${COMPILE_DIR} && make clean && make GPU=$GPU OMP=$OMP FULL_GPU=$FULL_GPU DEBUG=$DEBUG ENERGY_CPU=$ENERGY_CPU ENERGY_GPU=$ENERGY_GPU EXEC=$EXECUTABLE

# === Prepare Execution === #
PARAM_FILE="${SUBMIT_DIR}/$1"    # Assume first argument is the param file name
OUTPUTS_FILE="${SUBMIT_DIR}/$2"  # Assume second argument is the outputs file name

if [[ -z "$1" || -z "$2" ]]; then
    echo "Error: Parameter file or outputs file not provided!"
    exit 1
fi

# === Nested Loops: Box Size & Grid Size === #
for BOX_SIZE in "${BOX_SIZES[@]}"; do
    for N_PARTICLES in "${N_PARTICLES[@]}"; do
        
        echo "Preparing execution for BoxSize=${BOX_SIZE}, GridSize=${N_PARTICLES}"  # Debugging info

        # Create unique directory name with execution type, grid, box, MPI, and OMP settings
        OUTPUT_DIR="${BASE_DIR}/Test_SBI_NewCLusterMock/${RUN_TYPE}_grid_${N_PARTICLES}_box_${BOX_SIZE}_MPI_${SLURM_NTASKS}_OMP_${SLURM_CPUS_PER_TASK}_test_0"
        mkdir -p $OUTPUT_DIR  # Create directory for this configuration
        if [[ ! -d "$OUTPUT_DIR" ]]; then
            echo "Error: Failed to create output directory $OUTPUT_DIR"
            continue
        fi

        # Copy base parameter and outputs file
        PARAM_FILE_RUN="${OUTPUT_DIR}/param_file_${N_PARTICLES}_${BOX_SIZE}"
        cp $PARAM_FILE $PARAM_FILE_RUN  
        cp $OUTPUTS_FILE "${OUTPUT_DIR}/outputs"
         
        # Check if parameter file was copied successfully
        if [[ ! -f "$PARAM_FILE_RUN" ]]; then
            echo "Error: Failed to copy parameter file to $PARAM_FILE_RUN"
            continue
        fi

        # Modify GridSize and BoxSize in the copied parameter file
        sed -i "s/^GridSize\s\+[0-9]\+/GridSize                   ${N_PARTICLES}/" $PARAM_FILE_RUN
        sed -i "s/^BoxSize\s\+[0-9]\+/BoxSize                     ${BOX_SIZE}/" $PARAM_FILE_RUN
        sed -i "s/^MaxMem\s\+[0-9]\+/MaxMem                       ${MAX_MEM}/" $PARAM_FILE_RUN
        sed -i "s/^MaxMemPerParticle\s\+[0-9]\+/MaxMemPerParticle ${MAX_MEM_PER_PARTICLE}/" $PARAM_FILE_RUN

        # Verify modifications were applied
        grep -E "^GridSize|^BoxSize" $PARAM_FILE_RUN

        # Move to RUN_DIR and execute
        cd $OUTPUT_DIR
        export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
        echo "Running: mpirun -n $SLURM_NTASKS $COMPILE_DIR/$EXECUTABLE $PARAM_FILE_RUN"
        mpirun -n $SLURM_NTASKS --map-by ppr:$SLURM_NTASKS_PER_NODE:node:pe=$OMP_NUM_THREADS $COMPILE_DIR/$EXECUTABLE $PARAM_FILE_RUN |& tee log_run_grid_${N_PARTICLES}__box_${BOX_SIZE}.txt
    done
done