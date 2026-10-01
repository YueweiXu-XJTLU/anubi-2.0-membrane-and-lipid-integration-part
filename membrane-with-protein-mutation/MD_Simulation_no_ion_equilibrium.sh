#!/bin/bash
#SBATCH -J  gpu_number_test
#SBATCH --partition=gpu4090
#SBATCH --qos=4gpus
#SBATCH --cpus-per-task=4
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4     # 1-4
#SBATCH --gres=gpu:4            # 1-4
#SBATCH -e  job.%j.err
#SBATCH -o  job.%j.out
#SBATCH --mail-user=Yuewei.Xu23@student.xjtlu.edu.cn
#SBATCH --mail-type=ALL

set -euo pipefail

ml load gromacs/2023.2-gcc-9.5.0-jzxesel

# EM
gmx_mpi grompp -f minim.mdp -c membrane_amber_lipid_neutral.gro -p topol.top -o minim.tpr -maxwarn 1
gmx_mpi mdrun -deffnm minim -v

# NVT
gmx_mpi grompp -f NVT.mdp -c minim.gro -r minim.gro -p topol.top -o nvt.tpr -maxwarn 1
srun -n ${SLURM_NTASKS} \
gmx_mpi mdrun -deffnm nvt -v \
-ntomp 4 -nb gpu -pme gpu -bonded gpu -update gpu

# NPT
gmx_mpi grompp -f NPT.mdp -c nvt.gro -r nvt.gro -p topol.top -o npt.tpr -maxwarn 1
srun -n ${SLURM_NTASKS} \
gmx_mpi mdrun -deffnm npt -v \
-ntomp 4 -nb gpu -pme gpu -bonded gpu -update gpu

# Production
gmx_mpi grompp -f md.mdp -c npt.gro -r npt.gro -p topol.top -o md.tpr -maxwarn 1
srun -n ${SLURM_NTASKS} \
gmx_mpi mdrun -deffnm md -v \
-ntomp 4 -nb gpu -pme gpu -bonded gpu -update gpu