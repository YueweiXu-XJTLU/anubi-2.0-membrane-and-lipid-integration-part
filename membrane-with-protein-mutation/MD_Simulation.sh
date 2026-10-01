#!/bin/bash
#SBATCH -J  membrane_protein_mutation
#SBATCH --partition=gpu4090         # Partition name
#SBATCH -n 1                         # Number of tasks (typically 1 if not using MPI)
#SBATCH --qos=1gpu               # Quality of Service (QOS) setting
#SBATCH --cpus-per-task=4       # 1:4 GPU:CPU ratio
#SBATCH --gres=gpu:1
#SBATCH -e  job.%j.err
#SBATCH -o  job.%j.out
#SBATCH --mail-user=Yuewei.Xu23@student.xjtlu.edu.cn
#SBATCH --mail-type=ALL

set -euo pipefail

ml load gromacs/2023.2-gcc-9.5.0-jzxesel

# Equilibrate Ions
gmx_mpi grompp -f minim.mdp -c membrane_amber_lipid.gro -p topol.top -o neutral.tpr -maxwarn 1
printf "SOL\n" | gmx_mpi genion -s neutral.tpr -o membrane_amber_lipid_neutral.gro -p topol.top -pname K -nname CL -neutral

# EM
gmx_mpi grompp -f minim.mdp -c membrane_amber_lipid_neutral.gro -p topol.top -o minim.tpr -maxwarn 1
gmx_mpi mdrun -deffnm minim -v

# NVT
gmx_mpi grompp -f NVT.mdp -c minim.gro -r minim.gro -p topol.top -o nvt.tpr -maxwarn 1
gmx_mpi mdrun -deffnm nvt -v -gpu_id 0 -nb gpu -pme gpu -bonded gpu -update gpu

# NPT
gmx_mpi grompp -f NPT.mdp -c nvt.gro -r nvt.gro -p topol.top -o npt.tpr -maxwarn 1
gmx_mpi mdrun -deffnm npt -v -gpu_id 0 -nb gpu -pme gpu -bonded gpu -update gpu

# Production
gmx_mpi grompp -f md.mdp -c npt.gro -r npt.gro -p topol.top -o md.tpr -maxwarn 1
gmx_mpi mdrun -deffnm md -v -gpu_id 0 -nb gpu -pme gpu -bonded gpu -update gpu