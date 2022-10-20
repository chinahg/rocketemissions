#!/bin/bash
                                            

#SBATCH --time=48:00:00
#SBATCH --job-name="2000x5000"
#SBATCH --mail-user=chinahg@mit.edu
#SBATCH --mail-type=BEGIN,END
#SBATCH -o slurm-%j.out
#xSBATCH -e slurm-%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --partition=normal
#SBATCH --mem=56000MB
#####################################

n=$11000 #2000
s=$25000 #5000
psi_init=$30.001 #0.001
phi_init=$40.001 #0.01
psi_mult=$51.01 #1.01
phi_mult=$61.01 #1.01
job_id=$7#SLURM_JOB_ID
echo "Setting s = $s and n = $n"
echo "psi_init = $psi_init, phi_init = $phi_init"
echo "psi_mult = $psi_mult, phi_init = $phi_mult"

now=$(date +"%T")
echo "Start time : $now"

echo "Number of CPUs per task: $SLURM_CPUS_PER_TASK"
julia --project=/home/chinahg/GCresearch/rocketemissions/ --threads $SLURM_CPUS_PER_TASK /home/chinahg/GCresearch/rocketemissions/plumeSSME.jl $n $s $psi_init $phi_init $psi_mult $phi_mult $job_id

now=$(date +"%T")
echo "End time : $now"

#