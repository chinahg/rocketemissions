#!/bin/bash
                                            

#SBATCH --time=48:00:00
#SBATCH --job-name="1000x2000"
#SBATCH --mail-user=chinahg@mit.edu
#SBATCH --mail-type=BEGIN,END
#SBATCH -o slurm-%j.out
#xSBATCH -e slurm-%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=normal
#SBATCH --mem=50000MB
#####################################

n=$11000
s=$22000
psi_init=$30.00001
phi_init=$40.00001
psi_mult=$51.01
phi_mult=$61.01
job_id=$SLURM_JOBID
set_T=$8imported #imported or ambient
set_u=$9imported #imported or ambient

echo "Cantera"
echo "Job ID: $job_id"
echo "Setting s = $s and n = $n"
echo "psi_init = $psi_init, phi_init = $phi_init"
echo "psi_mult = $psi_mult, phi_init = $phi_mult"
echo "Temperature = $set_T, Velocity = $set_u"

now=$(date +"%T")
echo "Start time : $now"

echo "Number of CPUs per task: $SLURM_CPUS_PER_TASK"
julia --project=/home/chinahg/GCresearch/rocketemissions/  /home/chinahg/GCresearch/rocketemissions/plumeSSME.jl $n $s $psi_init $phi_init $psi_mult $phi_mult $job_id $set_T $set_u

now=$(date +"%T")
echo "End time : $now"

#--threads $SLURM_CPUS_PER_TASK