#!/bin/bash
#SBATCH --job-name=scProcess
#SBATCH --account pedrini.edoardo
#SBATCH --mem=64GB  # amout of RAM in MB required (and max ram available).
#SBATCH --time=INFINITE  ## OR #SBATCH --time=10:00 means 10 minutes OR --time=01:00:00 means 1 hour
#SBATCH --ntasks=6  # number of required cores
#SBATCH --nodes=1  # not really useful for not mpi jobs
#SBATCH --mail-type=FAIL ## BEGIN, END, FAIL or ALL
#SBATCH --mail-user=pedrini.edoardo@hsr.it
#SBATCH --error="02_integrationSkip_harmony.err"
#SBATCH --output="02_integrationSkip_harmony.out"

echo "my job strart now" > 02_integrationSkip_harmony.log;

date >> 02_integrationSkip_harmony.log;

. /home/pedrini.edoardo/miniconda3/bin/activate;
conda activate env_R4.2_recover2;

Rscript /beegfs/scratch/ric.cosr/pedrini.edoardo/scRNAseq_Brainsphere_Absinta/BS_run_04/analysis/R_analysis/scr/02_integrationSkip_harmony.R

date >> 02_integrationSkip_harmony.log;
echo "all done!!" >> 02_integrationSkip_harmony.log