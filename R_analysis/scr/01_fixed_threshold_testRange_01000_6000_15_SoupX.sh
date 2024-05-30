#!/bin/bash
#SBATCH --job-name=scProcess
#SBATCH --account pedrini.edoardo
#SBATCH --mem=128GB  # amout of RAM in MB required (and max ram available).
#SBATCH --time=INFINITE  ## OR #SBATCH --time=10:00 means 10 minutes OR --time=01:00:00 means 1 hour
#SBATCH --ntasks=8  # number of required cores
#SBATCH --nodes=1  # not really useful for not mpi jobs
#SBATCH --mail-type=FAIL ## BEGIN, END, FAIL or ALL
#SBATCH --mail-user=pedrini.edoardo@hsr.it
#SBATCH --error="01_fixed_threshold_testRange_01000_6000_15_SoupX.err"
#SBATCH --output="01_fixed_threshold_testRange_01000_6000_15_SoupX.out"

echo "my job strart now" > 01_fixed_threshold_testRange_01000_6000_15_SoupX.log;

date >> 01_fixed_threshold_testRange_01000_6000_15_SoupX.log;

. /home/pedrini.edoardo/miniconda3/bin/activate;
conda activate env_R4.2_recover2;

Rscript /beegfs/scratch/ric.cosr/pedrini.edoardo/scRNAseq_Brainsphere_Absinta/BS_run_04/analysis/R_analysis/scr/01_fixed_threshold_testRange_01000_6000_15_SoupX.R

date >> 01_fixed_threshold_testRange_01000_6000_15_SoupX.log;
echo "all done!!" >> 01_fixed_threshold_testRange_01000_6000_15_SoupX.log