#!/bin/bash
#SBATCH --job-name=integ
#SBATCH --account pedrini.edoardo
#SBATCH --mem=128GB  # amout of RAM in MB required (and max ram available).
#SBATCH --time=INFINITE  ## OR #SBATCH --time=10:00 means 10 minutes OR --time=01:00:00 means 1 hour
#SBATCH --ntasks=20  # number of required cores
#SBATCH --nodes=1  # not really useful for not mpi jobs
#SBATCH --mail-type=FAIL ## BEGIN, END, FAIL or ALL
#SBATCH --mail-user=pedrini.edoardo@hsr.it
#SBATCH --error="03_label_transfer_BS_run01_02.err"
#SBATCH --output="03_label_transfer_BS_run01_02.out"

echo "my job strart now" > 03_label_transfer_BS_run01_02.log;

date >> 03_label_transfer_BS_run01_02.log;

. /home/pedrini.edoardo/miniconda3/bin/activate;
conda activate env_R4.2_renv;

Rscript /beegfs/scratch/ric.cosr/pedrini.edoardo/scRNAseq_Brainsphere_Absinta/BS_batch_2023/analysis/R_analysis/scr/03_label_transfer_BS_run01_02.R

date >> 03_label_transfer_BS_run01_02.log;
echo "all done!!" >> 03_label_transfer_BS_run01_02.log