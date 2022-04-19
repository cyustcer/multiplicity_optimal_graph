#!/bin/bash
#BSUB -q short
#BSUB -J R_job[1-295]
#BSUB -o ../outputs/logs/grid-%I.out
#BSUB -e ../outputs/logs/grid-%I.err
#BSUB -W 480
#BSUB -M 10000

# source /CHBS/apps/busdev_apps/init.sh
export PATH="/CHBS/apps/EB/software/R/4.1.0-gomkl-2019a/lib64/R/bin/:$PATH"

# Rscript 01_simulate_data.R --number ${LSB_JOBINDEX} --margin_power "0.5, 0.8, 0.9" --w_step 0.01 --G_step 0.2 --w_range "0, 1" --batch_size 4000

Rscript 01_simulate_data_G.R --number ${LSB_JOBINDEX} --margin_power "0.7, 0.8, 0.9" --w_specified "0.2147974, 0.3144578, 0.4707447" --G_step 0.01 --batch_size 3500
