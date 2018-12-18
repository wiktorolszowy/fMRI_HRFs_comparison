#!/bin/bash

#! Wiktor Olszowy

module load AFNI/AFNI_18.0.11/

set -o errexit

#-Folder where output and errors from running 'sbatch' will be saved
mkdir out_err

part=1; export part; sbatch --array=1-772 slurm_submit.array.hphi; sleep 20
part=2; export part; sbatch --array=1-772 slurm_submit.array.hphi; sleep 20
part=3; export part; sbatch --array=1-772 slurm_submit.array.hphi; sleep 20
part=4; export part; sbatch --array=1-772 slurm_submit.array.hphi; sleep 20
