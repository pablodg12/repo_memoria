#!/bin/bash
#SBATCH -o generator25_o_out-%J.txt
#SBATCH -e generator25_e_out-%J.txt
#SBATCH --partition=general
#SBATCH --time=240:00:00
#SBATCH -n 1
#SBATCH --ntasks-per-node=1
cd $PBS_O_WORKDIR
module load Python/3.7.3 
python -u /home/pibarra/generador_cadenas/repo_memoria/generador_cadenas/main_restart.py ${1} ${2} ${3} ${4} ${5}
