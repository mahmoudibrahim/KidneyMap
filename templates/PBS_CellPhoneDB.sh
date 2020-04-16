#!/bin/sh
# Resources, e.g. a total time of 48 hours...
#PBS -l walltime=672:00:00
# CPU time
#PBS -l cput=3000:00:00
# Resources, ... and one node with 4 processors:
#PBS -l nodes=1:ppn=16
# Ensures that the Linux environment for the job is the same as the one we're working in:
#PBS -V
# stderr redirection
#PBS -e stderr.txt
# stdout redirection
#PBS -o stdout.txt

# Change directory to the project
if [ ! -z ${PBS_O_WORKDIR+x} ];then
  cd ${PBS_O_WORKDIR};
fi

source ~/.bashrc
conda activate CellPhoneDB
cellphonedb method statistical_analysis cellphonedb_meta.txt cellphonedb_count.txt
