#!/bin/bash
#SBATCH -n 8
#SBATCH --job-name=kim1
#SBATCH --output=kim1.out
#SBATCH --error=kim1.err

java -jar parsePhyloBayes.jar 140926 167 kim1 2 root.outgroup.txt 