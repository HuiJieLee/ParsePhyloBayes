#!/bin/bash
#SBATCH -n 8
#SBATCH --job-name=kim2
#SBATCH --output=kim2.out
#SBATCH --error=kim2.err

java -jar parsePhyloBayes.jar 140926 167 kim2 2 root.outgroup.txt 