#! /bin/bash
#-----------------------------------------
#SBATCH --job-name=Creating_kinship_matrix_weighted

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=3
#SBATCH --mem=128GB
#SBATCH --time=1-00:00:00
#SBATCH --partition=bigmem

#-----------------------------------------

java -Xmx2g -jar /work/long_lab/Rushani/EDLMM/Codes/jawamix5.jar kinshipwe -ig /work/long_lab/Rushani/IEDLMM/data/BPD/BDO_BARD_GRU_merge.match.b38.hdf5 -o /work/long_lab/Rushani/IEDLMM/data/BPD/BPD_8_kinship -wg /work/long_lab/Rushani/IEDLMM/weight_files/BPD/BPD_8_weights_adjusted.txt -m RRM

