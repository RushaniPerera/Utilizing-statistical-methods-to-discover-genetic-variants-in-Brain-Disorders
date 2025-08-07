#! /bin/bash
#-----------------------------------------
#SBATCH --job-name=Asscoiations

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=128GB
#SBATCH --time=1-00:00:00
#SBATCH --partition=bigmem

#-----------------------------------------

java -Xmx2g -jar /work/long_lab/Rushani/EDLMM/Codes/jawamix5.jar emmax -ig /work/long_lab/Rushani/IEDLMM/data/BPD/BDO_BARD_GRU_merge.match.b38.hdf5 -ip /work/long_lab/Rushani/IEDLMM/data/BPD/BDO_BARD_GRU_merge.match.tsv -o /work/long_lab/Rushani/IEDLMM/data/BPD -ik /work/long_lab/Rushani/IEDLMM/data/BPD/BPD_8_kinship.RRMwe