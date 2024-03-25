#!/bin/bash                         
#$ -S /bin/bash
#$ -cwd
#$ -r n
#$ -j y
#$ -o log
#$ -N gmm_test
#$ -l h_rt=48:00:00
#$ -l h_vmem=24G

python /wynton/home/sali/aplatham/NPC_Assembly/Tempkin_rerun/pm_assembly_code/src/tools/pdb_to_isd_gmm.py 5min_150.pdb 5min_150.gmm
python /wynton/home/sali/aplatham/NPC_Assembly/Tempkin_rerun/pm_assembly_code/src/tools/pdb_to_isd_gmm.py 6min_150.pdb 6min_150.gmm
python /wynton/home/sali/aplatham/NPC_Assembly/Tempkin_rerun/pm_assembly_code/src/tools/pdb_to_isd_gmm.py 8min_150.pdb 8min_150.gmm
python /wynton/home/sali/aplatham/NPC_Assembly/Tempkin_rerun/pm_assembly_code/src/tools/pdb_to_isd_gmm.py 10min_150.pdb 10min_150.gmm
python /wynton/home/sali/aplatham/NPC_Assembly/Tempkin_rerun/pm_assembly_code/src/tools/pdb_to_isd_gmm.py 15min_150.pdb 15min_150.gmm
python /wynton/home/sali/aplatham/NPC_Assembly/Tempkin_rerun/pm_assembly_code/src/tools/pdb_to_isd_gmm.py mature_150.pdb mature_150.gmm
python /wynton/home/sali/aplatham/NPC_Assembly/Tempkin_rerun/pm_assembly_code/src/tools/pdb_to_isd_gmm.py mature_no_channel_150.pdb mature_no_channel_150.gmm


date
hostname

qstat -j $JOB_ID
