"""
Script that writes a series of alignment enumerations by stating the degree of translation.
The directory where output files will be written must be specified as main_dir. It should already exist.
"""
import sys
import numpy as np
import os

# Directories
main_dir='/.../align_NPC/enumerate/run/'

# parameters
# take this data from the included Nups directory. Includes times of each simulation and number of states for each
times=['5min','6min','8min','10min','15min','mature']
# max_trans is the largest translation to make
# ntrans is the number of translations to make
ntrans=101
max_trans=50
translation=np.zeros((ntrans,1))
for i in range(0,ntrans):
    stride=max_trans*2/(ntrans-1)
    translation[i]=-1*max_trans+i*stride


def write_job(filename,time,trans):
    new=open(filename,'w')
    new.write('#!/bin/bash\n#$ -S /bin/bash\n#$ -cwd\n#$ -r n\n#$ -j y\n#$ -o log\n#$ -N pm_enum\n#$ -pe smp 1\n#$ -l h_rt=335:00:00\n\n')
    new.write('# required to limit scikit learn attempt to utilize multithreading\nexport OMP_NUM_THREADS=$NSLOTS\n\n')
    new.write('module load mpi/openmpi-x86_64\n\n')
    new.write('mpirun -np $NSLOTS python align_npc.py '+time+' '+str(trans)+'\n\n')
    new.write('date\nhostname\n\nqstat -j $JOB_ID')
    new.close()

    return

for j in range(0,len(times)):
    for i in range(0,len(translation)):
        # make a directory for this simulation
        time=times[j]
        trans=translation[i][0]
        time_dir=main_dir+'/'+time
        os.chdir(time_dir)
        # write job script
        job_file='job'+str(trans)+'.sh'
        write_job(job_file,time,trans)
        os.system('qsub '+job_file)
