"""
Script that checks that alignment enumerations completed successfully. Checks in 2 folders align_NPC/enumerate/run or align_NPC/enumerate/run_rot.
"""

import sys
import numpy as np
import os

# Directories
main_dir='/.../align_NPC/enumerate/run/'
main_dir2='/.../align_NPC/enumerate/run_rot/'

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

count = 0
count2 = 0
rot=[0,5,10,15,20,25,30,35,40]
for j in range(0,len(times)):
    for i in range(0,len(translation)):
        # make a directory for this simulation
        time=times[j]
        trans=translation[i][0]
        time_dir=main_dir+'/'+time
        os.chdir(time_dir)
        # write job script
        results_file=main_dir+times[j]+"/cc_mat"+str(translation[i,0])+".txt"
        # Check if results file exists
        if not os.path.exists(results_file):
            time_dir2 = main_dir2 + '/' + time
            os.chdir(time_dir2)
            count = count + 1
            check=0
            for k in range(0,len(rot)):
                results_file2=time_dir2 = main_dir2 + '/' + time+"/cc_mat" + str(translation[i,0]) + "_" + str(rot[k]) + ".txt"
                if os.path.exists(results_file2):
                    check=check+1
            if check!=len(rot):
                count2=count2+1


print('Number of initial jobs not finished: '+str(count))
print('Number of new jobs not finished: '+str(count2))
