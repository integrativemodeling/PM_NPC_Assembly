"""
Script that reads in the translation / rotation enumerations and reads out the best cross-correlation for each translation between the density from the mature NPC structure and the ET map.
Looks for data in main_dir or main_dir2
"""


import sys
import os
import numpy as np

times=['5min','6min','8min','10min','15min','mature']

main_dir='/.../align_NPC/enumerate/run/'
main_dir2='/.../align_NPC/enumerate/run_rot/'

nrot = 45
ntrans=101
max_trans=50


translation=np.zeros((ntrans,1+len(times)))
for i in range(0,ntrans):
    stride=max_trans*2/(ntrans-1)
    translation[i]=-1*max_trans+i*stride

# cc_max is an array of maximum cross correlation values. The first column is the cross correlation. The 2nd is the translation used to achieve that cross correlation and the 3rd column is the rotation to achieve that correlation.
rot=[0,5,10,15,20,25,30,35,40]
for i in range(0,len(times)):
    for j in range(0,ntrans):
        results_file=main_dir+times[i]+"/cc_mat"+str(translation[j,0])+".txt"
        #print(os.path.exists(results_file))
        if os.path.exists(results_file):
            cc_mat=np.loadtxt(results_file)
            translation[j,i+1]=np.max(cc_mat[:,1])
        else:
            for index in range(0,len(rot)):
                results_file2 = time_dir2 = main_dir2 + '/' + times[i] + "/cc_mat" + str(translation[j, 0]) + "_" + str(rot[index]) + ".txt"
                if index==0:
                    cc_mat = np.loadtxt(results_file2)
                elif index==len(rot)-1:
                    cc_temp=np.loadtxt(results_file2)
                    cc_mat = np.concatenate((cc_mat, cc_temp))
                    translation[j, i + 1] = np.max(cc_mat[:, 1])
                else:
                    cc_temp = np.loadtxt(results_file2)
                    cc_mat = np.concatenate((cc_mat,cc_temp))
np.savetxt('cc_vs_trans.txt',translation)