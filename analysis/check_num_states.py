"""
Scipt to check the total number of simulations, the number of total sampled states, and the number of good scoring states selected for a spatiotemporal model.
Path for main_dir should point to the results of the refinement simulations.
"""
import sys
# global imports
import itertools
import os
import numpy as np
import glob
from scipy import stats
import pandas as pd

main_dir='/.../simulations_round2/'

energies_dir=main_dir+'Refined_energies_1model_480/'
output_dir=energies_dir+'filtered_noNup188/'
if not os.path.exists(output_dir):
            os.mkdir(output_dir)

times=['5min','6min','8min','10min','15min','mature']
N_sims={'5min':19,'6min':24,'8min':24,'10min':13,'15min':8,'mature':1}

# function to filter energies in input_dir, write them to output_dir. Filtering is done by selecting all structures in the top X%.
def read_energies(main_dir,input_dir,output_dir,times,N_sims):

    num_sims=0
    num_configs=0


    for j in range(0,len(times)):
        time=times[j]
        N_states = N_sims[time]
        for k in range(0,N_states):
            # Name of directory. Includes the state (k) and the time
            prefix=str(k+1)+'_'+time
            os.chdir(main_dir)
            energies1=np.loadtxt(input_dir+prefix+'_energies1.txt')
            energies2=np.loadtxt(input_dir+prefix+'_energies2.txt')
            # save filtered energies to text
            energies1_filtered=np.loadtxt(output_dir+prefix+'_scores1.log')
            energies2_filtered=np.loadtxt(output_dir+prefix+'_scores2.log')
            N1=len(energies1)/8.0
            N2=len(energies2)/8.0
            N1_filtered=len(energies1_filtered)
            N2_filtered=len(energies2_filtered)
            if N1+N2<200.0:
                print('Error!!!')
                print(prefix)
            """else:
                print(prefix)
                print(N1+N2)"""
            if N1_filtered+N2_filtered < 800:
                print('Error!!!')
                print(prefix)
            else:
                print(prefix)
                print(N1_filtered+N2_filtered)
            num_sims=num_sims+N1+N2
            num_configs = num_configs + N1_filtered + N2_filtered
    print('Number of simulations: '+str(num_sims))
    print('Total number of states: '+str(num_sims*(16*1000+8*200+8*200)))
    print('Number of good scoring models: '+str(num_configs))
    return



# main function calls. This script prepares the data in a format for imp_spatiotemporal --------------------------------
# Filter energies
read_energies(main_dir,energies_dir,output_dir,times,N_sims)




