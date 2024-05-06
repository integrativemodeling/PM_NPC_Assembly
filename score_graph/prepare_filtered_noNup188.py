"""
Script to filter scores from snapshot models. Scores less than or equal to filter_fraction are kept. Also checks which states pass the KS test.
Paths for main_dir (directory of the refinement simulations), directory with the FCS data (fluor_data_dir), and which gathered energies to filter (energies_dir) need to be specified.
Scales energies by kT at the geometeric mean of the sampling temperatures.
Prepares the FCS data in the format used by IMP.spatiotemporal.
"""
import sys
# global imports
import itertools
import os
import numpy as np
import glob
from scipy import stats
import pandas as pd
import IMP.mpi
rem = IMP.mpi.ReplicaExchange()

# geometeric mean of the temperatures at which the simulations were performed and Boltzmann constant
kB=0.001987204259
min_temp=300*kB
max_temp=1500*kB
# create array of temperatures, in geometric progression
boltz_weights = rem.create_temperatures(min_temp, max_temp, 8)
temp_array=boltz_weights/kB
Temperature=stats.mstats.gmean(temp_array)
print('Temperature is: '+str(Temperature))

main_dir='/.../simulations_round2/'

energies_dir=main_dir+'Refined_energies_1model_200/'
output_dir=energies_dir+'filtered_noNup188/'
if not os.path.exists(output_dir):
            os.mkdir(output_dir)

times=['5min','6min','8min','10min','15min','mature']
N_sims={'5min':15,'6min':20,'8min':22,'10min':18,'15min':15,'mature':1}

# Directories with existing data
included_nups_dir=main_dir+'included_nups/'
fluor_data_dir='/.../data/qfluor_data/'

# for debugging:
#times=['5min','6min','8min','10min','15min','mature']
#N_sims={'5min':2,'6min':2,'8min':2,'10min':2,'15min':2,'mature':1}

# Inputs for interpretting qflour data. Number corresponds to the mean copy number per pore
prots = ["Nup107", "Nup62", "Nup93", "Nup205"]
final_cn = {"Nup107": 31.3454,"Nup62": 46.4844,"Nup205": 16.4733,"Nup93": 47.2745}


# function to filter energies in input_dir, write them to output_dir. Filtering is done by selecting all structures in the top X%.
def filter_energies(main_dir,input_dir,output_dir,times,N_sims,kT=1.0,filter_fraction=0.5):

    converged1=0
    converged2=0
    count=0


    for j in range(0,len(times)):
        time=times[j]
        N_states = N_sims[time]
        for k in range(0,N_states):
            # blank values for output
            p_value=[]
            statistic=[]
            # Name of directory. Includes the state (k) and the time
            prefix=str(k+1)+'_'+time
            os.chdir(main_dir)
            energies1=np.loadtxt(input_dir+prefix+'_energies1.txt')
            energies2=np.loadtxt(input_dir+prefix+'_energies2.txt')
            # determine cutoff by concatenating both samplings and selecting the value for which a filter_fraction percent of structures score better than that value
            tot_energies=np.concatenate((energies1,energies2))
            tot_energies=np.sort(tot_energies)
            N_structures=len(tot_energies)
            E_cut=tot_energies[int(filter_fraction*N_structures)]
            # filter energies1 by the filter_fraction of both distributions
            energies1_filtered=[]
            for i in range(0,len(energies1)):
                if energies1[i] <= E_cut:
                    energies1_filtered.append(energies1[i])
            energies1_filtered=np.asarray(energies1_filtered)
            # filter energies2 by the filter_fraction of both distributions
            energies2_filtered = []
            for i in range(0, len(energies2)):
                if energies2[i] <= E_cut:
                    energies2_filtered.append(energies2[i])
            energies2_filtered = np.asarray(energies2_filtered)
            # convert energy units if necessary
            energies1_filtered=energies1_filtered/kT
            energies2_filtered=energies2_filtered/kT
            # save filtered energies to text
            np.savetxt(output_dir+prefix+'_scores1.log',energies1_filtered)
            np.savetxt(output_dir+prefix+'_scores2.log',energies2_filtered)
            np.savetxt(output_dir+prefix+'_scores_tot.log',np.concatenate((energies1_filtered,energies2_filtered)))
            # evaluate the difference between the 2 energy lists
            result=stats.kstest(energies1_filtered,energies2_filtered)
            statistic.append(result.statistic)
            p_value.append(result.pvalue)
            # write new file. first column is the time 2nd column is the KS test result
            # name of output_file
            new_file = output_dir + prefix + '_KS.txt'
            new=open(new_file,'w')
            new.write('# KS-statistic\tp-value\n')
            new.write(str(statistic[0])+'\t'+str(p_value[0])+'\n')
            if statistic[0]<0.3:
                converged1=converged1+1
            if p_value[0]>0.05:
                converged2=converged2+1
            else:
                print(prefix)
                print('P: '+str(p_value[0])+' D: '+str(statistic[0]))
            count=count+1
            new.close()
    print('Total states: '+str(count))
    print('States converged by D: '+str(converged1))
    print('States converged by p-value: '+str(converged2))
    return

# Function to copy config files to the appropriate location
def copy_config(included_nups_dir,output_dir):
    for j in range(0,len(times)):
        time=times[j]
        N_states = N_sims[time]
        for k in range(0,N_states):
            # Name of directory. Includes the state (k) and the time
            prefix=str(k+1)+'_'+time
            # set up old and new files
            config_file=included_nups_dir+prefix+'.config'
            output_file=output_dir+prefix+'.config'
            old=open(config_file,'r')
            new=open(output_file,'w')
            line=old.readline()
            while line:
                if len(line)>1:
                    for i in range(1,9):
                        if "yc" in line:
                            new.write(str(i) + '_' + 'Nup107_' + line)
                        elif "ir_core_1" in line:
                            new.write(str(i) + '_' + 'Nup205_' + line)
                            new.write(str(i) + '_' + 'Nup93_' + line)
                        elif "ir_core_3" in line:
                            new.write(str(i) + '_' + 'Nup205_' + line)
                            new.write(str(i) + '_' + 'Nup93_' + line)
                        elif "ir_core_2" in line:
                            new.write(str(i) + '_' + 'Nup93_' + line)
                        elif "ir_core_4" in line:
                            new.write(str(i) + '_' + 'Nup93_' + line)
                        elif "ir_chan" in line:
                            new.write(str(i) + '_' + 'Nup62_' + line)
                        else:
                            new.write(str(i)+'_'+line)
                line=old.readline()
            old.close()
            new.close()
    return

# Function to convert a series of time based intensities to
def convert_qflour(fluor_data_dir,nups,final_cn,times,output_dir):
    # read in available csv files
    nd = {nup: pd.read_csv(fluor_data_dir + "total_" + nup + "_homoZ.csv") for nup in nups}
    sigmas = {key: [] for key in nd.keys()}
    means = {key: [] for key in nd.keys()}
    # Read in data. Use 30 min for mature
    for t in [5, 6, 8, 10, 15, 30]:
        # for each protein
        for key in nd.keys():
            m5_data = nd[key].loc[nd[key]["Time (min)"] == t]
            for m5_key in m5_data.keys():
                for value in m5_data.loc[:, m5_key]:
                    # replace negative intensity with 0
                    if value < 0:
                        m5_data = m5_data.replace(value, 0.0)
            # Read in mean and std_dev from the data
            traces = nd[key].columns[2:]  # from the CSV, the row index is stored as well
            cn = m5_data[traces] * final_cn[key]
            means[key].append(float(cn.mean(1)))
            sigmas[key].append(float(cn.std(1)))
    # convert means and std dictionaries into pd dataframes and save to csv
    for key in nd.keys():
        nup_cn=pd.DataFrame(index=times,columns=['Time','mean','std'])
        for i in range(len(times)):
            nup_cn['Time'][times[i]]=times[i]
            nup_cn['mean'][times[i]]=means[key][i]
            nup_cn['std'][times[i]] = sigmas[key][i]
        # save to csv
        output_fn=output_dir+'exp_comp_'+key+'.csv'
        nup_cn.to_csv(output_fn)
    return

# main function calls. This script prepares the data in a format for imp_spatiotemporal --------------------------------
# Filter energies
filter_energies(main_dir,energies_dir,output_dir,times,N_sims,kT=kB*Temperature,filter_fraction=0.5)
# save config files to the new location
copy_config(included_nups_dir,output_dir)
# prepare FCS data
convert_qflour(fluor_data_dir,prots,final_cn,times,output_dir)




