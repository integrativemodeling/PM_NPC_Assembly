"""
Script to gather the lowest scoring structural state from each replica of each replica exchange simulation for each snapshot model. Only gathers states from snapshots that are not yet converged.
Need to set the path for the main refinement simulation directory (main_dir), the directory from the previously gathered energies (prev_energies_dirs) and the energies for the new independent sampling runs (output_dir)
Here, we demonstrate the process for 220 independent sampling runs, based on the results of 200 independent sampling runs.
"""
import sys
# global imports
import itertools
import os
import numpy as np
import glob
from scipy import stats

# inputs
results_list=['results1','results2','results3','results4','results5','results6','results7','results8','results9','results10','results11','results12','results13','results14','results15','results16','results17','results18','results19','results20','results21','results22','results23','results24','results25','results26','results27','results28','results29','results30','results31','results32','results33','results34','results35','results36','results37','results38','results39','results40','results41','results42','results43','results44','results45','results46','results47','results48','results49','results50','results51','results52','results53','results54','results55','results56','results57','results58','results59','results60','results61','results62','results63','results64','results65','results66','results67','results68','results69','results70','results71','results72','results73','results74','results75','results76','results77','results78','results79','results80','results81','results82','results83','results84','results85','results86','results87','results88','results89','results90','results91','results92','results93','results94','results95','results96','results97','results98','results99','results100','results101','results102','results103','results104','results105','results106','results107','results108','results109','results110','results111','results112','results113','results114','results115','results116','results117','results118','results119','results120','results121','results122','results123','results124','results125','results126','results127','results128','results129','results130','results131','results132','results133','results134','results135','results136','results137','results138','results139','results140','results141','results142','results143','results144','results145','results146','results147','results148','results149','results150','results151','results152','results153','results154','results155','results156','results157','results158','results159','results160','results161','results162','results163','results164','results165','results166','results167','results168','results169','results170','results171','results172','results173','results174','results175','results176','results177','results178','results179','results180','results181','results182','results183','results184','results185','results186','results187','results188','results189','results190','results191','results192','results193','results194','results195','results196','results197','results198','results199','results200','results201','results202','results203','results204','results205','results206','results207','results208','results209','results210','results211','results212','results213','results214','results215','results216','results217','results218','results219','results220']

main_dir='/.../simulations_round2/'

prev_energies_dirs=['/.../simulations_round2/Refined_energies_1model_200/']

output_dir=main_dir+'Refined_energies_1model_220/'
if not os.path.exists(output_dir):
            os.mkdir(output_dir)

times=['5min','6min','8min','10min','15min','mature']
N_sims={'5min':15,'6min':20,'8min':22,'10min':18,'15min':15,'mature':1}
# For debugging
#times=['5min','15min']
#N_sims={'5min':3,'15min':3}
Nrep=8

# loop over all times
for j in range(0,len(times)):
    time=times[j]
    if time=='5min':
        time_dir=time+'_v3_s7g2/'
    else:
        time_dir=time+'/'
    N_states = N_sims[time]
    # loop over the number of snapshots at each time point
    for k in range(0,N_states):
        count=0

        # Name of directory. Includes the state (k) and the time
        prefix=str(k+1)+'_'+time

        # Check if simulation is already converged (done=1)
        done = 0
        for prev_energies_dir in prev_energies_dirs:
            KS_test = np.loadtxt(prev_energies_dir + 'filtered_noNup188/'+ prefix + '_KS.txt')
            if KS_test[1] > 0.05:
                done = 1
                count += 1
        # if the simulation does not pass the KS test
        if done == 0:
            print(prefix)
            # blank values for output
            p_value=[]
            statistic=[]
            energy_list1=[]
            energy_list2=[]
            # Loop over all independent samplings
            for i in range(len(results_list)):
                # read in the log file and extract the energy
                results_dir=main_dir + results_list[i] +'/'+ time_dir + prefix
                print(results_dir)
                os.chdir(results_dir)
                for rep in range(0,Nrep):
                    # Name of log file, with energies
                    log=prefix+'_step2_'+str(rep)+'.log'
                    dat=np.loadtxt(log)
                    if len(dat)>0:
                        # read in all energies, starting at 1/2 way through the simulation (Note: 3rd column in this case)
                        energy = dat[:, 2]
                        # Check energy list is the correct size
                        if len(energy)!=200:
                            print('Error!!! Check energy at state:')
                            print(results_dir+'/'+log)
                        min_E=100000000000000
                        for sim_step in range(0,200):
                            temp_E=energy[199-sim_step]
                            if temp_E<min_E:
                                min_E=temp_E
                        # append the energy to the appropriate energy_list
                        if i<len(results_list)/2:
                            energy_list1.append(min_E)
                        else:
                            energy_list2.append(min_E)
            # convert each energy_list into 2 1-D arrays
            os.chdir(main_dir)
            energy_list1=np.asarray(energy_list1)
            energy_list1=energy_list1.flatten()
            energy_list2 = np.asarray(energy_list2)
            energy_list2=energy_list2.flatten()
            np.savetxt(output_dir+prefix+'_energies1.txt',energy_list1)
            np.savetxt(output_dir+prefix+'_energies2.txt',energy_list2)
            # evaluate the difference between the 2 energy lists
            result=stats.kstest(energy_list1,energy_list2)
            statistic.append(result.statistic)
            p_value.append(result.pvalue)
    
            # write new file. first column is the time 2nd column is the KS test result
            # name of output_file
            new_file = output_dir + prefix + '_KS.txt'
            new=open(new_file,'w')
            new.write('# KS-statistic\tp-value\n')
            new.write(str(statistic[0])+'\t'+str(p_value[0])+'\n')
            print('P: '+str(p_value[0])+' D: '+str(statistic[0]))
            if statistic[0]>0.3:
                print('Warning!!! Simulation '+str(prefix)+' is not converged!!!')
                print(statistic[0])
            new.close()
        else:
            # If simulations are converged already, copy files over from most recent run
            in_file1 = prev_energies_dirs[len(prev_energies_dirs)-1] + prefix+'_energies1.txt'
            out_file1=output_dir+prefix+'_energies1.txt'
            os.system('cp '+in_file1+' '+out_file1)
            in_file2 = prev_energies_dirs[len(prev_energies_dirs)-1] + prefix+'_energies2.txt'
            out_file2=output_dir+prefix+'_energies2.txt'
            os.system('cp '+in_file2+' '+out_file2)
