"""
Script to write one rmf for each half of samplings of a snapshot model.
Removes formatting from the rmfs to prepare them for sampcon.
Takes in the path number (sys.argv[1]). The script looks for a file named path*.txt, where * is the path number.
Ensure directories for the refined energies (energies_dir), output model (model_dir), and refinement simulations (main_dir) are correctly specified.
"""
# set locations for input and output data
import sys
# global imports
import itertools
import os
import numpy as np
import glob
from scipy import stats
import IMP
import RMF
import IMP.rmf

# function to write the rmf for one frame. It removes all non-protein particles
def strip_rmf(rmf_fname,frame,outfile):
    # load imp model
    m = IMP.Model()
    rmf_fh = RMF.open_rmf_file_read_only(rmf_fname)
    hc = IMP.rmf.create_hierarchies(rmf_fh, m)[0]

    # load the desired frame
    IMP.rmf.load_frame(rmf_fh, RMF.FrameID(frame))
    destroy_count=0
    indeces=m.get_particle_indexes()
    # loop over all indeces. Destroy all rigid bodies
    for i in range(0,len(indeces)):
        if IMP.core.RigidBody.get_is_setup(m.get_particle(indeces[i])):
            p=m.get_particle(indeces[i])
            #print(p.get_name())
            IMP.core.RigidBody.teardown_particle(IMP.core.RigidBody(p))
            IMP.atom.destroy(p)
            destroy_count=destroy_count+1
    # loop over all indeces again. Destroy all non-protein atoms
    indeces = m.get_particle_indexes()
    for i in range(0, len(indeces)):
        if "Data_density" in m.get_particle_name(indeces[i]):
            p = m.get_particle(indeces[i])
            IMP.atom.destroy(p)
        elif "Sigma" in m.get_particle_name(indeces[i]):
            p = m.get_particle(indeces[i])
            IMP.atom.destroy(p)
        elif "bond" in m.get_particle_name(indeces[i]):
            p = m.get_particle(indeces[i])
            IMP.atom.destroy(p)
        elif "Density" in m.get_particle_name(indeces[i]):
            p = m.get_particle(indeces[i])
            IMP.atom.destroy(p)
        elif "Gaussian" in m.get_particle_name(indeces[i]):
            p = m.get_particle(indeces[i])
            IMP.atom.destroy(p)
    m.update()
    # save output file
    rmf_out_fh = RMF.create_rmf_file(outfile)
    IMP.rmf.add_hierarchy(rmf_out_fh, hc)
    IMP.rmf.save_frame(rmf_out_fh,str(0))
    return

# Inputs
results_list=['results1','results2','results3','results4','results5','results6','results7','results8','results9','results10','results11','results12','results13','results14','results15','results16','results17','results18','results19','results20','results21','results22','results23','results24','results25','results26','results27','results28','results29','results30','results31','results32','results33','results34','results35','results36','results37','results38','results39','results40','results41','results42','results43','results44','results45','results46','results47','results48','results49','results50','results51','results52','results53','results54','results55','results56','results57','results58','results59','results60','results61','results62','results63','results64','results65','results66','results67','results68','results69','results70','results71','results72','results73','results74','results75','results76','results77','results78','results79','results80','results81','results82','results83','results84','results85','results86','results87','results88','results89','results90','results91','results92','results93','results94','results95','results96','results97','results98','results99','results100','results101','results102','results103','results104','results105','results106','results107','results108','results109','results110','results111','results112','results113','results114','results115','results116','results117','results118','results119','results120','results121','results122','results123','results124','results125','results126','results127','results128','results129','results130','results131','results132','results133','results134','results135','results136','results137','results138','results139','results140','results141','results142','results143','results144','results145','results146','results147','results148','results149','results150','results151','results152','results153','results154','results155','results156','results157','results158','results159','results160','results161','results162','results163','results164','results165','results166','results167','results168','results169','results170','results171','results172','results173','results174','results175','results176','results177','results178','results179','results180','results181','results182','results183','results184','results185','results186','results187','results188','results189','results190','results191','results192','results193','results194','results195','results196','results197','results198','results199','results200']
main_dir='/.../simulations_round2/'
path_num=sys.argv[1]
energies_dir='/.../simulations_round2/Refined_energies_1model_200/'
model_dir=['filtered_noNup188/total','filtered_noNup188/sampling1','filtered_noNup188/sampling2']
Nrep=8

filter_fraction=0.5

# determine which snapshot models to analyze
prefix_list=[]
path_fn=energies_dir+model_dir[0]+'/'+'path'+path_num+'.txt'
f=open(path_fn,'r')
line=f.readline()
while line:
    prefix_list.append(line[:-1])
    line=f.readline()

# Loop over all snapshot models from the chosen pathway
for prefix in prefix_list:
    ensemble_rmf1=prefix+'_sampling1.rmf'
    ensemble_rmf2=prefix+'_sampling2.rmf'
    if os.path.exists(energies_dir + model_dir[1]+'/'+ensemble_rmf1):
        print(ensemble_rmf1+' already exists. Skipping...')
        pass
    elif os.path.exists(energies_dir + model_dir[2]+'/'+ensemble_rmf2):
        print(ensemble_rmf2+' already exists. Skipping...')
        pass
    else:
        prefix_split=prefix.split('_')
        print(prefix_split)
        time=prefix_split[1]
        if time == '5min':
            time_dir = time + '_v3_s7g2/'
        else:
            time_dir = time + '/'
        print(time_dir)
        # Load in the raw energies. Recall the cutoff energy
        energy1=np.loadtxt(energies_dir+prefix+'_energies1.txt')
        energy2=np.loadtxt(energies_dir+prefix+'_energies2.txt')
        tot_energies=np.concatenate((energy1,energy2))
        tot_energies = np.sort(tot_energies)
        N_structures = len(tot_energies)
        filtered_energies=[]
        E_cut = tot_energies[int(filter_fraction * N_structures)]
        # variables for output
        count=0
        write2_flag=0
        to_cat=''
        # get length of filtered energies
        energy1_filtered = np.loadtxt(energies_dir + 'filtered_noNup188/' + prefix + '_scores1.log')
        energy2_filtered = np.loadtxt(energies_dir + 'filtered_noNup188/' + prefix + '_scores2.log')
        L1 = len(energy1_filtered)
        L2 = len(energy2_filtered)
        # loop over all results
        for i in range(len(results_list)):
            results_dir = main_dir + results_list[i] + '/' + time_dir + prefix + '/'
            # check if results belong to sampling1 or sampling 2
            if count < L1:
                os.chdir(energies_dir + model_dir[1])
                # for each replica
                for rep in range(0, Nrep):
                    log = results_dir + prefix + '_step2_' + str(rep) + '.log'
                    rmf = results_dir + prefix + '_step2_' + str(rep) + '.rmf'
                    dat = np.loadtxt(log)
                    energy = dat[:, 2]
                    if len(energy) != 200:
                        print('Error!!! Check energy at state:')
                        print(log)
                    min_E = 100000000000000
                    sim_index = -1
                    # select the minimum energy
                    for sim_step in range(0, 200):
                        temp_E = energy[199 - sim_step]
                        if temp_E < min_E:
                            min_E = temp_E
                            sim_index = 199 - sim_step
                    # Select the minimum energy state only if it is less than the cutoff
                    if min_E <= E_cut:
                        #print(rmf)
                        #print(sim_index)
                        # write temporary file for that frame
                        out_rmf = 'temp' + str(count) + '.rmf'
                        to_cat = to_cat + out_rmf + ' '
                        strip_rmf(rmf, sim_index, out_rmf)
                        count += 1
            # make sure results are used in the model
            elif count < L1+L2:
                os.chdir(energies_dir + model_dir[2])
                # for each replica
                for rep in range(0, Nrep):
                    log = results_dir + prefix + '_step2_' + str(rep) + '.log'
                    rmf = results_dir + prefix + '_step2_' + str(rep) + '.rmf'
                    dat = np.loadtxt(log)
                    energy = dat[:, 2]
                    if len(energy) != 200:
                        print('Error!!! Check energy at state:')
                        print(log)
                    min_E = 100000000000000
                    sim_index = -1
                    # select the minimum energy
                    for sim_step in range(0, 200):
                        temp_E = energy[199 - sim_step]
                        if temp_E < min_E:
                            min_E = temp_E
                            sim_index = 199 - sim_step
                    # Select the minimum energy state only if it is less than the cutoff
                    if min_E<=E_cut:
                        # write temporary file for that frame
                        out_rmf = 'temp' + str(count) + '.rmf'
                        to_cat = to_cat + out_rmf + ' '
                        strip_rmf(rmf,sim_index,out_rmf)
                        count+=1
            else:
                pass
            # concatonate only if at the length of the 1st or 2nd list of scores
            if count==L1:
                # concatonate all frames and clean up the temporary frame rmfs
                print(count)
                print('rmf_cat --force ' + to_cat + ' '+ensemble_rmf1)
                os.system('rmf_cat --force ' + to_cat + ' '+ensemble_rmf1)
                print('rm '+str(to_cat))
                os.system('rm '+str(to_cat))
                # reset string if end of 1st (resets early in loop after 2nd list)
                to_cat=''
            # write2_flag denotes whether the 2nd rmf has been written (1) or not (0)
            elif count==L1+L2 and write2_flag==0:
                # concatonate all frames and clean up the temporary frame rmfs
                print(count)
                print('rmf_cat --force ' + to_cat + ' '+ensemble_rmf2)
                os.system('rmf_cat --force ' + to_cat + ' '+ensemble_rmf2)
                print('rm '+str(to_cat))
                os.system('rm '+str(to_cat))
                write2_flag=1
