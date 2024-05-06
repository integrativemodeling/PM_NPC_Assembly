"""
Script to select the lowest energy structure from each independent sampling of each snapshot model. This structure serves as the starting point for refinement simulations.
The independent sampling runs are specified by the results_list, which should be directories located in main_dir, and structures are output to lowestE_dir
"""

# %%
# set locations for input and output data
import sys

# global imports
import os
import numpy as np
import subprocess
import glob
import pandas as pd
import RMF
import IMP
import IMP.pmi.macros
import IMP.rmf
import IMP.core

results_list=['results1','results2','results3','results4','results5','results6','results7','results8','results9','results10','results11','results12','results13','results14','results15','results16','results17','results18','results19','results20','results21','results22','results23','results24','results25','results26','results27','results28','results29','results30','results31','results32','results33','results34','results35','results36','results37','results38','results39','results40','results41','results42','results43','results44','results45','results46','results47','results48','results49','results50','results51','results52','results53','results54','results55','results56','results57','results58','results59','results60','results61','results62','results63','results64','results65','results66','results67','results68','results69','results70','results71','results72','results73','results74','results75','results76','results77','results78','results79','results80','results81','results82','results83','results84','results85','results86','results87','results88','results89','results90','results91','results92','results93','results94','results95','results96','results97','results98','results99','results100','results101','results102','results103','results104','results105','results106','results107','results108','results109','results110','results111','results112','results113','results114','results115','results116','results117','results118','results119','results120','results121','results122','results123','results124','results125','results126','results127','results128','results129','results130','results131','results132','results133','results134','results135','results136','results137','results138','results139','results140','results141','results142','results143','results144','results145','results146','results147','results148','results149','results150','results151','results152','results153','results154','results155','results156','results157','results158','results159','results160','results161','results162','results163','results164','results165','results166','results167','results168','results169','results170','results171','results172','results173','results174','results175','results176','results177','results178','results179','results180','results181','results182','results183','results184','results185','results186','results187','results188','results189','results190','results191','results192','results193','results194','results195','results196','results197','results198','results199','results200']

# make sure directory parameters formatted with trailing /
lowestE_dir='/.../simulations_round2/starting_rmfs/'

# times for simulations
times=['5min','6min','8min','10min','15min','mature']
N_sims={'5min':15,'6min':20,'8min':22,'10min':18,'15min':15,'mature':1}

# results directory
main_dir='/.../simulations/'

replicas=8
sim_frames=1000

for time in times:
    # make directory for each time point
    if time=='5min':
        time_dir=lowestE_dir+time+'_v3_s7g2/'
    else:
        time_dir = lowestE_dir + time +'/'
    if not os.path.exists(time_dir):
        os.makedirs(time_dir)
    os.chdir(time_dir)
    # for each simulation
    for result_folder in results_list:
        # for each Nup combination
        for sim_index in range(1,N_sims[time]+1):
            prefix=str(sim_index)+'_'+time
            score = 10000000000000
            frame = -1
            replica = -1
            # for each replica
            for rep in range(0, replicas):
                # Name of log file, with energies
                if time == '5min':
                    log = main_dir + result_folder + '/' + time + '_v3_s7g2/' + str(sim_index) + '_' + time + '/' + str(sim_index) + '_' + time + '_' + str(rep) + '.log'
                else:
                    log = main_dir + result_folder + '/' + time + '/' + str(sim_index) + '_' + time + '/' + str(sim_index) + '_' + time + '_' + str(rep) + '.log'
                dat = np.loadtxt(log)
                if len(dat) > 0:
                    # read in all energies, starting at 1/2 way through the simulation (Note: 3rd column in this case)
                    energy = dat[:, 2]
                    # Check energy list is the correct size
                    if len(energy) != sim_frames:
                        print('Error!!! Check energy at state:')
                        print(results_dir + '/' + log)
                    for sim_step in range(0, sim_frames):
                        temp_E = energy[(sim_frames-1) - sim_step]
                        if temp_E < score:
                            score=temp_E
                            frame=(sim_frames-1) - sim_step
                            replica=rep
            # Load the rmf file. Write the frame of interest to a new file
            if time=='5min':
                rmf_fn=main_dir+result_folder+'/'+time+'_v3_s7g2/'+str(sim_index)+'_'+time+'/'+str(sim_index)+'_'+time+'_'+str(replica)+'.rmf'
            else:
                rmf_fn=main_dir+result_folder+'/'+time+'/'+str(sim_index)+'_'+time+'/'+str(sim_index)+'_'+time+'_'+str(replica)+'.rmf'
            outfile=result_folder + '_' + str(sim_index)+'_'+time+'.rmf'
            print('rmf: '+rmf_fn)
            print('frame: '+str(frame))
            print('score: '+str(score))
            print('Writing rmf '+outfile+' ...')
            model = IMP.Model()
            rmf_fh = RMF.open_rmf_file_read_only(rmf_fn)
            hc = IMP.rmf.create_hierarchies(rmf_fh, model)[0]
            IMP.rmf.load_frame(rmf_fh, RMF.FrameID(frame))
            rmf_out_fh = RMF.create_rmf_file(time_dir + outfile)
            IMP.rmf.add_hierarchy(rmf_out_fh, hc)
            IMP.rmf.save_frame(rmf_out_fh, "0")
            print('Done.')
