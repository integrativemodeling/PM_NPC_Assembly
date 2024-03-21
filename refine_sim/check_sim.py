"""
Script that checks the length of refinement simulations by examining the log file.
"""
import sys
import os
import numpy as np

# Directories
main_dir='/.../simulations_round2/'
# list of snapshot model directories to check
results=['results1','results2','results3','results4','results5','results6','results7','results8','results9','results10','results11','results12','results13','results14','results15','results16','results17','results18','results19','results20','results21','results22','results23','results24','results25','results26','results27','results28','results29','results30','results31','results32','results33','results34','results35','results36','results37','results38','results39','results40','results41','results42','results43','results44','results45','results46','results47','results48','results49','results50','results51','results52','results53','results54','results55','results56','results57','results58','results59','results60','results61','results62','results63','results64','results65','results66','results67','results68','results69','results70','results71','results72','results73','results74','results75','results76','results77','results78','results79','results80','results81','results82','results83','results84','results85','results86','results87','results88','results89','results90','results91','results92','results93','results94','results95','results96','results97','results98','results99','results100','results101','results102','results103','results104','results105','results106','results107','results108','results109','results110','results111','results112','results113','results114','results115','results116','results117','results118','results119','results120','results121','results122','results123','results124','results125','results126','results127','results128','results129','results130','results131','results132','results133','results134','results135','results136','results137','results138','results139','results140','results141','results142','results143','results144','results145','results146','results147','results148','results149','results150','results151','results152','results153','results154','results155','results156','results157','results158','results159','results160','results161','results162','results163','results164','results165','results166','results167','results168','results169','results170','results171','results172','results173','results174','results175','results176','results177','results178','results179','results180','results181','results182','results183','results184','results185','results186','results187','results188','results189','results190','results191','results192','results193','results194','results195','results196','results197','results198','results199','results200']

# take this data from the included Nups directory. Includes times of each simulation and number of states for each
times=['5min','6min','8min','10min','15min','mature']
N_sims={'5min':19,'6min':24,'8min':24,'10min':13,'15min':8,'mature':1}



# loop over all times
done1=0
done2=0
notdone1=0
notdone2=0
for k in range(0,len(results)):
    results_dir = main_dir + results[k]
    for i in range(0,len(times)):
        time=times[i]
        N_states=N_sims[time]
        # loop over all states
        for j in range(0,N_states):
            index=j+1
            # make a directory for this simulation
            if time=='5min':
                time_dir=results_dir+'/'+time+'_v3_s7g2'
            else:
                time_dir=results_dir+'/'+time
            #os.chdir(time_dir)
            state_dir=time_dir+'/'+str(index)+'_'+time
            os.chdir(state_dir)
            log_file1 = str(index) + '_' + time + '_step1_0.log'
            log_file2 = str(index) + '_' + time + '_step2_0.log'
            if os.path.exists(log_file1):
                #print(state_dir)
                log1=np.loadtxt(log_file1)
                if len(log1)==200:
                    done1=done1+1
                else:
                    print(state_dir)
                    notdone1=notdone1+1
            else:
                print('Error, file1 does not exist:')
                print(state_dir)
                notdone1=notdone1+1
            if os.path.exists(log_file2):
                log2=np.loadtxt(log_file2)
                if len(log2)==200:
                    done2=done2+1
                else:
                    print(state_dir)
                    notdone2=notdone2+1
            else:
                print('Error, file2 does not exist:')
                print(state_dir)
                notdone2=notdone2+1
print('Done step1: '+str(done1))
print('Not done step1: '+str(notdone1))
print('Done step2: '+str(done2))
print('Not done step2: '+str(notdone2))

