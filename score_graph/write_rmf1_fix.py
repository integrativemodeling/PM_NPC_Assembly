"""
Script to rewrite the final frame of each snapshot model in case it is not read in correctly by rmf_cat initially (check # of frames in the output of write_rmf1.py, as it should be the same as the energy files).
Takes in the path number (sys.argv[1]). The script looks for a file named path*.txt, where * is the path number.
Ensure directories for the refined energies (energies_dir), output model (model_dir[0]), and refinement simulations (main_dir) are correctly specified.
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

# Inputs
results_list=['results1','results2','results3','results4','results5','results6','results7','results8','results9','results10','results11','results12','results13','results14','results15','results16','results17','results18','results19','results20','results21','results22','results23','results24','results25','results26','results27','results28','results29','results30','results31','results32','results33','results34','results35','results36','results37','results38','results39','results40','results41','results42','results43','results44','results45','results46','results47','results48','results49','results50','results51','results52','results53','results54','results55','results56','results57','results58','results59','results60','results61','results62','results63','results64','results65','results66','results67','results68','results69','results70','results71','results72','results73','results74','results75','results76','results77','results78','results79','results80','results81','results82','results83','results84','results85','results86','results87','results88','results89','results90','results91','results92','results93','results94','results95','results96','results97','results98','results99','results100','results101','results102','results103','results104','results105','results106','results107','results108','results109','results110','results111','results112','results113','results114','results115','results116','results117','results118','results119','results120','results121','results122','results123','results124','results125','results126','results127','results128','results129','results130','results131','results132','results133','results134','results135','results136','results137','results138','results139','results140','results141','results142','results143','results144','results145','results146','results147','results148','results149','results150','results151','results152','results153','results154','results155','results156','results157','results158','results159','results160','results161','results162','results163','results164','results165','results166','results167','results168','results169','results170','results171','results172','results173','results174','results175','results176','results177','results178','results179','results180','results181','results182','results183','results184','results185','results186','results187','results188','results189','results190','results191','results192','results193','results194','results195','results196','results197','results198','results199','results200','results201','results202','results203','results204','results205','results206','results207','results208','results209','results210','results211','results212','results213','results214','results215','results216','results217','results218','results219','results220','results221','results222','results223','results224','results225','results226','results227','results228','results229','results230','results231','results232','results233','results234','results235','results236','results237','results238','results239','results240','results241','results242','results243','results244','results245','results246','results247','results248','results249','results250','results251','results252','results253','results254','results255','results256','results257','results258','results259','results260','results261','results262','results263','results264','results265','results266','results267','results268','results269','results270','results271','results272','results273','results274','results275','results276','results277','results278','results279','results280','results281','results282','results283','results284','results285','results286','results287','results288','results289','results290','results291','results292','results293','results294','results295','results296','results297','results298','results299','results300','results301','results302','results303','results304','results305','results306','results307','results308','results309','results310','results311','results312','results313','results314','results315','results316','results317','results318','results319','results320','results321','results322','results323','results324','results325','results326','results327','results328','results329','results330','results331','results332','results333','results334','results335','results336','results337','results338','results339','results340','results341','results342','results343','results344','results345','results346','results347','results348','results349','results350','results351','results352','results353','results354','results355','results356','results357','results358','results359','results360','results361','results362','results363','results364','results365','results366','results367','results368','results369','results370','results371','results372','results373','results374','results375','results376','results377','results378','results379','results380','results381','results382','results383','results384','results385','results386','results387','results388','results389','results390','results391','results392','results393','results394','results395','results396','results397','results398','results399','results400','results401','results402','results403','results404','results405','results406','results407','results408','results409','results410','results411','results412','results413','results414','results415','results416','results417','results418','results419','results420','results421','results422','results423','results424','results425','results426','results427','results428','results429','results430','results431','results432','results433','results434','results435','results436','results437','results438','results439','results440','results441','results442','results443','results444','results445','results446','results447','results448','results449','results450','results451','results452','results453','results454','results455','results456','results457','results458','results459','results460','results461','results462','results463','results464','results465','results466','results467','results468','results469','results470','results471','results472','results473','results474','results475','results476','results477','results478','results479','results480']
main_dir='/.../simulations_round2/'
path_num=sys.argv[1]
energies_dir='/.../simulations_round2/Refined_energies_1model_480/'
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

os.chdir(energies_dir + model_dir[0])
# Loop over all snapshot models from the chosen pathway
for prefix in prefix_list:
    ensemble_rmf=prefix+'_ensemble.rmf'
    print(ensemble_rmf)
    if os.path.exists(ensemble_rmf):
        print(ensemble_rmf+' exists. Skipping...')
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
        to_cat=''
        # get length of filtered energies
        energy1_filtered = np.loadtxt(energies_dir + 'filtered_noNup188/' + prefix + '_scores1.log')
        energy2_filtered = np.loadtxt(energies_dir + 'filtered_noNup188/' + prefix + '_scores2.log')
        L1 = len(energy1_filtered)
        L2 = len(energy2_filtered)
        # loop over all results
        for i in range(len(results_list)):
            results_dir = main_dir + results_list[i] + '/' + time_dir + prefix + '/'
            # make sure results are used in the model
            if count < L1+L2:
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
                        if count<L1+L2-2:
                            pass
                        else:
                            print(count)
                            print(rmf)
                            print(sim_index)
                            # load IMP model with hierarchy and geometry
                            model = IMP.Model()
                            rmf_fh = RMF.open_rmf_file_read_only(rmf)
                            hc = IMP.rmf.create_hierarchies(rmf_fh, model)[0]
                            geom = IMP.rmf.create_geometries(rmf_fh)
                            IMP.rmf.load_frame(rmf_fh, RMF.FrameID(sim_index))
                            # write temporary file for that frame
                            out_rmf = prefix + '_frame' + str(count) + '.rmf'
                            rmf_out_fh = RMF.create_rmf_file(out_rmf)
                            IMP.rmf.add_hierarchy(rmf_out_fh, hc)
                            IMP.rmf.add_geometries(rmf_out_fh, geom)
                            IMP.rmf.save_frame(rmf_out_fh, "0")
                        count+=1
            else:
                pass