"""
Script that starts one round of initial simulations for each snapshot model.
The round of simulation is in the results_dir, first characters of the seed, and job script.
Paths for main_dir and code_dir should be adjusted for the location of the outer directory.
Simulations are launched by copying a template directory (results_dir + template), which contains main.py
Job scripts to run simulations are automatically written and submitted to the queue.
"""

import sys
import os
import numpy as np

# Directories
main_dir='/.../simulations/'
code_dir='/.../code/'
sim_num='200'

results_dir=main_dir+'results'+sim_num

# parameters
# take this data from the included Nups directory. Includes times of each simulation and number of states for each
times=['5min','6min','8min','10min','15min','mature']
# top 3 combos:
# N_sims={'5min':14,'6min':16,'8min':16,'10min':14,'15min':9,'mature':1}
# top 4 combos:
N_sims={'5min':15,'6min':20,'8min':22,'10min':18,'15min':15,'mature':1}
# top 5 combos:
# N_sims={'5min':19,'6min':24,'8min':26,'10min':24,'15min':16,'mature':1}
# random seed to ensure simulations are different
seed=4*sim_num+'52945'
# number of mpi slots to use for replica exchange
mpi_np=str(8)
# parameters for membrane model
membrane_R={'5min':'515.3','6min':'583.9','8min':'727.4','10min':'845.9','15min':'798.3','mature':'870.0'}
membrane_H={'5min':'427.2','6min':'424.9','8min':'429.9','10min':'405.4','15min':'342.5','mature':'300.0'}




def write_config(filename,main_dir,code_dir,sim_index,time,membrane_R,membrane_H):
    new=open(filename,'w')

    # write model components
    new.write('# toggles for components of score function\nmodel_components:\n')
    new.write('        go_model: True  # go-like terms\n')
    new.write('        truncated: False  # use truncated harmonic in go-model (This version is not implemented)\n')
    new.write('        useBonds: False   # toggle for using Bonds implementation\n')
    new.write('        useTruncated: False # use truncated harmonic in go-model\n')
    new.write('        useGaussian: False # use Gaussian potential for go-model\n')
    new.write('        sterics: True  # include protein-protein sterics\n')
    new.write('        gmm_restraint: True  # include density restraints\n')
    new.write('        membrane: True  # include membrane\n')
    new.write('        OverallPos: True # Restraint on the position of every Nup relative to its position in the fully assembled complex\n')
    new.write('        XYRadial: False # Restraint on the postion of the nr y-complex to be the same distance from the X-axis as in the assembled complex\n')
    new.write('        XYAxial: False # Restraint on the postion of the nr y-complex to the same X and Y positions as in the assembled complex\n')
    new.write('        Zscore: False # restraint the upper and lower Z-position of all leaves\n')
    new.write('        randomize_start: True   # randomize initial configuration\n')
    new.write('        randomize_inside: False  # when true, Nups are only placed on the nucleoplasm side\n')
    new.write('        brownian_dynamics: False  # toggle to sample with BD, if false, uses MC\n')
    new.write('        steepest_descent: False  # perform steepest descent before sampling\n\n')

    # paramters for nup representation
    new.write('# parameters for nup representation\nnup_representation:\n')
    new.write('        cg_model_rmf: \"'+code_dir+'data/cg_models/fitted/'+time+'_fitted.rmf\"\n')
    new.write('        remove_nups: True\n')
    new.write('        included_nups_file: \"'+main_dir+'included_nups/'+sim_index+'_'+time+'.config\"\n\n')

    # boundary parameters
    new.write('# boundary parameters\n')
    new.write('boundary:\n')
    new.write('        boundary_k: 10.0\n')
    new.write('        boundary_zo: 0.0\n')
    new.write('        L: 2000.0  # simulation box edge length\n\n')

    # parameters for membrane model
    new.write('# parameters for membrane model\n')
    new.write('membrane:\n')
    new.write('        mem_type: "toroid"\n')
    new.write('        R: '+membrane_R[time]+' # inner diameter for toroidal pore model, major raduis is (R+h)/2, added 100A to account for membrane thickness\n')
    new.write('        h: '+membrane_H[time]+' # slab thickness in z for toroidal pore model\n')
    new.write('        k_mem: 0.001\n')
    new.write('        zo_mem: 0.0\n')
    new.write('        zo_upper: 37.0\n')
    new.write('        zo_lower: 0.0\n')
    new.write('        ev_mem: 0.01\n\n')

    # Go ff parameters
    new.write('# Go ff parameters\n')
    new.write('go_model:\n')
    new.write('        go_k: 0.01\n')
    new.write('        go_cutoff: 50.0\n')
    new.write('        sterics_k: 0.01\n')
    new.write('        sterics_slack: 100.0\n')
    new.write('        trunc_cut: 10\n')
    new.write('        weight: 1.0  # can set this to the string "recip" to set weight to 1 / num_go_bonds\n\n')

    # em score parameters
    new.write('# em score parameters\n')
    new.write('em_score:\n')
    new.write('        gmm_target_density: \"'+code_dir+'data/fit_etdata_with_gmm/Andrew_run/data/'+time+'_150.gmm\"\n')
    new.write('        gmm_sig_mass: 1.0\n')
    new.write('        gmm_scaling_factor: 10\n')
    new.write('        gmm_model_cutoff: 250.0\n')
    new.write('        gmm_data_cutoff: 250.0\n')
    new.write('        gmm_slope: 0.0\n')
    new.write('        omp_parallel: False\n')
    new.write('        gmm_weight: 1000\n\n')

    # position score
    new.write('# position score parameters\n')
    new.write('pos_score:\n')
    new.write('        pos_sigma: 1000000\n')
    new.write('        pos_tol: 0.0\n\n')

    # simulation parameters
    new.write('# simulation parameters\n')
    new.write('mc:\n')
    new.write('        steps: 1000000\n')
    new.write('        wfrq: 1000\n')
    new.write('        temperature_min: 300\n')
    new.write('        temperature_max: 450\n')
    new.write('        mc_max_trans: 1\n')
    new.write('        mc_max_rot: 0.01\n')
    new.write('        output_filename: \"'+sim_index+'_'+time+'\"\n\n')
    new.write('steepest_descent_parameters:\n')
    new.write('        sd_nsteps: 1000  # number of steepest descent steps\n\n')
    new.write('brownian_dynamics_parameters:\n')
    new.write('        time_step: 50.0  # in fs\n')
    new.write('        temperature: 300.0\n')
    new.close()

    return

def write_config_5min(filename,main_dir,code_dir,sim_index,time,membrane_R,membrane_H):
    new=open(filename,'w')

    # write model components
    new.write('# toggles for components of score function\nmodel_components:\n')
    new.write('        go_model: True  # go-like terms\n')
    new.write('        truncated: False  # use truncated harmonic in go-model (This version is not implemented)\n')
    new.write('        useBonds: False   # toggle for using Bonds implementation\n')
    new.write('        useTruncated: False # use truncated harmonic in go-model\n')
    new.write('        useGaussian: False # use Gaussian potential for go-model\n')
    new.write('        sterics: True  # include protein-protein sterics\n')
    new.write('        gmm_restraint: True  # include density restraints\n')
    new.write('        membrane: True  # include membrane\n')
    new.write('        OverallPos: True # Restraint on the position of every Nup relative to its position in the fully assembled complex\n')
    new.write('        XYRadial: False # Restraint on the postion of the nr y-complex to be the same distance from the X-axis as in the assembled complex\n')
    new.write('        XYAxial: False # Restraint on the postion of the nr y-complex to the same X and Y positions as in the assembled complex\n')
    new.write('        Zscore: False # restraint the upper and lower Z-position of all leaves\n')
    new.write('        randomize_start: True   # randomize initial configuration\n')
    new.write('        randomize_inside: False  # when true, Nups are only placed on the nucleoplasm side\n')
    new.write('        brownian_dynamics: False  # toggle to sample with BD, if false, uses MC\n')
    new.write('        steepest_descent: False  # perform steepest descent before sampling\n\n')

    # paramters for nup representation
    new.write('# parameters for nup representation\nnup_representation:\n')
    new.write('        cg_model_rmf: \"'+code_dir+'data/cg_models/fitted/'+time+'_fitted_v2.rmf\"\n')
    new.write('        remove_nups: True\n')
    new.write('        included_nups_file: \"'+main_dir+'included_nups/'+sim_index+'_'+time+'.config\"\n\n')

    # boundary parameters
    new.write('# boundary parameters\n')
    new.write('boundary:\n')
    new.write('        boundary_k: 10.0\n')
    new.write('        boundary_zo: 0.0\n')
    new.write('        L: 2000.0  # simulation box edge length\n\n')

    # parameters for membrane model
    new.write('# parameters for membrane model\n')
    new.write('membrane:\n')
    new.write('        mem_type: "toroid"\n')
    new.write('        R: '+membrane_R[time]+' # inner diameter for toroidal pore model, major raduis is (R+h)/2, added 100A to account for membrane thickness\n')
    new.write('        h: '+membrane_H[time]+' # slab thickness in z for toroidal pore model\n')
    new.write('        k_mem: 0.001\n')
    new.write('        zo_mem: 0.0\n')
    new.write('        zo_upper: 37.0\n')
    new.write('        zo_lower: 0.0\n')
    new.write('        ev_mem: 0.01\n\n')

    # Go ff parameters
    new.write('# Go ff parameters\n')
    new.write('go_model:\n')
    new.write('        go_k: 0.01\n')
    new.write('        go_cutoff: 50.0\n')
    new.write('        sterics_k: 0.01\n')
    new.write('        sterics_slack: 100.0\n')
    new.write('        trunc_cut: 10\n')
    new.write('        weight: 1.0  # can set this to the string "recip" to set weight to 1 / num_go_bonds\n\n')

    # em score parameters
    new.write('# em score parameters\n')
    new.write('em_score:\n')
    new.write('        gmm_target_density: \"'+code_dir+'data/fit_etdata_with_gmm/Andrew_run/data/'+time+'_150.gmm\"\n')
    new.write('        gmm_sig_mass: 1.0\n')
    new.write('        gmm_scaling_factor: 10\n')
    new.write('        gmm_model_cutoff: 250.0\n')
    new.write('        gmm_data_cutoff: 250.0\n')
    new.write('        gmm_slope: 0.0\n')
    new.write('        omp_parallel: False\n')
    new.write('        gmm_weight: 100\n\n')

    # position score
    new.write('# position score parameters\n')
    new.write('pos_score:\n')
    new.write('        pos_sigma: 10000000\n')
    new.write('        pos_tol: 0.0\n\n')

    # simulation parameters
    new.write('# simulation parameters\n')
    new.write('mc:\n')
    new.write('        steps: 1000000\n')
    new.write('        wfrq: 1000\n')
    new.write('        temperature_min: 300\n')
    new.write('        temperature_max: 450\n')
    new.write('        mc_max_trans: 1\n')
    new.write('        mc_max_rot: 0.01\n')
    new.write('        output_filename: \"'+sim_index+'_'+time+'\"\n\n')
    new.write('steepest_descent_parameters:\n')
    new.write('        sd_nsteps: 1000  # number of steepest descent steps\n\n')
    new.write('brownian_dynamics_parameters:\n')
    new.write('        time_step: 50.0  # in fs\n')
    new.write('        temperature: 300.0\n')
    new.close()

    return

# write job script as a single job
def write_job(filename,config_file,seed,sim_num,nslots):
    new=open(filename,'w')
    new.write('#!/bin/bash\n#$ -S /bin/bash\n#$ -cwd\n#$ -r n\n#$ -j y\n#$ -o log\n#$ -N p'+sim_num+'m'+config_file[:-4]+'\n#$ -pe smp '+nslots+'\n#$ -l h_rt=168:00:00\n\n')
    new.write('# required to limit scikit learn attempt to utilize multithreading\nexport OMP_NUM_THREADS=$NSLOTS\n\n')
    new.write('module load Sali\nmodule load boost\nmodule load mpi/openmpi-x86_64\n\n')
    new.write('mpirun -np $NSLOTS python main.py '+config_file+' '+seed+'\n\n')
    new.write('date\nhostname\n\nqstat -j $JOB_ID')
    new.close()
    return

# write job script as an array
def write_array(filename,sim_time,nsims,seed,nslots,sim_num):
    new=open(filename,'w')
    new.write('#!/bin/bash\n#$ -S /bin/bash\n#$ -cwd\n#$ -r n\n#$ -j y\n#$ -N pm'+sim_num+'_'+sim_time+'\n#$ -pe smp '+nslots+'\n#$ -l h_rt=335:00:00\n#$ -t 1-'+str(nsims)+'\n\n')
    new.write('# required to limit scikit learn attempt to utilize multithreading\nexport OMP_NUM_THREADS=$NSLOTS\n\n')
    new.write('module load mpi/openmpi-x86_64\n\n')
    new.write('cd ${SGE_TASK_ID}_'+time+'\n\n')
    new.write('mpirun -np $NSLOTS python main.py '+time+'_${SGE_TASK_ID}.yml '+seed+'\n\n')
    new.write('date\nhostname\n\nqstat -j $JOB_ID')
    new.close()
    return


# loop over all times
for i in range(0,len(times)):
    time=times[i]
    N_states=N_sims[time]
    # loop over all states
    for j in range(0,N_states):
        index=j+1

        # make a directory for this simulation
        os.chdir(results_dir)
        if time=='5min':
            time_dir = results_dir + '/' + time+'_v3_s7g2'
        else:
            time_dir=results_dir+'/'+time
        if not os.path.exists(time_dir):
            os.mkdir(time_dir)
        os.chdir(time_dir)
        state_dir=time_dir+'/'+str(index)+'_'+time
        template_dir=results_dir+'/template'
        if not os.path.exists(state_dir):
            os.system('cp -r '+template_dir+' '+state_dir)
        os.chdir(state_dir)
        # write config file
        config_file=time+'_'+str(index)+'.yml'
        if time=='5min':
            write_config_5min(filename=config_file,main_dir=main_dir,code_dir=code_dir,sim_index=str(index),time=time,membrane_R=membrane_R,membrane_H=membrane_H)
        else:
            write_config(filename=config_file,main_dir=main_dir,code_dir=code_dir,sim_index=str(index),time=time,membrane_R=membrane_R,membrane_H=membrane_H)
        # write job script
        job_file='job.sh'
        write_job(job_file,config_file,seed,sim_num,mpi_np)
        #os.system('qsub '+job_file)
    os.chdir(time_dir)
    array_file = 'array.sh'
    write_array(array_file, time, N_states, seed, mpi_np,sim_num)
    os.system('qsub ' + array_file)



