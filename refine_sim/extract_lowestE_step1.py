"""
Script that selects the lowest energy structure from a given state and save that structure as a single RMF to continue the next refinement step.
"""

# %%
# set locations for input and output data
import sys

# global imports
import os
import numpy as np
import RMF
import IMP
import IMP.pmi.macros
import IMP.rmf
import IMP.core

replicas=8

file_start=sys.argv[1]

sim_frames=200

score = 10000000000000
frame = -1
replica = -1
# for each replica
for rep in range(0, replicas):
    # Name of log file, with energies
    log =  file_start + '_step1_' + str(rep) + '.log'
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
rmf_fn=file_start+'_step1_'+str(replica)+'.rmf'
outfile=file_start+'_step1_final.rmf'
print('rmf: '+rmf_fn)
print('frame: '+str(frame))
print('score: '+str(score))
print('Writing rmf '+outfile+' ...')
# write selected rmf
model = IMP.Model()
rmf_fh = RMF.open_rmf_file_read_only(rmf_fn)
hc = IMP.rmf.create_hierarchies(rmf_fh, model)[0]
IMP.rmf.load_frame(rmf_fh, RMF.FrameID(frame))
rmf_out_fh = RMF.create_rmf_file(outfile)
IMP.rmf.add_hierarchy(rmf_out_fh, hc)
IMP.rmf.save_frame(rmf_out_fh, "0")
print('Done.')
