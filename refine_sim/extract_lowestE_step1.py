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

score=1000000000000

# loop over all simulation replicas
for i in range(replicas):
    # load the scores
    scores_fn=file_start+'_step1_'+str(i)+'.log'
    score_list = np.loadtxt(scores_fn)
    N = len(score_list)
    # find the lowest score
    for j in range(N):
        if score_list[j][2] < score:
            score = score_list[j][2]
            frame = j
            replica = i
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
