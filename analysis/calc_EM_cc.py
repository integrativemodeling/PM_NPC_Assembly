"""Script that calculates the cross correlation between a snapshot model and the ET data used to score that snapshot model.
Requires an RMF file (sys.argv[1]), the output mrc_file (sys.argv[2]), and the time at which the simulations were performed (sys.argv[3])"""

import sys
import os
import math
import numpy as np
import RMF
import IMP
import IMP.rmf
import IMP.core
import IMP.em

# function that extracts all protein leaves
def extract_particles(npc,rmf_fh,frameid):
    # Load each frame and extract coordinates of Nups
    print("Extracting trajectory coordinates...")
    IMP.rmf.load_frame(rmf_fh, RMF.FrameID(frameid))
    particles = []
    for child in npc.get_children():
        if child.get_name() not in ["Data_density", "Sigma"]:
            for leaf in IMP.core.get_leaves(child):
                sc=leaf.get_parent()
                if sc.get_name()=="Density":
                    continue
                # return each XYZR particle
                #particles.append(IMP.core.XYZR(leaf))
                particles.append(leaf)
    print("Done.")
    return particles

# function to read in an rmf and write the averaged protein density for that rmf
def write_mrc(sim_file,mrc_file,MRCresolution=20.0,voxel=10.0):
    # Read in RMF file
    rmf_fh = RMF.open_rmf_file_read_only(sim_file)
    model = IMP.Model()
    npc = IMP.rmf.create_hierarchies(rmf_fh, model)[0]
    last_frame_id = rmf_fh.get_number_of_frames()
    count_frames=0
    for i in range(last_frame_id):
    #for i in range(2):
        ps=extract_particles(npc,rmf_fh,i)
        print('Number of leafs:')
        print(len(ps))
        # calculate density map
        dmap = IMP.em.SampledDensityMap(ps, MRCresolution, voxel)
        dmap.calcRMS()
        dmap.set_was_used(True)
        # dmap2 stores the overall density map for the system. dmap is the density map at this timestep. dmap3 is a temporary variable for combining the current dmap with dmap2
        # if new, the summed density map is the initial density map
        if count_frames==0:
            dmap2=dmap
        # otherwise, the density map
        else:
            bbox1 = IMP.em.get_bounding_box(dmap2)
            bbox2 = IMP.em.get_bounding_box(dmap)
            bbox1 += bbox2
            dmap3 = IMP.em.create_density_map(bbox1, voxel)
            dmap3.set_was_used(True)
            dmap3.add(dmap)
            dmap3.add(dmap2)
            dmap2 = dmap3
        count_frames=count_frames+1
        # write mrc file for this timestep - for debugging
        #temp_mrc='frame_'+str(i)+'.mrc'
        #IMP.em.write_map(dmap, temp_mrc, IMP.em.MRCReaderWriter())
    # normalize density and write overall density map
    if count_frames==0:
        print('Warning: RMF empty, no frames were read')
    else:
        dmap2.multiply(1. / count_frames)
        IMP.em.write_map(dmap2, mrc_file, IMP.em.MRCReaderWriter())
    return dmap2

# inputs
rmf_file=sys.argv[1]
sim_mrc_file=sys.argv[2]
time=sys.argv[3]
# experimental density files
exp_mrc1='/.../data/mask_densities/'+time+'_masked_0.mrc'
exp_mrc2='/.../data/fit_etdata_with_gmm/Andrew_run/data/'+time+'_150.mrc'

# calculate mrc
model_density=write_mrc(rmf_file,sim_mrc_file)

# calculate correlation to experimental densities
exp_density1 = IMP.em.read_map(exp_mrc1, IMP.em.MRCReaderWriter())
exp_density2 = IMP.em.read_map(exp_mrc2, IMP.em.MRCReaderWriter())
cc1 = IMP.em.get_coarse_cc_coefficient(exp_density1, model_density,0,True)
cc2 = IMP.em.get_coarse_cc_coefficient(exp_density2, model_density,0,True)

# Write output
new=open('EM_cc.txt','w')
new.write('# masked_CC\tgmm_CC\n')
new.write(str(cc1)+'\t'+str(cc2))
new.close()


