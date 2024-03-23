"""
Script that reads in the translation / rotation enumerations and writes out the mature NPC structure that best fits the ET map at each time point.
Looks for data in main_dir or main_dir2. Looks for the mature structure at model_rmf.
"""


import sys
import os
import IMP
import RMF
import IMP.rmf
import IMP.em
import IMP.npcassembly
import IMP.algebra
import numpy as np

times=['5min','6min','8min','10min','15min','mature']

model_rmf = '/.../data/cg_models/10/npc_cg.rmf'
main_dir='/.../align_NPC/enumerate/run/'
main_dir2='/.../align_NPC/enumerate/run_rot/'

nrot = 45
ntrans=101
max_trans=50


translation=np.zeros((ntrans,1))
for i in range(0,ntrans):
    stride=max_trans*2/(ntrans-1)
    translation[i]=-1*max_trans+i*stride

# cc_max is an array of maximum cross correlation values. The first column is the cross correlation. The 2nd is the translation used to achieve that cross correlation and the 3rd column is the rotation to achieve that correlation.
cc_max=np.zeros((len(times),3))
cc_max[:,0]=cc_max[:,0]-1
rot=[0,5,10,15,20,25,30,35,40]
for i in range(0,len(times)):
    for j in range(0,ntrans):
        results_file=main_dir+times[i]+"/cc_mat"+str(translation[j,0])+".txt"
        #print(os.path.exists(results_file))
        if os.path.exists(results_file):
            cc_mat=np.loadtxt(results_file)
            for k in range(0,nrot):
                cc_temp=cc_mat[k,1]
                if cc_temp>cc_max[i,0]:
                    cc_max[i,0]=cc_mat[k,1] # The second column of cc_mat is the cross-correlation
                    cc_max[i, 1]=translation[j]
                    cc_max[i, 2] = cc_mat[k,0] # The first column of cc_mat is the rotation in degrees
        else:
            for index in range(0,len(rot)):
                results_file2 = time_dir2 = main_dir2 + '/' + times[i] + "/cc_mat" + str(translation[j, 0]) + "_" + str(rot[index]) + ".txt"
                cc_mat = np.loadtxt(results_file2)
                for k in range(0, int(nrot/len(rot))):
                    cc_temp = cc_mat[k, 1]
                    if cc_temp > cc_max[i, 0]:
                        cc_max[i, 0] = cc_mat[k, 1]  # The second column of cc_mat is the cross-correlation
                        cc_max[i, 1] = translation[j]
                        cc_max[i, 2] = cc_mat[k, 0]  # The first column of cc_mat is the rotation in degrees
print(cc_max)

# loop over times. At each time read in the model_rmf and apply the translation to that causes the highest cross-correlation
for i in range(0,len(times)):
    # Determine transformation
    m=IMP.Model()
    angle=cc_max[i, 2]*3.1415926535897931/180.0 # rotation is in degrees, convert to radians
    dist=cc_max[i, 1]
    rot = IMP.algebra.get_rotation_about_axis(IMP.algebra.Vector3D(0.0, 0.0, 1.0), angle)
    transform = IMP.algebra.Transformation3D(rot, IMP.algebra.Vector3D(0.0, 0.0, dist))
    print(transform)
    trns = IMP.core.Transform(transform)
    # transform all leaves
    npcfn = RMF.open_rmf_file_read_only(model_rmf)
    npc = IMP.rmf.create_hierarchies(npcfn, m)[0]
    IMP.rmf.load_frame(npcfn, RMF.FrameID(0))
    for leaf in IMP.atom.get_leaves(npc):
        if not IMP.core.Gaussian.get_is_setup(leaf):
            trns.apply_index(m, leaf.get_particle_index())
        if IMP.core.Gaussian.get_is_setup(leaf):
            lfrm = IMP.core.Gaussian(leaf).get_reference_frame()
            newfrm = IMP.algebra.get_transformed(lfrm, transform)
            IMP.core.Gaussian(leaf).set_reference_frame(newfrm)
    # Save transformed model to new file
    outfile = RMF.create_rmf_file(times[i] + "_fitted.rmf")
    outfile.set_description("RMF fitted to ET data")
    # use a "display" hierarchy here to avoid writing gaussian particles to rmf.
    IMP.rmf.add_hierarchy(outfile, npc)
    IMP.rmf.save_frame(outfile, str(0))








