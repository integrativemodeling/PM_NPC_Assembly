"""
Script that computes the best fit (via cross correlation coefficient) between the forward model of the mature pore and an ET map by enumerating over multiple rotations at a single translation.
This code is designed for a single time (sys.argv[1]) and translation (maxtrans, sys.argv[2])
"""


import time
import sys
import IMP
import RMF
import IMP.rmf
import IMP.em
import IMP.npcassembly
import IMP.algebra
import numpy as np
sys.path.insert(0, "/.../tools")
from gmm_util import gmm2map_parallel

def compute_rot_ary(m, mdl, tmap, nrot,ntrans,maxtrans):
    """Return angle and translation that maximizes CC between two model densities.

    Rotates the first model first assuming octagonal rotational symmetries.

    Only translates the model by a maximum of 100 Angstroms along Z
    """

    rot_ary = np.zeros((nrot,2))
    stride = 45.0 / nrot
    stride_trans = maxtrans / ntrans

    angle = stride * 3.1415926535897931/180.0 # one degree rotation coverted to radians
    dist = stride_trans

    bb = IMP.em.get_bounding_box(tmap)
    o = tmap.get_origin()
    model_den = [IMP.core.Gaussian(leaf) for leaf in IMP.atom.get_leaves(mdl) if IMP.core.Gaussian.get_is_setup(leaf)]

    # used to track cross-correlation. Not used here
    #max_cc=-1
    #indexi=-1
    #indexj=-1

    for j in range(ntrans):

        for i in range(nrot):

            # rotate model my one degree and compute CC
            if i>0:
                rot = IMP.algebra.get_rotation_about_axis(IMP.algebra.Vector3D(0.0, 0.0, 1.0), angle)
                transform = IMP.algebra.Transformation3D(rot, IMP.algebra.Vector3D(0.0, 0.0, 0.0))
                rot_ary[i,0]=rot_ary[i-1,0]+stride
            elif i==0: # if i==0, translate as well as rotate
                rot = IMP.algebra.get_rotation_about_axis(IMP.algebra.Vector3D(0.0, 0.0, 1.0), 0)
                transform = IMP.algebra.Transformation3D(rot, IMP.algebra.Vector3D(0.0, 0.0, dist))
                rot_ary[i, 0]=0

            print(transform)
            trns = IMP.core.Transform(transform)

            # apply 1 deg rotation
            for leaf in IMP.atom.get_leaves(mdl):
                if not IMP.core.Gaussian.get_is_setup(leaf):
                    trns.apply_index(m, leaf.get_particle_index())
                if IMP.core.Gaussian.get_is_setup(leaf):
                    lfrm = IMP.core.Gaussian(leaf).get_reference_frame()
                    newfrm = IMP.algebra.get_transformed(lfrm, transform)
                    IMP.core.Gaussian(leaf).set_reference_frame(newfrm)


            # rasterize gmm into density map object
            mdl_mrc = gmm2map_parallel(model_den,voxel_size=10,bounding_box=bb,origin=o,fast=True)

            cc = IMP.em.bayesem3d_get_cross_correlation_coefficient(tmap, mdl_mrc)

            rot_ary[i,1] = cc

    return rot_ary


start_time = time.time()
state=sys.argv[1]
# compute CCC heat map against masked density with no central channel signal
model_rmf = '/.../data/cg_models/10/npc_cg.rmf'
data_mrc = '/.../data/fit_etdata_with_gmm/Andrew_run/data/'+state+'_150.mrc'

nrot = 45
ntrans=1
maxtrans=float(sys.argv[2])
scaling = 10.0


m = IMP.Model()

# rvmat = np.zeros(nrot)

# read MP structure
npcfn = RMF.open_rmf_file_read_only(model_rmf)
npc = IMP.rmf.create_hierarchies(npcfn, m)[0]
IMP.rmf.load_frame(npcfn, RMF.FrameID(0))

# handle for model GMM
# model_gmm = [leaf for leaf in IMP.atom.get_leaves(npc)
#              if IMP.core.Gaussian.get_is_setup(leaf)]

# read target density GMM
tmap = IMP.em.read_map(data_mrc, IMP.em.MRCReaderWriter())

#npc = IMP.atom.create_clone(npc)

model_gmm = [IMP.core.Gaussian(leaf)
             for leaf in IMP.atom.get_leaves(npc)
             if IMP.core.Gaussian.get_is_setup(leaf)]

for leaf in model_gmm:
    leaf.set_variances(leaf.get_variances()*scaling)

# compute rotation array against target
rvmat = compute_rot_ary(m, npc, tmap, nrot, ntrans,maxtrans)

np.savetxt("cc_mat"+sys.argv[2]+".txt", rvmat)

end_time = time.time()

print("timing: " + str(end_time - start_time))

print("SUCCESSFUL TERMINATION")
