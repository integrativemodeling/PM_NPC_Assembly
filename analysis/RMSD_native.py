"""Script to calculate the RMSD between each frame in an rmf file and the native NPC structure.
Structures are not aligned, and comparisions are made to the pre-aligned rmf of the full native structure to the ET map at that time point.
Requires an rmf_file (sys.argv[1]) and the time of the snapshot model (sys.argv[2]). The native_rmf variable should direct to the folder where the included_nups are located."""
import sys
import os
import math
import numpy as np
import RMF
import IMP
import IMP.pmi.macros
import IMP.rmf
import IMP.core

time=sys.argv[2]

native_rmf='/.../data/cg_models/fitted/'+time+'_fitted.rmf'

# Function that gets the RMSD between 2 states from 2 dictionaries
def RMSD_from_dict(dict1,dict2):
    Nups=list(dict1.keys())
    coords1=[]
    coords2=[]
    for i in range(0,len(Nups)):
        for j in range(0,len(dict1[Nups[i]])):
            coords1.append(dict1[Nups[i]][j])
            coords2.append(dict2[Nups[i]][j])
    from IMP.algebra import get_rmsd
    rmsd=get_rmsd(coords1,coords2)
    return rmsd

# Function that reads in protein leaves as dictionaries with the subcomplex name as the keys and a list XYZ coordinates of all leaves as the values.
def append_coords(coord_list,rmf_fh,frameid):

    # Load each frame and extract coordinates. Give these coordinates to the clustering object
    print("Extracting trajectory coordinates...")
    coords={}
    IMP.rmf.load_frame(rmf_fh, RMF.FrameID(frameid))
    for child in npc.get_children():
        child_coords=[]
        if child.get_name() not in ["Data_density", "Sigma"]:
            for leaf in IMP.core.get_leaves(child):
                sc=leaf.get_parent()
                if sc.get_name()=="Density":
                    continue
                child_coords.append(IMP.core.XYZ(leaf).get_coordinates())
            coords[child.get_name()]=child_coords
    coord_list.append(coords)
    print("Done.")

    return coord_list

# Read in native rmf
native_coords=[]
rmf_fh = RMF.open_rmf_file_read_only(native_rmf)
model = IMP.Model()
npc = IMP.rmf.create_hierarchies(rmf_fh, model)[0]
native_coords = append_coords(native_coords, rmf_fh, 0)

# comp_coords is a list over each simulation
sim_file=sys.argv[1]
rmf_fh = RMF.open_rmf_file_read_only(sim_file)
model = IMP.Model()
npc = IMP.rmf.create_hierarchies(rmf_fh, model)[0]
last_frame_id = rmf_fh.get_number_of_frames()
coord_list=[]
for j in range(0,last_frame_id):
    print(j)
    coord_list=append_coords(coord_list,rmf_fh,j)

rmsd=[]
for j in range(0,len(coord_list)):
        rmsd.append(RMSD_from_dict(coord_list[j],native_coords[0]))
rmsd=np.asarray(rmsd)
np.savetxt('native_rmsd.txt',rmsd)


