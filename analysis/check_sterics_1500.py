"""Script that calculates the violations in Nups excluded volume.
Violations are defined by comparing the spring constant to the expected standard deviation.
Requires an RMF file (sys.argv[1]) and the spring constant (sys.argv[2], in kJ/molA**2). Assumes a temperature of 1500 K."""

import sys
import os
import math
import numpy as np
import RMF
import IMP
import IMP.pmi.macros
import IMP.rmf
import IMP.core
import IMP.mpi
from scipy import stats

# Function that calculates all pairwise distances between different Nups
def calc_distances(dict1):
    Nups=list(dict1.keys())
    dist_list=[]
    for i in range(0,len(Nups)-1):
        for j in range(i+1,len(Nups)):
            Nup_i=dict1[Nups[i]]
            Nup_j = dict1[Nups[j]]
            for p_i in Nup_i:
                for p_j in Nup_j:
                    dist=IMP.core.get_distance(p_i,p_j)
                    dist_list.append(dist)
    dist_list=np.asarray(dist_list)
    return dist_list

# extracts coordinates from a hierarchy (npc), file (rmf_fh), and frameid (frameid)
# Returns a dictionary where keys are each subcomplex and each entry is a list of XYZR particles belonging to that subcomplex
"""def extract_coords_radius(npc,rmf_fh,frameid):
    coord_list=[]

    # Load each frame and extract coordinates of Nups
    print("Extracting trajectory coordinates...")
    coords={}
    R=[]
    IMP.rmf.load_frame(rmf_fh, RMF.FrameID(frameid))
    for child in npc.get_children():
        particles=[]
        if child.get_name() not in ["Data_density", "Sigma"]:
            for leaf in IMP.core.get_leaves(child):
                sc=leaf.get_parent()
                if sc.get_name()=="Density":
                    continue
                # return each XYZR particle
                particles.append(IMP.core.XYZR(leaf))
                R.append(IMP.core.XYZR(leaf).get_radius())
            coords[child.get_name()]=particles
    print("Done.")
    R = np.asarray(R)
    return coords,R"""

# extracts coordinates from a hierarchy (npc), file (rmf_fh), and frameid (frameid)
# Returns a dictionary where keys are each subcomplex and each entry is a list of XYZR particles belonging to that subcomplex
def extract_coords(npc,rmf_fh,frameid):
    # Load each frame and extract coordinates of Nups
    print("Extracting trajectory coordinates...")
    coords={}
    IMP.rmf.load_frame(rmf_fh, RMF.FrameID(frameid))
    for child in npc.get_children():
        particles=[]
        if child.get_name() not in ["Data_density", "Sigma"]:
            for leaf in IMP.core.get_leaves(child):
                sc=leaf.get_parent()
                if sc.get_name()=="Density":
                    continue
                # return each XYZR particle
                particles.append(IMP.core.XYZR(leaf))
            coords[child.get_name()]=particles
    print("Done.")
    return coords

# counts number of particles. Assumes the leafs not named Data_density or Sigma correspond to the particles of interest
def num_leafs(npc,rmf_fh):
    # Load each frame and extract coordinates of Nups
    print("Counting atoms ...")
    count=0
    IMP.rmf.load_frame(rmf_fh, RMF.FrameID(0))
    for child in npc.get_children():
        if child.get_name() not in ["Data_density", "Sigma"]:
            for leaf in IMP.core.get_leaves(child):
                sc=leaf.get_parent()
                if sc.get_name()=="Density":
                    continue
                # count each XYZR particle
                count=count+1
    print("Done.")
    return count

# Function to check EV restraints by looping over all frames in the RMF, extracting the coordinates at each frame, and calculating the distance between residues.
# Takes in an rmf (sim_file), and the strength of the harmonic restraint (k_EV) in kcal/mol Ang^2
# Returns number of violations per
def check_EV(sim_file,k_EV):
    # Read in RMF file
    rmf_fh = RMF.open_rmf_file_read_only(sim_file)
    model = IMP.Model()
    npc = IMP.rmf.create_hierarchies(rmf_fh, model)[0]
    Nleafs=num_leafs(npc,rmf_fh)
    print('Number of atoms: '+str(Nleafs))
    last_frame_id = rmf_fh.get_number_of_frames()
    # calculate standard deviation from spring constant
    # Determine temperature to use
    rem = IMP.mpi.ReplicaExchange()
    kB=0.001987204259
    min_temp=300*kB
    max_temp=1500*kB
    # create array of temperatures, in geometric progression
    Temperature=max_temp/kB
    print(Temperature)
    # calculate standard deviation from spring constant
    kT = kB * Temperature
    # Assumes Boltzmann weighting. Needs to be adjusted for your system
    k_EV = k_EV / kT
    sigma_EV = np.sqrt(1 / k_EV)
    print('sigma: ' + str(sigma_EV))
    viol_list=[]
    norm_viol_list=[]
    # Loop over all frames. Get pairwise distances and compare to the expected standard deviation
    for i in range(last_frame_id):
        print(i)
        viol=0
        coords=extract_coords(npc,rmf_fh,i)
        distances=calc_distances(coords)
        for dist in distances:
            if dist<sigma_EV*-3:
                viol=viol+1
        viol_list.append(viol)
        norm_viol_list.append(viol/Nleafs)
    return viol_list,norm_viol_list


k=float(sys.argv[2])
EV_viol,EV_viol_norm=check_EV(sys.argv[1],k)
EV_viol=np.asarray(EV_viol)
EV_viol_norm=np.asarray(EV_viol_norm)
np.savetxt('EV_viol_1500.txt',EV_viol)
np.savetxt('EV_viol_norm_1500.txt',EV_viol_norm)
