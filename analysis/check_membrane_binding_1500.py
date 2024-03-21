"""Script that calculates the violations in Nups binding to the membrane.
Violations are defined by comparing the spring constant to the expected standard deviation.
Requires an RMF file (sys.argv[1]), the spring constant (sys.argv[2], in kJ/molA**2) the time (sys.argv[3]), which is used to select the membrane parameters form the selection below. Assumes a temperature of 1500 K."""
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

# dictionaries for membrane diameter and height
membrane_D={'5min':515.3,'6min':583.9,'8min':727.4,'10min':845.9,'15min':798.3,'mature':870.0}
membrane_H={'5min':427.2,'6min':424.9,'8min':429.9,'10min':405.4,'15min':342.5,'mature':300.0}

def membrane_distance(particle, mem_D, mem_H):
    """Calculates the distance to the membrane for a given particle
    """
    # Radius is 1/2 of the diameter
    mem_R=mem_D/2

    # determine position of particle
    Rxy = np.sqrt(particle[0]*particle[0] + particle[1]*particle[1])
    Rz = particle[2]

    # Determine distance from the toroid center to the origin
    Toroid_R=mem_R+mem_H/2

    # if particle is outside pore region, compute distance in Z-dimension from the membrane
    if (Rxy > (Toroid_R)):
        # if particle is below the membrane
        if Rz<0:
            dist = -Rz - (mem_H / 2.0)
        # if particle is above the membrane
        else:
            dist = Rz - (mem_H / 2.0)
        return dist

    else:
        # compute distance in the XY-plane between the point and the toroidal center
        XY_dist=Toroid_R-Rxy
        # compute the distance in Z to the toroidal center
        # 1) if particle is below the membrane
        if Rz<0:
            Z_dist = -particle[2]
        # 2) if particle is above the membrane
        else:
            Z_dist = particle[2]
        # compute the overall distance, subtract out the minor radius of the torus
        dist=np.sqrt(XY_dist**2 + Z_dist**2) - (mem_H / 2.0)
        return dist


# Function that gets the distance between particles and a membrane
def calc_distances(dict1,mem_D,mem_H,include_radius=True):
    Nups=list(dict1.keys())
    dist_list=[]
    for i in range(0,len(Nups)):
        Nup_i=dict1[Nups[i]]
        for p_i in Nup_i:
            if include_radius:
                dist_temp=membrane_distance(p_i.get_coordinates(),mem_D,mem_H)-p_i.get_radius()
            else:
                dist_temp=membrane_distance(p_i.get_coordinates(),mem_D,mem_H)
            dist_list.append(dist_temp)
    dist_list=np.asarray(dist_list)
    return dist_list

# extracts coordinates from a hierarchy (npc), file (rmf_fh), and frameid (frameid)
# Returns a dictionary where keys are each subcomplex and each entry is a list of XYZR particles belonging to that subcomplex
def extract_binding_coords(npc,rmf_fh,frameid):
    # Load each frame and extract coordinates of Nups
    print("Extracting trajectory coordinates...")
    coords={}
    IMP.rmf.load_frame(rmf_fh, RMF.FrameID(frameid))
    for child in npc.get_children():
        particles=[]
        for nup in child.get_children():
            if "Nup155" in nup.get_name():
                for frg in nup.get_children():
                    if bool(set(IMP.atom.Fragment(frg).get_residue_indexes()).intersection(range(262, 271))):
                        particles.append(IMP.core.XYZR(frg))
            if "Nup160" in nup.get_name():
                for frg in nup.get_children():
                    if bool(set(IMP.atom.Fragment(frg).get_residue_indexes()).intersection(range(260, 268))):
                        particles.append(IMP.core.XYZR(frg))
        if len(particles)>0:
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
        for nup in child.get_children():
            if "Nup155" in nup.get_name():
                for frg in nup.get_children():
                    if bool(set(IMP.atom.Fragment(frg).get_residue_indexes()).intersection(range(262, 271))):
                        # count each XYZR particle
                        count = count + 1
            if "Nup160" in nup.get_name():
                for frg in nup.get_children():
                    if bool(set(IMP.atom.Fragment(frg).get_residue_indexes()).intersection(range(260, 268))):
                        # count each XYZR particle
                        count=count+1
    print("Done.")
    return count

# Function to check EV restraints by looping over all frames in the RMF, extracting the coordinates at each frame, and calculating the distance between residues.
# Takes in an rmf (sim_file), and the strength of the harmonic restraint (k_EV) in kcal/mol Ang^2
# Returns number of violations and the number of violations per particle that binds the membrane
def check_membrane_EV(sim_file,k_EV,mem_D,mem_H):
    # Read in RMF file
    rmf_fh = RMF.open_rmf_file_read_only(sim_file)
    model = IMP.Model()
    npc = IMP.rmf.create_hierarchies(rmf_fh, model)[0]
    Nleafs=num_leafs(npc,rmf_fh)
    print('Number of atoms: '+str(Nleafs))
    last_frame_id = rmf_fh.get_number_of_frames()
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
    k_EV = k_EV / kT
    sigma_EV = np.sqrt(1 / k_EV)
    print('sigma: ' + str(sigma_EV))
    viol_list=[]
    norm_viol_list=[]
    # Loop over all frames. Get pairwise distances and compare to the expected standard deviation
    #for i in range(last_frame_id):
    for i in range(last_frame_id):
        viol=0
        coords=extract_binding_coords(npc,rmf_fh,i)
        distances=calc_distances(coords,mem_D,mem_H)
        print(np.max(distances))
        for dist in distances:
            if dist>sigma_EV*3:
                viol=viol+1
                print(dist)
        viol_list.append(viol)
        norm_viol_list.append(viol/Nleafs)
    return viol_list,norm_viol_list


k=float(sys.argv[3])
time=sys.argv[2]
EV_viol,EV_viol_norm=check_membrane_EV(sys.argv[1],k,membrane_D[time],membrane_H[time])
EV_viol=np.asarray(EV_viol)
EV_viol_norm=np.asarray(EV_viol_norm)
np.savetxt('bind_viol_membrane_1500.txt',EV_viol)
np.savetxt('bind_viol_norm_membrane_1500.txt',EV_viol_norm)
