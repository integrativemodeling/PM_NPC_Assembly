"""
Script that restarts a simulation by aligning the rmf to an existing rmf and updating the scoring function to scores_pm_refine2.
Model parameters are read in from a yml file. See the sample yml: 1_mature_step2.yml
yml files are unique to each snapshot model and are generated by start_refine.py
"""
# set random seed
import random
import sys
random.seed(int(sys.argv[2]))


# ------------------------------------------
#          IMPORTS AND HEADERS
# ------------------------------------------

import IMP
import RMF
import time
import yaml
import numpy as np
import IMP.rmf
import IMP.core
import IMP.npc
import IMP.mpi
import IMP.isd
import IMP.pmi.macros
import IMP.algebra
import IMP.container
sys.path.insert(0, "/.../tools")
import scores_pm_refine2

# path to configuration file
config_file = sys.argv[1]

# read parameters
with open(config_file, "r") as ymlfile:
    prm = yaml.load(ymlfile, Loader=yaml.FullLoader)

# handle for model components
mdlc = prm["model_components"]

# simulation parameters
simp = prm["mc"]

# initialize Replica Exchange class
rem = IMP.mpi.ReplicaExchange()
# get number of replicas
kB=0.001987204259
nproc = rem.get_number_of_replicas()
min_temp=simp["temperature_min"]*kB
max_temp=simp["temperature_max"]*kB
# create array of temperatures, in geometric progression
temp = rem.create_temperatures(min_temp, max_temp, nproc)
# get replica index
myindex = rem.get_my_index()
# set initial value of the parameter (temperature) to exchange
rem.set_my_parameter("temp", [temp[myindex]])

#IMP.set_number_of_threads(8)

# ------------------------------------------
#     INITIALIZE MODEL REPRESENTATION
# ------------------------------------------
print("Creating IMP model...")
m = IMP.Model()
IMP.set_log_level(IMP.NONE)
print(" Done.")

# representation parameters
nup_rep = prm["nup_representation"]

print("Loading input RMF from \n" + nup_rep["cg_model_rmf"]  + " ...")
# read in prepared rmf file for CG NPC structure
hc = IMP.atom.Hierarchy(m, IMP.Particle(m))
rmf_file = RMF.open_rmf_file_read_only(nup_rep["cg_model_rmf"])
hc = IMP.rmf.create_hierarchies(rmf_file, m)[0]
IMP.rmf.load_frame(rmf_file, RMF.FrameID(0))
print(" Done.")

print("Setting up nup copy number...")
# organize into subcomplex units

included_nups_file = nup_rep["included_nups_file"]
with open(included_nups_file, "r") as fh:
    included_nups = [ l.rstrip("\n") for l in fh.readlines()]

m.update()
if nup_rep["remove_nups"]:
    for sc in hc.get_children():
        if any([x in sc.get_name() for x in included_nups]):
            continue
        # bc Guassians are their own RB's, teardown before destroying
        for nup in sc.get_children():
            if "Density" in nup.get_name():
                for p in nup.get_children():
                    IMP.core.RigidBody.teardown_particle(IMP.core.RigidBody(p))
        IMP.atom.destroy(sc)

print("Done.")

print("Building boundary...")
# boundary params
bndp = prm["boundary"]
L_2 = bndp["L"] / 2  # half box length
boundary = IMP.algebra.BoundingBox3D(IMP.algebra.Vector3D(-L_2, -L_2, -L_2),
                                     IMP.algebra.Vector3D(L_2, L_2, L_2))
print(" Done.")

# ------------------------------------------
#          INITIALIZE SCORE FUNCTION
# ------------------------------------------

print("Initializing restraint terms...")
rs = scores_pm_refine2.initialize_restraints(m, hc, boundary, prm)
print("Done.")

print("Initializing score function...")
sf = IMP.core.RestraintsScoringFunction(rs)
m.update()
print(" Done.")

# translate particles to match other rmf
if mdlc["start_from_rmf"]:
    print('Determining correct transformations...')
    restart_rmf = RMF.open_rmf_file_read_only(mdlc["start_pos"])
    npc_old = IMP.rmf.create_hierarchies(restart_rmf, m)[0]
    IMP.rmf.load_frame(restart_rmf, RMF.FrameID(0))
    goal_coords={}
    starting_coords={}
    for child in npc_old.get_children():
        child_coords=[]
        if child.get_name() not in ["Data_density", "Sigma"]:
            for leaf in IMP.core.get_leaves(child):
                child_coords.append(IMP.core.XYZ(leaf).get_coordinates())
            goal_coords[child.get_name()]=child_coords
    for child in hc.get_children():
        child_coords=[]
        if child.get_name() not in ["Data_density", "Sigma"]:
            for leaf in IMP.core.get_leaves(child):
                child_coords.append(IMP.core.XYZ(leaf).get_coordinates())
            starting_coords[child.get_name()]=child_coords
    keys=goal_coords.keys()
    transform_dict={}
    for key in keys:
        from IMP.algebra import get_transformation_aligning_first_to_second
        transform_dict[key]=get_transformation_aligning_first_to_second(starting_coords[key],goal_coords[key])
    print('Done.')

# initialize rigid bodies
print("Initializing rigid bodies...")
# only make subcomplexes rigid
rigbods = [IMP.atom.create_rigid_body(mol) for mol in hc.get_children() if "spoke" in mol.get_name()]

for rb in rigbods:
    rbd = IMP.atom.RigidBodyDiffusion.setup_particle(rb)
    rbd.set_coordinates_are_optimized(True)
print("Done.")

if mdlc["start_from_rmf"]:
    print('Transforming particles to aligned positions...')
    for i in range(len(rigbods)):
        # Check that rigbods lines up with transform_list using the dictionary
        name=rigbods[i].get_name()
        name=name.replace(' rigid body','')
        IMP.core.transform(rigbods[i], transform_dict[name])
    print('Done')

# randomize initial coordinates
if mdlc["randomize_start"]:
    print("Randomizing initial configuration...")

    for rb in rigbods:

        # get a random translation and rotation for each rb
        translation = IMP.algebra.get_random_vector_in(boundary)
        rotation = IMP.algebra.get_random_rotation_3d()

        # if randomize_inside, reflect across across xy plane if above INM
        if mdlc["randomize_inside"]:
            if translation[2] > 0.0:
                translation[2] = translation[2] * -1.0

        # apply the rotation and translation to each rb
        transformation = IMP.algebra.Transformation3D(rotation, translation)
        rb.set_reference_frame(IMP.algebra.ReferenceFrame3D(transformation))

    print(" Done.")

# ------------------------------------------
#          INITIALIZE SAMPLING
# ------------------------------------------

# make boundary and membrane geometries for display
if mdlc["membrane"]:
    memp = prm["membrane"]
    if memp["mem_type"]=="toroid":
        mem_geom = IMP.npc.SlabWithToroidalPoreWireGeometry(memp["h"], (memp["R"] + memp["h"])/2.0, bndp["L"])
    elif memp["mem_type"]=="indent":
        mem_geom = IMP.npc.SlabWithSphericalIndentGeometry(memp["R"], bndp["L"], memp["h"])

# write initial coordinates as a check
myindex = rem.get_my_index()
if myindex == 0:
    print("Writing inital rmf...")
    inital_f = RMF.create_rmf_file("start_step2.rmf")
    inital_f.set_description("RMF trajectory with full score data")
    # use a "display" hierarchy here to avoid writing gaussian particles to rmf.
    IMP.rmf.add_hierarchy(inital_f, hc)
    IMP.rmf.add_restraints(inital_f, rs)
    if mdlc["membrane"]:
        IMP.rmf.add_geometry(inital_f, mem_geom)
    IMP.rmf.add_geometry(inital_f, IMP.display.BoundingBoxGeometry(boundary))
    IMP.rmf.save_frame(inital_f, str(0))
    print(" Done.")


if prm["model_components"]["steepest_descent"]:
    print("Starting steepest descent...")
    o = IMP.core.SteepestDescent(m)
    o.set_scoring_function(sf)
    o.optimize(prm["steepest_descent_parameters"]["sd_nsteps"])
    print(" Done.")


print("Opening rmf outfile...")
# write rmf with full simulation details 
outfile = RMF.create_rmf_file(simp["output_filename"] + "_" + str(rem.get_my_index()) + ".rmf")
outfile.set_description("RMF trajectory with full score data")
# use a "display" hierarchy here to avoid writing gaussian particles to rmf.
IMP.rmf.add_hierarchy(outfile, hc)
IMP.rmf.add_restraints(outfile, rs)
if mdlc["membrane"]:
    IMP.rmf.add_geometry(outfile, mem_geom)
IMP.rmf.add_geometry(outfile, IMP.display.BoundingBoxGeometry(boundary))
print(" Done.")

# set up optimizers and sampling
if prm["model_components"]["brownian_dynamics"]:
    bdprm = prm["brownian_dynamics_parameters"]
    # set optimizer states
    sos = IMP.rmf.SaveOptimizerState(m, outfile)

    o_mc = IMP.atom.BrownianDynamics(m)
    o_mc.set_scoring_function(sf)
    o_mc.set_maximum_time_step(bdprm["time_step"])  # in fs
    o_mc.set_temperature(bdprm["temperature"])

    sos.set_simulator(o_mc)
    sos.set_period(simp["wfrq"])
    o_mc.add_optimizer_state(sos)
    sos.update_always("initial configuration")

    # launch simulation
    start_time = time.time()
    o_mc.optimize(simp["steps"])
    end_time = time.time()


else:
    print("Starting MC sampling...")
    o_mc = IMP.core.MonteCarlo(m)
    o_mc.set_scoring_function(sf)
    o_mc.set_kt(temp[myindex])  # set kt in kcal/mol
    print('Temp is: '+str(temp[myindex]))
    movers = [IMP.core.RigidBodyMover(m, rb, simp["mc_max_trans"], simp["mc_max_rot"]) for rb in rigbods]
    o_mc.add_movers(movers)

    nstp = int(simp["steps"] / simp["wfrq"])
    start_time = time.time()

    # create log filename to write index record
    logname = simp["output_filename"]
    if logname.endswith(".rmf"):
        logname = logname[:-4]

    logname += "_"+str(rem.get_my_index())+".log"
    log = open(logname, "w")

    count_exchange=np.zeros((2,1))

    for step in range(nstp):
        # run wfrq number of MC steps
        score = o_mc.optimize(simp["wfrq"])

        IMP.rmf.save_frame(outfile, str(step))

        # do exchange proposal
        myindex = rem.get_my_index()
        mytemp = rem.get_my_parameter("temp")[0]
        # temperature scaled score
        myscore = score / mytemp

        log.write("%4d %2d %6.3f %6.3f\n" % (step, myindex, score, mytemp))

        findex = rem.get_friend_index(step)
        ftemp = rem.get_friend_parameter("temp", findex)[0]
        fscore = score / ftemp

        # try exchange
        exchanged = rem.do_exchange(myscore, fscore, findex)
        if (exchanged):
            o_mc.set_kt(ftemp)
            count_exchange[1]=count_exchange[1]+1
        count_exchange[0]=count_exchange[0]+1

    # get MC move statistics
    proposed = o_mc.get_number_of_proposed_steps()
    accepted = o_mc.get_number_of_accepted_steps()
    print("MC acceptance ratio "+str(myindex)+': '+str( float(accepted) / proposed))

#print('Attempted REM exchanges '+str(myindex)+': '+str(count_exchange[0]))
#print('Accepted REM exchanges '+str(myindex)+': '+str(count_exchange[1]))
print('REM exchange rate '+str(myindex)+': '+str(float(count_exchange[1]/count_exchange[0])))

print(" Done.")
end_time = time.time()
print("timing: " + str(end_time - start_time))

print("SUCCESSFUL TERMINATION")
