#!/usr/bin/env python

import unittest
import os
import sys
import IMP.spatiotemporal as spatiotemporal
from IMP.spatiotemporal import analysis
import IMP.test
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


# General paths
TOPDIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
data_path = os.path.join(TOPDIR, 'simulations_round2', 'Refined_energies_1model_460', 'filtered_noNup188')
old_pdf_path = os.path.join(TOPDIR, 'simulations_round2', 'Refined_energies_1model_460', 'filtered_noNup188', 'total')

sys.path.insert(0, TOPDIR+"/tools")
import scores_pm_refine2

def run_sim(config_file):
    """Function that restarts a simulation by aligning the rmf to an existing rmf and updating the scoring function to scores_pm_refine2"""
    # read parameters
    with open(config_file, "r") as ymlfile:
        prm = yaml.load(ymlfile, Loader=yaml.SafeLoader)

    # handle for model components
    mdlc = prm["model_components"]

    # simulation parameters
    simp = prm["mc"]

    # initialize Replica Exchange class
    rem = IMP.mpi.ReplicaExchange()
    # get number of replicas
    kB = 0.001987204259
    nproc = rem.get_number_of_replicas()
    min_temp = simp["temperature_min"] * kB
    max_temp = simp["temperature_max"] * kB
    # create array of temperatures, in geometric progression
    temp = rem.create_temperatures(min_temp, max_temp, nproc)
    # get replica index
    myindex = rem.get_my_index()
    # set initial value of the parameter (temperature) to exchange
    rem.set_my_parameter("temp", [temp[myindex]])

    # IMP.set_number_of_threads(8)

    # ------------------------------------------
    #     INITIALIZE MODEL REPRESENTATION
    # ------------------------------------------
    print("Creating IMP model...")
    m = IMP.Model()
    IMP.set_log_level(IMP.NONE)
    print(" Done.")

    # representation parameters
    nup_rep = prm["nup_representation"]

    print("Loading input RMF from \n" + nup_rep["cg_model_rmf"] + " ...")
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
        included_nups = [l.rstrip("\n") for l in fh.readlines()]

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
        goal_coords = {}
        starting_coords = {}
        for child in npc_old.get_children():
            child_coords = []
            if child.get_name() not in ["Data_density", "Sigma"]:
                for leaf in IMP.core.get_leaves(child):
                    child_coords.append(IMP.core.XYZ(leaf).get_coordinates())
                goal_coords[child.get_name()] = child_coords
        for child in hc.get_children():
            child_coords = []
            if child.get_name() not in ["Data_density", "Sigma"]:
                for leaf in IMP.core.get_leaves(child):
                    child_coords.append(IMP.core.XYZ(leaf).get_coordinates())
                starting_coords[child.get_name()] = child_coords
        keys = goal_coords.keys()
        transform_dict = {}
        for key in keys:
            from IMP.algebra import get_transformation_aligning_first_to_second
            transform_dict[key] = get_transformation_aligning_first_to_second(starting_coords[key], goal_coords[key])
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
            name = rigbods[i].get_name()
            name = name.replace(' rigid body', '')
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
        if memp["mem_type"] == "toroid":
            mem_geom = IMP.npc.SlabWithToroidalPoreWireGeometry(memp["h"], (memp["R"] + memp["h"]) / 2.0, bndp["L"])
        elif memp["mem_type"] == "indent":
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
        print('Temp is: ' + str(temp[myindex]))
        movers = [IMP.core.RigidBodyMover(m, rb, simp["mc_max_trans"], simp["mc_max_rot"]) for rb in rigbods]
        o_mc.add_movers(movers)

        nstp = int(simp["steps"] / simp["wfrq"])
        start_time = time.time()

        # create log filename to write index record
        logname = simp["output_filename"]
        if logname.endswith(".rmf"):
            logname = logname[:-4]

        logname += "_" + str(rem.get_my_index()) + ".log"
        log = open(logname, "w")

        count_exchange = np.zeros((2, 1))

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
                count_exchange[1] = count_exchange[1] + 1
            count_exchange[0] = count_exchange[0] + 1

        # get MC move statistics
        proposed = o_mc.get_number_of_proposed_steps()
        accepted = o_mc.get_number_of_accepted_steps()
        log.close()
        print("MC acceptance ratio " + str(myindex) + ': ' + str(float(accepted) / proposed))

    print('REM exchange rate ' + str(myindex) + ': ' + str(float(count_exchange[1] / count_exchange[0])))

    print(" Done.")
    end_time = time.time()
    print("timing: " + str(end_time - start_time))

    print("SUCCESSFUL TERMINATION")

def write_included_Nups(fn='1_mature.config'):
    subcomplexes=['yc_inner_cr','yc_inner_nr','yc_outer_cr','yc_outer_nr','ir_core_1','ir_core_3','ir_core_2','ir_core_4','ir_chan_1','ir_chan_2','ir_chan_3','ir_chan_4','conn_1','conn_2']
    f=open(fn, 'w')
    for subcomplex in subcomplexes:
        f.write(subcomplex + "\n")
    f.close()

def write_config_step2(filename,code_dir,steps='200000',wfrq='1000',sim_index='1',time='mature',membrane_R='870.0',membrane_H='300.0'):
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
    new.write('        randomize_start:  False  # randomize initial configuration\n')
    new.write('        randomize_inside: False  # when true, Nups are only placed on the nucleoplasm side\n')
    new.write('        start_from_rmf: True  # when true, starts from an existing rmf file\n')
    new.write('        start_pos: \"'+code_dir+'/data/cg_models/fitted/'+time+'_fitted.rmf\"\n')
    new.write('        brownian_dynamics: False  # toggle to sample with BD, if false, uses MC\n')
    new.write('        steepest_descent: False  # perform steepest descent before sampling\n\n')

    # paramters for nup representation
    new.write('# parameters for nup representation\nnup_representation:\n')
    if time == '5min':
        new.write('        cg_model_rmf: \"' + code_dir + '/data/cg_models/fitted/' + time + '_fitted_v2.rmf\"\n')
    else:
        new.write('        cg_model_rmf: \"'+code_dir+'/data/cg_models/fitted/'+time+'_fitted.rmf\"\n')
    new.write('        remove_nups: True\n')
    new.write('        included_nups_file: \"'+sim_index+'_'+time+'.config\"\n\n')

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
    new.write('        R: '+membrane_R+' # inner diameter for toroidal pore model, major raduis is (R+h)/2, added 100A to account for membrane thickness\n')
    new.write('        h: '+membrane_H+' # slab thickness in z for toroidal pore model\n')
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
    new.write('        sterics_k: 10.0\n')
    new.write('        sterics_slack: 10.0\n')
    new.write('        trunc_cut: 10\n')
    new.write('        weight: 1.0  # can set this to the string "recip" to set weight to 1 / num_go_bonds\n\n')

    # em score parameters
    new.write('# em score parameters\n')
    new.write('em_score:\n')
    new.write('        gmm_target_density: \"'+code_dir+'/data/fit_etdata_with_gmm/Andrew_run/data/'+time+'_150.gmm\"\n')
    new.write('        gmm_sig_mass: 1.0\n')
    new.write('        gmm_scaling_factor: 10\n')
    new.write('        gmm_model_cutoff: 250.0\n')
    new.write('        gmm_data_cutoff: 250.0\n')
    new.write('        gmm_slope: 0.0\n')
    new.write('        omp_parallel: False\n')
    if time=='5min':
        new.write('        gmm_weight: 100\n\n')
    else:
        new.write('        gmm_weight: 1000\n\n')

    # position score
    new.write('# position score parameters\n')
    new.write('pos_score:\n')
    if time == '5min':
        new.write('        pos_sigma: 10000000\n')
    else:
        new.write('        pos_sigma: 1000000\n')
    new.write('        pos_tol: 0.0\n\n')

    # simulation parameters
    new.write('# simulation parameters\n')
    new.write('mc:\n')
    new.write('        steps: '+steps+'\n')
    new.write('        wfrq: '+wfrq+'\n')
    new.write('        temperature_min: 300\n')
    new.write('        temperature_max: 1500\n')
    new.write('        mc_max_trans: 1\n')
    new.write('        mc_max_rot: 0.01\n')
    new.write('        output_filename: \"'+sim_index+'_'+time+'_step2\"\n\n')
    new.write('steepest_descent_parameters:\n')
    new.write('        sd_nsteps: 1000  # number of steepest descent steps\n\n')
    new.close()

    return


class TestSpatiotemporalDAG(unittest.TestCase):

    def test_modeling_script(self):
        """Test the main modeling script runs with 1_mature"""
        with IMP.test.temporary_directory() as tmpdir:
            # Full simulation length
            # steps = '200000'
            # wfrq = '1000'
            # Short simulation length for testing
            steps='10'
            wfrq='2'

            os.chdir(tmpdir)
            write_included_Nups()
            write_config_step2('mature_1.yml',TOPDIR,steps=steps,wfrq=wfrq)
            run_sim('mature_1.yml')
            rmf_path = os.path.join(tmpdir, "1_mature_step2_0.rmf")
            self.assertTrue(os.path.exists(rmf_path), "1_mature_step2_0.rmf should exist in the output directory")
            # Make sure output is logical
            rmf_fh = RMF.open_rmf_file_read_only(rmf_path)
            last_frame_id = rmf_fh.get_number_of_frames()
            expected_frames=(int(steps) / int(wfrq))
            self.assertAlmostEqual(last_frame_id, expected_frames, delta=1e-4)


if __name__ == '__main__':
    unittest.main()








