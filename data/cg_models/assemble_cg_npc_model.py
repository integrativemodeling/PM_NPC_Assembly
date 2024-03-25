"""
Script for assembling a full NPC model from atomic resolution subcomplex structures in RMF files.
"""

import IMP
import IMP.rmf
import IMP.core
import IMP.atom
import RMF
import sys
import os
from chain_annotations import mapping, sc_names
sys.path.insert(0, "/.../tools")
from pdb_to_isd_gmm import convert_pdb_to_isd
import gmm_util


# parameters
num_res = [1, 5, 10, 50, 100, 200]
#num_res = [50]
sc_groups = ["ir_channel",
             "ir_core1",
             "ir_core2",
             "Nup155",
             "yc"]
ddir = "/.../data/cg_models/templates"

sc_masses = {sc:0.0 for sc in sc_groups}

output_file_dir = "/.../data/cg_models/"

m = IMP.Model()

# read in the template library from file
subcomplexes = []
for sc in sc_groups:
    sc_file = RMF.open_rmf_file_read_only(ddir + "/" + sc + "_template.rmf")
    sc_rmf = IMP.rmf.create_hierarchies(sc_file, m)[0]
    IMP.rmf.load_frame(sc_file, RMF.FrameID(0))
    sc_rmf.set_name(sc)
    subcomplexes.append(sc_rmf)
    # compute absolute mass of template subcomplex
    sc_mass = 0.0
    for leaf in IMP.atom.get_leaves(sc_rmf):
        sc_mass += IMP.atom.Mass(leaf).get_mass()

    sc_masses[sc] = sc_mass


# now make CG templates at each resolution 
cg_reps = {}

for sc in subcomplexes:

    templates = {}
    for n in num_res:
        # build template CG models  
        cg_nups = IMP.atom.Hierarchy(m, IMP.Particle(m))
        cg_nups.set_name(sc.get_name())
        for nup in sc.get_children():
            cgn = IMP.atom.create_simplified_along_backbone(nup, n)
            cg_nups.add_child(cgn)
        templates[n] = cg_nups

    cg_reps[sc.get_name()] = templates

    # import the gmm file for the corresponding template
    # convert pdb to isd file
    pdbfn = ddir + "/" + sc.get_name() + "_template_gmm_1.pdb"
    isdfn = ddir + "/" + sc.get_name() + "_template.isd"
    convert_pdb_to_isd(pdbfn, isdfn)

    ps_model = []
    gmm_util.decorate_gmm_from_text(isdfn,
                                    ps_model,
                                    m)

    """
    for leaf in IMP.atom.get_leaves(sc):
        points.append(IMP.core.XYZR(leaf).get_coordinates())
        mass += IMP.atom.Mass(leaf).get_mass()

    score, akaike = fit_gmm_to_points(np.array(points),
                                      25,
                                      m,
                                      ps = ps_model,
                                      min_covar = 40.0,
                                      mass_multiplier=mass)
    """
    density = IMP.atom.Hierarchy(m, IMP.Particle(m))
    density.set_name("Density")

    for indx,p in enumerate(ps_model):
        p.set_name("Gaussian_" + str(indx))
        density.add_child(p)

    for n in num_res:
        templates[n].add_child(density)

# read in full NPC structure as template
npc = IMP.atom.Hierarchy(m, IMP.Particle(m))
print("Reading full NPC structure.")
# read in prepared rmf file
npc_file = RMF.open_rmf_file_read_only("/.../data/cg_models/rmf/npc_ir_yc.rmf")
npc = IMP.rmf.create_hierarchies(npc_file, m)[0]
IMP.rmf.load_frame(npc_file, RMF.FrameID(0))

print("Translating npc model.")
# center the reference protein rmf model
com = IMP.atom.CenterOfMass.setup_particle(m, npc, IMP.atom.get_leaves(npc))
translation = IMP.algebra.Transformation3D(com.get_coordinates()*-1.0)
trns = IMP.core.Transform(translation)
# also compute NPC mass
npc_mass = 0.0
for leaf in IMP.atom.get_leaves(npc):
    trns.apply_index(m, leaf.get_particle_index())
    npc_mass += IMP.atom.Mass(leaf).get_mass()

m.update()

# ------- BUILD THE CG MODELS AND ADD DENSITIES --------------
for n in num_res:

    print("Building CG with " + str(n) + " residues per bead.")

    hc = IMP.atom.Hierarchy(m, IMP.Particle(m))
    hc.set_name("npc")

    for spoke in npc.get_children():

        name = spoke.get_name()

        print("Making cg representation for " + name)

        # group input nups into subcomplexes
        ir, yc = spoke.get_children()

        scs = [IMP.atom.Hierarchy(m, IMP.Particle(m)) for i in sc_names]
        for indx,name in enumerate(sc_names):
            scs[indx].set_name(name)


        for nup in ir.get_children() + yc.get_children():
            name = mapping[nup.get_name()]
            scs[sc_names.index(name)].add_child(nup)

        for sc in scs:
            # compute correct transformation
            name = sc.get_name()
            if any(["core_1" in name, "core_3" in name]):
                tmp = "ir_core1"
                tag = "Nup205"
            if any(["core_2" in name, "core_4" in name]):
                tmp = "ir_core2"
                tag = "Nup188"
            if any(["inner" in name, "outer" in name]):
                tmp = "yc"
                tag = "Nup160"
            if "chan" in name:
                tmp = "ir_channel"
                tag = "p54"
            if "conn" in name:
                tmp = "Nup155"
                tag = "Nup155"

            template_marker = subcomplexes[sc_groups.index(tmp)]

            for nup in template_marker.get_children():
                if tag in nup.get_name():
                    source = []
                    for leaf in IMP.atom.get_leaves(nup):
                        source.append(IMP.core.XYZ(leaf).get_coordinates())
                    break

            for nup in sc.get_children():
                if tag in nup.get_name():
                    target = []
                    for leaf in IMP.atom.get_leaves(nup):
                        target.append(IMP.core.XYZ(leaf).get_coordinates())
                    break

            transformation = IMP.algebra.get_transformation_aligning_first_to_second(source, target)
            trns = IMP.core.Transform(transformation)

            # clone template hierarchy and apply transformation
            nup = IMP.atom.create_clone(cg_reps[tmp][n])
            nup.set_name(spoke.get_name() + "_" + sc.get_name())

            # apply transform
            for leaf in IMP.atom.get_leaves(nup):
                if not IMP.core.Gaussian.get_is_setup(leaf):
                    trns.apply_index(m, leaf.get_particle_index())
                if IMP.core.Gaussian.get_is_setup(leaf):
                    lfrm = IMP.core.Gaussian(leaf).get_reference_frame()
                    newfrm = IMP.algebra.get_transformed(lfrm, transformation)
                    IMP.core.Gaussian(leaf).set_reference_frame(newfrm)
                    # scale mass correcty to be proportional to weight 
                    scaling_factor = sc_masses[tmp] / npc_mass
                    IMP.atom.Mass(leaf).set_mass(IMP.atom.Mass(leaf).get_mass()*scaling_factor)

            #rb_nup = IMP.core.RigidBody.setup_particle(m, IMP.Particle(m), IMP.atom.get_leaves(nup))
            #IMP.core.transform(rb_nup, transformation)

            hc.add_child(nup)

    print("Writing full cg model")
    # check if output directory exists, create recursively if needed 
    if not os.path.isdir(output_file_dir + "/" + str(n)):
        os.makedirs(output_file_dir + "/" + str(n))

    outfile = RMF.create_rmf_file(output_file_dir + "/" + str(n) + "/npc_cg.rmf")
    IMP.rmf.add_hierarchy(outfile, hc)
    IMP.rmf.save_frame(outfile, "0")

    # write to MRC file
    model_den = [IMP.core.Gaussian(leaf) for leaf in IMP.atom.get_leaves(hc) if IMP.core.Gaussian.get_is_setup(leaf)]
    # em = IMP.em.read_map("15min_r50.mrc", IMP.em.MRCReaderWriter())
    # bb = IMP.em.get_bounding_box(em)
    # o = em.get_origin()

    # hard coded boundary 
    bb = IMP.algebra.BoundingBox3D(IMP.algebra.Vector3D(-1000, -1000, -1000),
                                   IMP.algebra.Vector3D(1000, 1000, 1000))
    o = IMP.algebra.Vector3D(-995, -995, -995)
    mdl_mrc = gmm_util.gmm2map(model_den,
                               voxel_size=10,
                               bounding_box=bb,
                               origin=o,
                               fast=True)

    IMP.em.write_map(mdl_mrc, output_file_dir + "/" + str(n) + "/npc_cg.mrc", IMP.em.MRCReaderWriter())

