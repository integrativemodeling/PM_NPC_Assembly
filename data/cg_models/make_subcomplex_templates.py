"""
Script for writing template subcomplex structures to file.
"""

import IMP
import IMP.core
import IMP.algebra
import IMP.atom
import IMP.display
import IMP.rmf
import IMP.em
import RMF
import os
import sys
sys.path.insert(0, "/.../tools")
from gmm_util import fit_gmm_to_points, gmm2map
from chain_annotations import mapping, sc_names
import numpy as np


# where to write output files
output_file_dir = "/.../data/cg_models/templates/"
# ensure outdir exists
if not os.path.isdir(output_file_dir):
    os.mkdir(output_file_dir)

m = IMP.Model()

spoke = IMP.atom.Hierarchy(m, IMP.Particle(m))

print("Reading spoke1 rmf.")
# read in prepared rmf file
spk_file = RMF.open_rmf_file_read_only("/.../data/cg_models/rmf/spokes/spoke1.rmf")
spoke = IMP.rmf.create_hierarchies(spk_file, m)[0]
IMP.rmf.load_frame(spk_file, RMF.FrameID(0))


"""
print("Translating npc model.")
# center the reference protein rmf model
com = IMP.atom.CenterOfMass.setup_particle(m, npc, IMP.atom.get_leaves(npc))
translation = IMP.algebra.Transformation3D(com.get_coordinates()*-1.0)
trns = IMP.core.Transform(translation)
for leaf in IMP.atom.get_leaves(npc):
    trns.apply_index(m, leaf.get_particle_index())

m.update()
"""

# build library of template CG models with densities
ir, yc = spoke.get_children()

subcomplexes = []

# group inner ring subcomplexes
ir_core1 = IMP.atom.Hierarchy(m, IMP.Particle(m))
ir_core1.set_name("ir_core1")
ir_core2 = IMP.atom.Hierarchy(m, IMP.Particle(m))
ir_core2.set_name("ir_core2")
ir_channel = IMP.atom.Hierarchy(m, IMP.Particle(m))
ir_channel.set_name("ir_channel")

subcomplexes.append(ir_core1)
subcomplexes.append(ir_core2)
subcomplexes.append(ir_channel)

core1 = ["Nup205_1", "Nup93_1", "Nup155_1"]
core2 = ["Nup93_2", "Nup188_2", "Nup155_2"]
channel = ["p54_1", "p58_p45_1", "p62_1"]


# group ir nups into rb subcomplexes 
for nup in ir.get_children():
    if nup.get_name() in core1:
        ir_core1.add_child(nup)

    if nup.get_name() in core2:
        ir_core2.add_child(nup)

    if nup.get_name() in channel:
        ir_channel.add_child(nup)


# group y-complex nups into y-complex

yc_1 = IMP.atom.Hierarchy(m, IMP.Particle(m))
yc_1.set_name("yc")
connector = IMP.atom.Hierarchy(m, IMP.Particle(m))
connector.set_name("Nup155")

subcomplexes.append(yc_1)
subcomplexes.append(connector)

ycomplex = ["Nup133_inner_cr",
            "Nup107_inner_cr",
            "Nup96_inner_cr",
            "SEC13_inner_cr",
            "SEH1_inner_cr",
            "Nup85_inner_cr",
            "Nup43_inner_cr",
            "Nup160_inner_cr",
            "Nup37_inner_cr"]


template_names = [sc.get_name() for sc in subcomplexes]

for nup in yc.get_children():
    if nup.get_name() in ycomplex:
        yc_1.add_child(nup)
    if nup.get_name() == "Nup155_conn_1":
        connector.add_child(nup)

# write atomic templates to file
print("Writing atomic rmf templates to file")
for sc in subcomplexes:
    outfile = RMF.create_rmf_file(output_file_dir + "/" + sc.get_name() + "_template.rmf")
    # write rmf template
    IMP.rmf.add_hierarchy(outfile, sc)
    IMP.rmf.save_frame(outfile, "0")
    # write template in PDB format 
    selc = IMP.atom.Selection(sc)
    IMP.atom.write_pdb(selc, output_file_dir + "/" + sc.get_name() + "_template.pdb")

print("Done.")

