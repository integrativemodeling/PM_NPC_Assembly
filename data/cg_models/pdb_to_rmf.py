"""
Script to generate RMF file of NPC structure by mixing 5ijo, 5a9q PDB sturctures.

This code aligns the two structures by aligning NUP155 between the two PDB files,
since this NUP is shared between the two structures.

It then writes out an RMF file containing the full cytoplasmic, inner and nuclear
rings.
"""


import IMP
import IMP.core
import IMP.algebra
import IMP.atom
import IMP.display
import IMP.rmf
import RMF
from chain_annotations import chn_5ijo,chn_5a9q
import sys
import os
import numpy as np

convert_names = True  # whether to rewrite individual nup names
outdir = "/.../data/cg_models/rmf/"
datadir = "/.../data/cg_models/pdb"

m = IMP.Model()

full_model = IMP.atom.Hierarchy(m, IMP.Particle(m))

if not os.path.isdir(outdir):
    os.mkdir(outdir)

if not os.path.isdir(outdir + "IR"):
    os.mkdir(outdir + "IR")

if not os.path.isdir(outdir + "YC"):
    os.mkdir(outdir + "YC")

if not os.path.isdir(outdir + "spokes"):
    os.mkdir(outdir + "spokes")


# process each spoke of the merged IR + YC model
for i in range(1,9):

    print("Starting Spoke " + str(i))
    # particle to organize NPC by spoke
    spoke = IMP.atom.Hierarchy(m, IMP.Particle(m))
    spoke.set_name("spoke_" + str(i))

    print("Reading IR " + str(i))
    # load IR component in this spoke
    npc_ir = IMP.atom.read_pdb(datadir + "/pdb/5ijo/0." + str(i) + ".pdb", m)
    npc_ir.set_name("ir_" + str(i))

    if convert_names:
        # rename Nup components
        for child in npc_ir.get_children():
            name = child.get_name()
            child.set_name(chn_5ijo[name])

    # prune Nup155 fragments from IR model
    child1 = npc_ir.get_children()[0]
    child2 = npc_ir.get_children()[1]

    IMP.atom.destroy(child1)
    IMP.atom.destroy(child2)

    print("Writing IR " + str(i))
    # write IR structure
    # check for IR out directory, make if needed
    outfile = RMF.create_rmf_file(outdir + "IR/ir_" + str(i) + ".rmf")
    IMP.rmf.add_hierarchy(outfile, npc_ir)
    IMP.rmf.save_frame(outfile, "0")

    print("Reading YC " + str(i))
    # load Ycomplex model
    npc_yc = IMP.atom.read_pdb(datadir + "/pdb/5a9q/0." + str(i) + ".pdb", m)
    npc_yc.set_name("yc_" + str(i))

    if convert_names:
        # convert names
        for child in npc_yc.get_children():
            name = child.get_name()
            child.set_name(chn_5a9q[name])

    print("Writing YC " + str(i))
    # write YC structures
    outfile = RMF.create_rmf_file(outdir + "YC/yc_" + str(i) + ".rmf")
    IMP.rmf.add_hierarchy(outfile, npc_yc)
    IMP.rmf.save_frame(outfile, "0")

    # merge IR + YC models
    spoke.add_child(npc_ir)
    spoke.add_child(npc_yc)

    print("Writing Spoke " + str(i))
    # write merged models
    outfile = RMF.create_rmf_file(outdir + "spokes/spoke" + str(i) + ".rmf")
    IMP.rmf.add_hierarchy(outfile, spoke)
    IMP.rmf.save_frame(outfile, "0")

    # add spoke to full model
    full_model.add_child(spoke)



# let's now center the full model
#com = IMP.atom.CenterOfMass.setup_particle(m, full_model, IMP.atom.get_leaves(full_model))
#translation = IMP.algebra.Transformation3D(com.get_coordinates()*-1.0)
#IMP.core.transformation(IMP.core.XYZ(full_model), translation)

print("Writing full model.")
# write merged models
outfile = RMF.create_rmf_file(outdir + "npc_ir_yc.rmf")
IMP.rmf.add_hierarchy(outfile, full_model)
IMP.rmf.save_frame(outfile, "0")


"""
Dividing the subcomplexes in the IR into individual modules is a key question.

It's not clear which subcomplexes should be isolated (indeed FCS data might be
required here.)

One grouping for the 5ijn structure might be:

Nup93-Nup205-Nup155 (x2)
Nup93-Nup188-Nup155 (x2)
Nup54-Nup58-Nup62 (x4)

Subcomplexes to group:

nup82
NR Y-complex
CR Y-complex
channel complex
inner pore ring commplex

Chain annotations:

5ijo:

A,B,E,K,Q,W -> nup155
C,I,O,U -> Nup93
D,P -> Nup205
F,L,R,X -> Nup p54
G,M,S,Y -> nup p58/p45
H,N,T,Z -> p62
J,V -> nup188

5ijn:

A,B,E,K,Q,W -> NUP155
C,I,O,U -> NUP93
D,J,P,V -> NUP205
F,L,R,X -> NUP54
G,M,S,Y -> NUP58
H,N,T,Z -> p62

5a9q:

Nup155 -> A,B
Nup85 -> 8,H,Q,Z
SEH1 -> 7,G,P,Y
SEC13 -> 6,F,O,X
Nup96 -> 5,E,N,W
Nup107 -> 4,D,M,V
Nup133 -> 3,C,L,U
Nup37 -> 2,K,T,b
Nup160 -> 1,J,S,a
Nup43 -> 0,9,I,R
"""
