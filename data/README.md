### This directory contains the input information for the NPC assembly model and code for pre-processing that information.

#### qfluor_data - series of csv files with FCS data for each protein used in model scoring. For each protein, the intensity as a function of time is listed for independent experiments.

#### fit_etdata_with_gmm - code to fit a gmm to experimental mrc profiles. Fitting is done using gmconvert.

#### cg_models - code to convert pdbs of portions of the NPC into rmfs of a coarse grained full assembled NPC structure.
- structures with various leveling of coarse graining are shown in directions (1, 5, 10, 50, 100, and 200), where the number represents the number of residues per coarse-grained bead.
- chain_annotations.py - File containing chain annotations for NPC PDB structure files.
- pdb_to_isd_gmm.py - Script that converts quasi-PDB file from gmconvert to gmm file used by IMP.isd

#### coarse graining of stuctures was performed in 4 steps:
1. pdb_to_rmf.py - Script to generate RMF file of NPC structure by mixing 5ijo, 5a9q PDB sturctures. (python pdb_to_rmf.py)
2. make_subcomplex_templates.py - Script for writing template subcomplex structures to file. (python make_subcomplex_templates.py)
3. fit_template_gmm_gmconvert.sh - Fit gmm to each subcomplex using gmconvert.
4. assemble_cg_npc_model.py - Script for assembling a full NPC model from atomic resolution subcomplex structures in RMF files. (python assemble_cg_npc_model.py)