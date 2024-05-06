### This directory contains code to analyze the spatiotemporal model of NPC assembly and the snapshot models along the pathway.
For each script, we give a brief description and sample input.

#### check_membrane_EV_1500.py - Script that calculates the violations in Nups excluded volume with the membrane.
python /.../analysis/check_membrane_EV_1500.py 2_5min_ensemble.rmf 5min 0.01

#### check_membrane_binding_1500.py - Script that calculates the violations in Nups binding to the membrane.
python /.../analysis/check_membrane_binding_1500.py 2_5min_ensemble.rmf 5min 0.001

#### check_sterics_1500.py - Script that calculates the violations in Nups excluded volume.
python /.../check_sterics_1500.py 2_5min_ensemble.rmf 10

#### prep_render.py - Script to render a more visually appealing version of our NPC model. Renders the last rmf frame. (coloring_apl.py - Python dictionary defining a standard color scheme for the Nup subcomplexes. Used by prep_render.py)
python /.../analysis/prep_render.py cluster_center_model.rmf3

#### analyze_copy_number.py - Script to analyze nup copy number of a model.
python /.../analysis/analyze_copy_number.py labeled_pdf.txt

#### calc_EM_cc.py - Script that calculates the cross correlation between a snapshot model and the ET data used to score that snapshot model.
python /.../analysis/calc_EM_cc.py 2_5min_ensemble.rmf 2_5min.mrc 5min

#### RMSD_native.py - Script to calculate the RMSD between each frame in an rmf file and the native NPC structure.
python /.../analysis/RMSD_native.py 2_5min_ensemble.rmf 5min

#### check_num_states.py - Scipt to check the total number of simulations, the number of total sampled states, and the number of good scoring states selected for a spatiotemporal model.
python /.../check_num_states.py