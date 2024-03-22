### This directory contains code to score spatiotemporal models from a series of snapshot models.
Code for both checking the convergence of the snapshot models and preparing the snapshot models for analysis is also included here.
For each script, we give a brief description and sample input.

#### gather_results_1model_200.py - Script to gather the lowest scoring structural state from each replica of each replica exchange simulation for each snapshot model. Gathers the first 200 independent sampling runs.
python /.../score_graph/gather_results_1model_200.py

#### gather_results_1model.py - Script to gather the lowest scoring structural state from each replica of each replica exchange simulation for each snapshot model. Only gathers states from snapshots that are not yet converged.
python /.../score_graph/gather_results_1model.py

#### prepare_filtered_noNup188.py - Script to filter scores from snapshot models. Scores less than or equal to filter_fraction are kept. Also checks which states pass the KS test.
python /.../score_graph/prepare_filtered_noNup188.py

#### score_trj_noNup188.py - Script to create a spatiotemporal model. Creates the full model. Excludes Nup188/Seh1 from the scoring function for independent validation.
python /.../score_graph/score_trj_noNup188.py

#### score_trj_notemp.py - Script to create a spatiotemporal model. Creates a model without transition scoring. Excludes Nup188/Seh1 from the scoring function for independent validation.
python /.../score_graph/score_trj_notemp.py

#### score_trj_vis.py - Script to recreate the full spatiotemporal model, and write a more visually appealing directed acyclic graph.
python /.../score_graph/score_trj_vis.py

#### score_trj_vis_notemp.py - Script to recreate the model without spatiotemporal constraints, and write a more visually appealing directed acyclic graph.
python /.../score_graph/score_trj_vis_notemp.py

#### write_rmf1.py - Script to write a single rmf with the heirarchy / geometery for all structural states in a snapshot model. Which snapshot models to write are indicated by the path number.
python /.../score_graph/write_rmf1.py 1

#### write_rmf1_fix.py - Script to rewrite the final frame of each snapshot model in case it is not read in correctly by rmf_cat initially (check # of frames in the output of write_rmf1.py, as it should be the same as the energy files).
python /.../score_graph/write_rmf1_fix.py 1

#### write_rmf2.py - Script to write one rmf for each half of samplings of a snapshot model.
python /.../score_graph/write_rmf2.py 1

#### purity.py - Calculates the model precision (purity) of the full assembly model.
python /.../score_graph/purity.py
