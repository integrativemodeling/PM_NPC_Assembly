### This directory contains code to run refinement simulations of snapshot models.
Refinement simulations are the final 2 steps of sampling. First, excluded volume is added for each protein. Second, excluded volume is added at the smallest coarse-grained bead.
For each script, we give a brief description and sample input.

#### check_sim.py - Script that checks the length of refinement simulations by examining the log files.
python /.../refine_sim/check_sim.py

#### start_refine.py -Script that starts one round of refinement simulations for each snapshot model.
python /.../refine_sim/start_refine.py

#### main_restart1.py - Script that restarts a simulation by aligning the rmf to an existing rmf and updating the scoring function to scores_pm_refine1.
mpirun -np 8 python main_restart1.py mature_1_step1.yml 78200200200200013958

#### extract_lowestE_step1.py - Script that selects the lowest energy structure from a given state and save that structure as a single RMF to continue the next refinement step.
python /.../refine_sim/extract_lowestE_step1.py 1_mature

#### main_restart2.py - Script that restarts a simulation by aligning the rmf to an existing rmf and updating the scoring function to scores_pm_refine2.
mpirun -np 8 python main_restart2.py mature_1_step2.yml 78200200200200013958