### This directory contains code to run initial simulations of snapshot models.
Initial simulations are the first step of sampling. First, no excluded volume is included.
For each script, we give a brief description and sample input.

#### start_sim.py - Script that starts one round of initial simulations for each snapshot model.
python /.../start_sim/start_sim.py

#### main.py - Script that starts a replica exchange simulation to generate a snapshot model.
mpirun -np 8 python main.py mature_1.yml 1234567

#### check_sim.py - Script that checks that replica exchange simulations ran to completion.
python /.../start_sim/check_sim.py

#### extract_lowestE.py - Script to select the lowest energy structure from each independent sampling of each snapshot model. This structure serves as the starting point for refinement simulations.
python /.../start_sim/extract_lowestE.py