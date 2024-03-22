### This directory contains code to run initial simulations of snapshot models.
Initial simulations are the first step of sampling. First, no excluded volume is included.
For each script, we give a brief description and sample input.

#### start_sim.py - Script that starts one round of initial simulations for each snapshot model. Only time points from 6 minutes to mature are run by this script.
python /.../start_sim/start_sim.py

#### start_sim_5min.py - Script that starts one round of initial simulations for each snapshot model. Only considers snapshots at 5min.
python /.../start_sim/start_sim_5min.py

#### main.py - Script that starts a replica exchange simulation to generate a snapshot model.
mpirun -np 16 python main.py mature_1.yml 1234567

#### check_sim.py - Script that checks that replica exchange simulations ran to completion. Only time points from 6 minutes to mature are checked by this script.
python /.../start_sim/check_sim.py

#### check_sim_5min.py - Script that checks that replica exchange simulations ran to completion. Only checks simulations at 5 min.
python /.../start_sim/check_sim_5min.py

#### extract_lowestE.py - Script to select the lowest energy structure from each independent sampling of each snapshot model. This structure serves as the starting point for refinement simulations.
python /.../start_sim/extract_lowestE.py