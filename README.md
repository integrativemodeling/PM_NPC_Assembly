\brief Code for spatiotemporal modeling of the NPC Assembly Pathway

# Info
IMP code for gathering information, model representation, model scoring, model sampling and filtering, and model anysis and validation of the NPC Assembly model. For use with IMP2.20.

_Author(s)_: Andrew Latham and Jeremy Tempkin

_Maintainer_: alatham13

_License_: LGPL. This library is free software; you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.

_Publications_:
- A Latham, et al., in preparation. (2024).
- S Otsuka, et al., A quantitative map of nuclear pore assembly reveals two distinct mechanisms. Nature 613, 575–581 (2023).

The code is divided into directories. They are referenced in order of the steps necessary to create the model of NPC assembly:
1. data - ...
2. align_npc - ...
3. included_nups - ...
4. start_sim - ...
5. tools - ...
6. refine_sim - Code to run refinement simulations of snapshot models.
7. score_graph - Code to score spatiotemporal models from a series of snapshot models.
8. analysis - Code to analyze the spatiotemporal model of NPC assembly and the snapshot models along the pathway.