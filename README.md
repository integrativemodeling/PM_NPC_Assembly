\brief Code for spatiotemporal modeling of the NPC Assembly Pathway

# Info
IMP code for gathering information, model representation, model scoring, model sampling and filtering, and model anysis and validation of the NPC Assembly model. For use with IMP2.20.

_Author(s)_: Andrew Latham and Jeremy Tempkin

_Maintainer_: alatham13

_License_: [CC-BY-SA-4.0](https://creativecommons.org/licenses/by-sa/4.0/legalcode).
This work is freely available under the terms of the Creative Commons
Attribution-ShareAlike 4.0 International License.

_Publications_:
- A Latham, et al., in preparation. (2024).
- S Otsuka, et al., A quantitative map of nuclear pore assembly reveals two distinct mechanisms. Nature 613, 575–581 (2023).

The code is divided into directories. Here, they are referenced in order of the steps necessary to create the model of NPC assembly:
1. data - Input information for the NPC assembly model and code for pre-processing that information.
2. align_npc - Code for the aligning the mature NPC structure to time-dependent ET maps.
3. included_nups - Code to determine which nup combinations will be used to generate snapshot models.
4. start_sim - Code to run initial simulations of snapshot models.
5. tools - Code for the scoring function of snapshot models.
6. refine_sim - Code to run refinement simulations of snapshot models.
7. score_graph - Code to score spatiotemporal models from a series of snapshot models.
8. analysis - Code to analyze the spatiotemporal model of NPC assembly and the snapshot models along the pathway.

The Zenodo version of this respository (XXX) also includes the model, in the simulations_round2 directory
