### This directory contains code for the scoring function of snapshot models.
These scripts are implement scoring functions called by their various sampling scripts.

#### gmm_util.py - Tools for handling gmms. From IMP.isd.gmm_tools. Used in the various scoring functions.

#### scores_pm.py - Implementations of scoring terms used in the structural sampling of NPC intermediates. These scores are for the initial round of structural sampling.

#### scores_pm_refine1.py - Implementations of scoring terms used in the structural sampling of NPC intermediates. These scores are for the first refinement round of structural sampling, with excluded volume implemented at the protein level.

#### scores_pm_refine2.py - Implementations of scoring terms used in the structural sampling of NPC intermediates. These scores are for the second refinement round of structural sampling, with excluded volume implemented at the lowest level of the heirarchy.