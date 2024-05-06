### This directory contains code for the aligning the mature NPC structure to time-dependent ET maps.
This step selects the initial NPC structure who's native interactions are aligned to the ET map at that timepoint.

#### align_npc.py - Script that computes the best fit (via cross correlation coefficient) between the forward model of the mature pore and an ET map by enumerating over multiple rotations at a single translation.
python /.../align_npc/align_npc.py 5min 0.0

#### write_align.py - Script that writes a series of alignment enumerations by stating the degree of translation.
python /.../align_npc/write_align.py

#### check_pm.py - Script that checks that alignment enumerations completed successfully. Checks in 2 folders align_NPC/enumerate/run or align_NPC/enumerate/run_rot.
python /.../align_npc/check_pm.py

#### read_scores.py - Script that reads in the translation / rotation enumerations and reads out the best cross-correlation for each translation between the density from the mature NPC structure and the ET map.
python /.../align_npc/read_scores.py

#### read_align.py - Script that reads in the translation / rotation enumerations and writes out the mature NPC structure that best fits the ET map at each time point.
python /.../align_npc/read_align.py