"""
File containing chain annotations for NPC PDB structure files.

Subcomplexes to group:

nup82
NR Y-complex
CR Y-complex
channel complex
inner pore ring commplex

Chain annotations:

5ijo:

A,B,E,K,Q,W -> nup155
C,I,O,U -> Nup93
D,P -> Nup205
F,L,R,X -> Nup p54
G,M,S,Y -> nup p58/p45
H,N,T,Z -> p62
J,V -> nup188

5ijn:

A,B,E,K,Q,W -> NUP155
C,I,O,U -> NUP93
D,J,P,V -> NUP205
F,L,R,X -> NUP54
G,M,S,Y -> NUP58
H,N,T,Z -> p62

5a9q:

Nup155 -> A,B
Nup85 -> 8,H,Q,Z
SEH1 -> 7,G,P,Y
SEC13 -> 6,F,O,X
Nup96 -> 5,E,N,W
Nup107 -> 4,D,M,V
Nup133 -> 3,C,L,U
Nup37 -> 2,K,T,b
Nup160 -> 1,J,S,a
Nup43 -> 0,9,I,R
"""

chn_5ijo = {"Chain A" : "Nup155",
            "Chain B" : "Nup155",
            "Chain C" : "Nup93_1",
            "Chain D" : "Nup205_1",
            "Chain E" : "Nup155_1",
            "Chain F" : "p54_1",
            "Chain G" : "p58_p45_1",
            "Chain H" : "p62_1",
            "Chain I" : "Nup93_2",
            "Chain J" : "Nup188_2",
            "Chain K" : "Nup155_2",
            "Chain L" : "p54_2",
            "Chain M" : "p58_p45_2",
            "Chain N" : "p62_2",
            "Chain O" : "Nup93_3",
            "Chain P" : "Nup205_3",
            "Chain Q" : "Nup155_3",
            "Chain R" : "p54_3",
            "Chain S" : "p58_p45_3",
            "Chain T" : "p62_3",
            "Chain U" : "Nup93_4",
            "Chain V" : "Nup188_4",
            "Chain W" : "Nup155_4",
            "Chain X" : "p54_4",
            "Chain Y" : "p58_p45_4",
            "Chain Z" : "p62_4"}

chn_5a9q = {"Chain 0" : "Nup43_outer_cr",
            "Chain 1" : "Nup160_outer_cr",
            "Chain 2" : "Nup37_outer_cr",
            "Chain 3" : "Nup133_inner_cr",
            "Chain 4" : "Nup107_inner_cr",
            "Chain 5" : "Nup96_inner_cr",
            "Chain 6" : "SEC13_inner_cr",
            "Chain 7" : "SEH1_inner_cr",
            "Chain 8" : "Nup85_inner_cr",
            "Chain 9" : "Nup43_inner_cr",
            "Chain A" : "Nup155_conn_2",
            "Chain B" : "Nup155_conn_1",
            "Chain C" : "Nup133_outer_nr",
            "Chain D" : "Nup107_outer_nr",
            "Chain E" : "Nup96_outer_nr",
            "Chain F" : "SEC13_outer_nr",
            "Chain G" : "SEH1_outer_nr",
            "Chain H" : "Nup85_outer_nr",
            "Chain I" : "Nup43_outer_nr",
            "Chain J" : "Nup160_outer_nr",
            "Chain K" : "Nup37_outer_nr",
            "Chain L" : "Nup133_inner_nr",
            "Chain M" : "Nup107_inner_nr",
            "Chain N" : "Nup96_inner_nr",
            "Chain O" : "SEC13_inner_nr",
            "Chain P" : "SEH1_inner_nr",
            "Chain Q" : "Nup85_inner_nr",
            "Chain R" : "Nup43_inner_nr",
            "Chain S" : "Nup160_inner_nr",
            "Chain T" : "Nup37_inner_nr",
            "Chain U" : "Nup133_outer_cr",
            "Chain V" : "Nup107_outer_cr",
            "Chain W" : "Nup96_outer_cr",
            "Chain X" : "SEC13_outer_cr",
            "Chain Y" : "SEH1_outer_cr",
            "Chain Z" : "Nup85_outer_cr",
            "Chain a" : "Nup160_inner_cr",
            "Chain b" : "Nup37_inner_cr"}


mapping = {"Nup133_outer_cr" : "yc_outer_cr",
           "Nup107_outer_cr" : "yc_outer_cr",
           "Nup96_outer_cr" : "yc_outer_cr",
           "SEC13_outer_cr" : "yc_outer_cr",
           "SEH1_outer_cr" : "yc_outer_cr",
           "Nup85_outer_cr" : "yc_outer_cr",
           "Nup43_outer_cr" : "yc_outer_cr",
           "Nup160_outer_cr" : "yc_outer_cr",
           "Nup37_outer_cr" : "yc_outer_cr",
           "Nup133_inner_cr" : "yc_inner_cr",
           "Nup107_inner_cr" : "yc_inner_cr",
           "Nup96_inner_cr" : "yc_inner_cr",
           "SEC13_inner_cr" : "yc_inner_cr",
           "SEH1_inner_cr" : "yc_inner_cr",
           "Nup85_inner_cr" : "yc_inner_cr",
           "Nup43_inner_cr" : "yc_inner_cr",
           "Nup160_inner_cr" : "yc_inner_cr",
           "Nup37_inner_cr" : "yc_inner_cr",
           "Nup133_outer_nr" : "yc_outer_nr",
           "Nup107_outer_nr" : "yc_outer_nr",
           "Nup96_outer_nr" : "yc_outer_nr",
           "SEC13_outer_nr" : "yc_outer_nr",
           "SEH1_outer_nr" : "yc_outer_nr",
           "Nup85_outer_nr" : "yc_outer_nr",
           "Nup43_outer_nr" : "yc_outer_nr",
           "Nup160_outer_nr" : "yc_outer_nr",
           "Nup37_outer_nr" : "yc_outer_nr",
           "Nup133_inner_nr" : "yc_inner_nr",
           "Nup107_inner_nr" : "yc_inner_nr",
           "Nup96_inner_nr" : "yc_inner_nr",
           "SEC13_inner_nr" : "yc_inner_nr",
           "SEH1_inner_nr" : "yc_inner_nr",
           "Nup85_inner_nr" : "yc_inner_nr",
           "Nup43_inner_nr" : "yc_inner_nr",
           "Nup160_inner_nr" : "yc_inner_nr",
           "Nup37_inner_nr" : "yc_inner_nr",
           "Nup155_conn_1" : "conn_1",
           "Nup155_conn_2" : "conn_2",
           "Nup93_1" : "ir_core_1",
           "Nup205_1" : "ir_core_1",
           "Nup155_1" : "ir_core_1",
           "p54_1" : "ir_chan_1",
           "p58_p45_1" : "ir_chan_1",
           "p62_1" : "ir_chan_1",
           "Nup93_2" : "ir_core_2",
           "Nup188_2" : "ir_core_2",
           "Nup155_2" : "ir_core_2",
           "p54_2" : "ir_chan_2",
           "p58_p45_2" : "ir_chan_2",
           "p62_2" : "ir_chan_2",
           "Nup93_3" : "ir_core_3",
           "Nup205_3" : "ir_core_3",
           "Nup155_3" : "ir_core_3",
           "p54_3" : "ir_chan_3",
           "p58_p45_3" : "ir_chan_3",
           "p62_3" : "ir_chan_3",
           "Nup93_4" : "ir_core_4",
           "Nup188_4" : "ir_core_4",
           "Nup155_4" : "ir_core_4",
           "p54_4" : "ir_chan_4",
           "p58_p45_4" : "ir_chan_4",
           "p62_4" : "ir_chan_4"}


sc_names = ["yc_inner_cr",
            "yc_outer_cr",
            "yc_inner_nr",
            "yc_outer_nr",
            "ir_core_1",
            "ir_core_2",
            "ir_core_3",
            "ir_core_4",
            "ir_chan_1",
            "ir_chan_2",
            "ir_chan_3",
            "ir_chan_4",
            "conn_1",
            "conn_2"]
