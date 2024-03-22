"""
Script to recreate the model without spatiotemporal constraints, and write a more visually appealing directed acyclic graph.
input_dir must be specified, and point to the filtered energies
"""
import sys
import IMP.spatiotemporal as spatiotemporal
from IMP.spatiotemporal import write_output
import os

input_dir='/wynton/group/sali/aplatham/NPC_assembly/Latham/pm_assembly/simulations_round2/Refined_energies_1model_200/filtered_noNup188/'
# go to input_dir
os.chdir(input_dir)

output_dir3=input_dir+'total_notemp_vis/'

# Input variables.
dict={'5min':19,'6min':24,'8min':24,'10min':13,'15min':8,'mature':1}
subcomplexes_original=['yc_inner_cr','yc_inner_nr','yc_outer_cr','yc_outer_nr','ir_core_1','ir_core_3','ir_core_2','ir_core_4','ir_chan_1','ir_chan_2','ir_chan_3','ir_chan_4','conn_1','conn_2']
subcomplexes=[]
for subcomplex in subcomplexes_original:
    for i in range(1,9):
        if "yc" in subcomplex:
            subcomplexes.append(str(i) + '_' + 'Nup107_' + subcomplex)
        elif "ir_core_1" in subcomplex:
            subcomplexes.append(str(i) + '_' + 'Nup205_' + subcomplex)
            subcomplexes.append(str(i) + '_' + 'Nup93_' + subcomplex)
        elif "ir_core_3" in subcomplex:
            subcomplexes.append(str(i) + '_' + 'Nup205_' + subcomplex)
            subcomplexes.append(str(i) + '_' + 'Nup93_' + subcomplex)
        elif "ir_core_2" in subcomplex:
            subcomplexes.append(str(i) + '_' + 'Nup93_' + subcomplex)
        elif "ir_core_4" in subcomplex:
            subcomplexes.append(str(i) + '_' + 'Nup93_' + subcomplex)
        elif "ir_chan" in subcomplex:
            subcomplexes.append(str(i) + '_' + 'Nup62_' + subcomplex)
        else:
            subcomplexes.append(str(i) + '_' + subcomplex)
exp_comp={'Nup107':'exp_comp_Nup107.csv','Nup205':'exp_comp_Nup205.csv','Nup62':'exp_comp_Nup62.csv','Nup93':'exp_comp_Nup93.csv'}

nodes3,graph3,graph_prob3,graph_scores3=spatiotemporal.create_DAG(dict,scorestr='_scores_tot.log',npaths=6,input_dir=input_dir,output_dir=output_dir3,spatio_temporal_rule=False,expected_subcomplexes=subcomplexes,score_comp=True,exp_comp_map=exp_comp,draw_dag=False,out_labeled_pdf=False)
write_output.draw_dag('dag_heatmap',nodes3,graph3,graph_prob3,dict.keys(),colormap='coolwarm',draw_label=False,penscale = 1.2, arrowsize = 2.4, height = '0.6', width = '0.6')
