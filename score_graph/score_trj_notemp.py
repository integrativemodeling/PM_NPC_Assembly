"""
Script to create a spatiotemporal model. Creates a model without transition scoring. Excludes Nup188/Seh1 from the scoring function for independent validation.
input_dir must be specified, and point to the filtered energies
"""
import sys
import IMP.spatiotemporal as spatiotemporal
from IMP.spatiotemporal import analysis
import os

input_dir='/.../simulations_round2/Refined_energies_1model_200/filtered_noNup188/'
# go to input_dir
os.chdir(input_dir)

output_dir1=input_dir+'notemp_sampling1/'
output_dir2=input_dir+'notemp_sampling2/'
output_dir3=input_dir+'notemp_total/'

full_model=input_dir+'total/'

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

# create model. Create each 1/2 to evaluate the temporal precision.
os.chdir(input_dir)
nodes1,graph1,graph_prob1,graph_scores1=spatiotemporal.create_DAG(dict,scorestr='_scores1.log',npaths=6,input_dir=input_dir,output_dir=output_dir1,spatio_temporal_rule=False,expected_subcomplexes=subcomplexes,score_comp=True,exp_comp_map=exp_comp,draw_dag=False)
os.chdir(input_dir)
nodes2,graph2,graph_prob2,graph_scores2=spatiotemporal.create_DAG(dict,scorestr='_scores2.log',npaths=6,input_dir=input_dir,output_dir=output_dir2,spatio_temporal_rule=False,expected_subcomplexes=subcomplexes,score_comp=True,exp_comp_map=exp_comp,draw_dag=False)
nodes3,graph3,graph_prob3,graph_scores3=spatiotemporal.create_DAG(dict,scorestr='_scores_tot.log',npaths=6,input_dir=input_dir,output_dir=output_dir3,spatio_temporal_rule=False,expected_subcomplexes=subcomplexes,score_comp=True,exp_comp_map=exp_comp)
analysis.temporal_precision(output_dir1+'labeled_pdf.txt',output_dir2+'labeled_pdf.txt')
analysis.temporal_precision(full_model+'labeled_pdf.txt',output_dir3+'labeled_pdf.txt',output_fn='full_model_comp.txt')
analysis.purity(output_dir3+'labeled_pdf.txt',output_fn = output_dir3+'purity.txt')