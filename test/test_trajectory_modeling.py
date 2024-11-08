#!/usr/bin/env python

import unittest
import os
import IMP.spatiotemporal as spatiotemporal
from IMP.spatiotemporal import analysis
import IMP.test

# General paths
TOPDIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
data_path = os.path.join(TOPDIR, 'simulations_round2', 'Refined_energies_1model_460', 'filtered_noNup188')
old_pdf_path = os.path.join(TOPDIR, 'simulations_round2', 'Refined_energies_1model_460', 'filtered_noNup188', 'total')



# Paths for expected and generated files
expected_pdf_path = os.path.join(old_pdf_path, "labeled_pdf.txt")
# parameters.
state_dict={'5min':15,'6min':20,'8min':22,'10min':18,'15min':15,'mature':1}
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

class TestSpatiotemporalDAG(unittest.TestCase):

    def test_create_dag_and_check_pdf(self):
        """Test if spatiotemporal.create_DAG creates a model correctly. Tests that the model matches the previous model, that the precision matches the previous precision, and that the temporal precision matches the previous temporal precision."""

        # create output directory
        with IMP.test.temporary_directory() as tmpdir:
            output = os.path.join(tmpdir, 'output/')

            # Call the function to generate output
            nodes, graph, graph_prob, graph_scores = spatiotemporal.create_DAG(
                state_dict,
                input_dir=data_path,
                scorestr='_scores_tot.log',
                output_dir=output,
                spatio_temporal_rule=True,
                expected_subcomplexes=subcomplexes,
                score_comp=True,
                exp_comp_map=exp_comp,
                draw_dag=False # there is no need for heatmap
            )
            # Check if the labeled_pdf.txt file was created
            generated_pdf_path = os.path.join(output, "labeled_pdf.txt")
            self.assertTrue(os.path.exists(generated_pdf_path), "labeled_pdf.txt should exist in the output directory")

            # Use temporal_precision to check between old and new model
            analysis.temporal_precision(expected_pdf_path,generated_pdf_path,output_fn=output+'trj_test.txt')
            f=open(output+'trj_test.txt','r')
            # First line is description
            f.readline()
            # 2nd line is the temporal precision
            trj_comp=float(f.readline())
            f.close()
            self.assertAlmostEqual(trj_comp, 1.0, delta=1e-5)
            # Test the precision of the trajecotry
            analysis.precision(generated_pdf_path,output_fn=output+'precision.txt')
            f = open(output + 'precision.txt', 'r')
            # First line is description
            f.readline()
            # 2nd line is the precision
            precision = float(f.readline())
            f.close()
            self.assertAlmostEqual(precision, 0.9988181487368183, delta=1e-5)

            os.system('rm '+output+'labeled_pdf.txt')

            # Check temporal precision
            # make model 1
            nodes, graph, graph_prob, graph_scores = spatiotemporal.create_DAG(
                state_dict,
                input_dir=data_path,
                scorestr='_scores1.log',
                output_dir=output,
                spatio_temporal_rule=True,
                expected_subcomplexes=subcomplexes,
                score_comp=True,
                exp_comp_map=exp_comp,
                out_cdf=False,
                draw_dag=False  # there is no need for heatmap
            )
            os.system('mv '+output+'labeled_pdf.txt '+output+'model1_pdf.txt ')
            # make model 2
            nodes, graph, graph_prob, graph_scores = spatiotemporal.create_DAG(
                state_dict,
                input_dir=data_path,
                scorestr='_scores2.log',
                output_dir=output,
                spatio_temporal_rule=True,
                expected_subcomplexes=subcomplexes,
                score_comp=True,
                exp_comp_map=exp_comp,
                out_cdf=False,
                draw_dag=False  # there is no need for heatmap
            )
            os.system('mv '+output+'labeled_pdf.txt '+output+'model2_pdf.txt ')
            analysis.temporal_precision(output+'model1_pdf.txt',output+'model2_pdf.txt',output_fn=output+'temporal_prevision.txt')
            f = open(output + 'temporal_prevision.txt', 'r')
            # First line is description
            f.readline()
            # 2nd line is the temporal precision
            trj_comp = float(f.readline())
            f.close()
            self.assertAlmostEqual(trj_comp, 0.9995202597385684, delta=1e-5)

if __name__ == '__main__':
    unittest.main()








