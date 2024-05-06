"""
Calculates the model precision (purity) of the full assembly model.
Make sure the model directory is correctly specified.
"""
# import relevant modules
import sys
import IMP.spatiotemporal as spatiotemporal
from IMP.spatiotemporal import analysis
import os

input_dir='/.../simulations_round2/Refined_energies_1model_460/filtered_noNup188/'
# go to input_dir
os.chdir(input_dir)

output_dir3=input_dir+'total/'

analysis.purity(output_dir3+'labeled_pdf.txt',output_fn = output_dir3+'purity.txt')
