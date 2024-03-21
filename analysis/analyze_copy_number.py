"""Script to analyze nup copy number of a model.
Copy numbers are computed as the weighted sum of Nup copy numbers from each compositional state at each time point.
Requires a labeled_pdf (sys.argv[1]). The Nup_folder variable should direct to the folder where the included_nups are located"""
import sys
import os
import math
import numpy

labeled_pdf=sys.argv[1]

# read in labeled_pdf as a dictionary, with keys being each trajectory and values being the probability of that trajecotry
def read_labeled_pdf(pdf_file):
    # create blank dictonary to store the results
    prob_dict = {}
    # read in labeled pdf file
    old = open(pdf_file, 'r')
    line = old.readline()
    # store the path through various nodes, as well as the probability of that path
    while line:
        line_split = line.split()
        # assumes the first string is the trajectory string, the second string is the probability
        if len(line_split) > 1:
            # use # for comments
            if line_split[0]=='#':
                pass
            else:
                trj = line_split[0]
                prob = float(line_split[1])
                # store in dictionary
                prob_dict[trj] = prob
        line = old.readline()
    old.close()
    return prob_dict

# For a trajectory, returns an array of Nup copy numbers as a function of time. The first value loops over the time, the second value loops over the Nup
def copy_number_from_state(state):
    # Directory for included_nups
    Nup_folder = '/.../included_nups/'

    state_list=state.split('|')
    state_list=state_list[:-1]

    N = len(state_list)
    # Map from index to nup: 0 - n107, 1- n93, 2- n205, 3- n62, 4-Seh1, 5- n188
    Nups = numpy.zeros((N, 6))

    # Grab nups from included_nups file
    for i in range(0, N):
        nup_file = Nup_folder + state_list[i] + '.config'
        to_read = open(nup_file, 'r')
        line = to_read.readline()
        while line:
            if "yc" in line:
                Nups[i, 0] = Nups[i, 0] + 8
            if "ir_core" in line:
                Nups[i, 1] = Nups[i, 1] + 8
            if ("ir_core_1" in line) or ("ir_core_3" in line):
                Nups[i, 2] = Nups[i, 2] + 8
            if "ir_chan" in line:
                Nups[i, 3] = Nups[i, 3] + 8
            if "yc" in line:
                Nups[i, 4] = Nups[i, 4] + 8
            if ("ir_core_2" in line) or ("ir_core_4" in line):
                Nups[i, 5] = Nups[i, 5] + 8
            line = to_read.readline()

    return Nups,N

# Read in labeled_pdf file into a dictionary. Each trajectory is listed as a dictionary, with keys as the trajectory and the values as the probability of that trajectory
prob_dict=read_labeled_pdf(labeled_pdf)


# Loop over the full dictionary. Create a list with 2 values: 1) the probability of the state, 2) the Nup copy number of that state.
key_list=prob_dict.keys()
Nup_prob=[]
for key in key_list:
    CN,N_times=copy_number_from_state(key)
    Nup_prob.append([prob_dict[key],CN])
# Calculate the weight normalization terms
V1=0
V2=0
for state in Nup_prob:
    V1+=state[0]
    V2+=state[0]*state[0]
# Caclulate the copy number of each Nup as a function of time
# Nup107  - index 0 ----------------------------------------------------------------------------------------------------
copy_number=numpy.zeros((N_times,2))
index=0
# calculate mean
for state in Nup_prob:
    for i in range(N_times):
        copy_number[i,0]+=state[0]*state[1][i][index]
# calculate std
for state in Nup_prob:
    for i in range(N_times):
        copy_number[i,1]+=state[0]*((state[1][i][index]-copy_number[i,0])**2)
    copy_number[i,1]=copy_number[i,1]/(V1-(V2/V1))
numpy.savetxt('Nup107_copy_number.txt',copy_number)

# Nup107  - index 1 ----------------------------------------------------------------------------------------------------
copy_number=numpy.zeros((N_times,2))
index=1
# calculate mean
for state in Nup_prob:
    for i in range(N_times):
        copy_number[i,0]+=state[0]*state[1][i][index]
# calculate std
for state in Nup_prob:
    for i in range(N_times):
        copy_number[i,1]+=state[0]*((state[1][i][index]-copy_number[i,0])**2)
    copy_number[i,1]=copy_number[i,1]/(V1-(V2/V1))
numpy.savetxt('Nup93_copy_number.txt',copy_number)

# Nup205  - index 2 ----------------------------------------------------------------------------------------------------
copy_number=numpy.zeros((N_times,2))
index=2
# calculate mean
for state in Nup_prob:
    for i in range(N_times):
        copy_number[i,0]+=state[0]*state[1][i][index]
# calculate std
for state in Nup_prob:
    for i in range(N_times):
        copy_number[i,1]+=state[0]*((state[1][i][index]-copy_number[i,0])**2)
    copy_number[i,1]=copy_number[i,1]/(V1-(V2/V1))
numpy.savetxt('Nup205_copy_number.txt',copy_number)

# Nup62  - index 3 ----------------------------------------------------------------------------------------------------
copy_number=numpy.zeros((N_times,2))
index=3
# calculate mean
for state in Nup_prob:
    for i in range(N_times):
        copy_number[i,0]+=state[0]*state[1][i][index]
# calculate std
for state in Nup_prob:
    for i in range(N_times):
        copy_number[i,1]+=state[0]*((state[1][i][index]-copy_number[i,0])**2)
    copy_number[i,1]=copy_number[i,1]/(V1-(V2/V1))
numpy.savetxt('Nup62_copy_number.txt',copy_number)

# Seh1  - index 4 ----------------------------------------------------------------------------------------------------
copy_number=numpy.zeros((N_times,2))
index=4
# calculate mean
for state in Nup_prob:
    for i in range(N_times):
        copy_number[i,0]+=state[0]*state[1][i][index]
# calculate std
for state in Nup_prob:
    for i in range(N_times):
        copy_number[i,1]+=state[0]*((state[1][i][index]-copy_number[i,0])**2)
    copy_number[i,1]=copy_number[i,1]/(V1-(V2/V1))
numpy.savetxt('Seh1_copy_number.txt',copy_number)

# Nup188  - index 5 ----------------------------------------------------------------------------------------------------
copy_number=numpy.zeros((N_times,2))
index=5
# calculate mean
for state in Nup_prob:
    for i in range(N_times):
        copy_number[i,0]+=state[0]*state[1][i][index]
# calculate std
for state in Nup_prob:
    for i in range(N_times):
        copy_number[i,1]+=state[0]*((state[1][i][index]-copy_number[i,0])**2)
    copy_number[i,1]=copy_number[i,1]/(V1-(V2/V1))
numpy.savetxt('Nup188_copy_number.txt',copy_number)


