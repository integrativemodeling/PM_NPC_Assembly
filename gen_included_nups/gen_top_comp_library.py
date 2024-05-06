""" Script that determines which nup compositions are good scoring,
and generates files to prepare these simulations for structural sampling.

This file computes a likelihood for each possible composition and writes a list
of the top scoring compositions to a directory which is in turn used as a
library to seed subsequent configurational sampling simulations.

The likelihood model is the following. The data we use is from FCS
measurements of the bulk copy number. We make two assumptions about the
distribution of this data. 1) that the single pore copy number is Gaussian
distributed around it's mean and 2) that the estimate of the mean is subject to
a Gaussian uncertainty represented by the form of the prior on that mean with a
spread estimated from the individual traces.

Using these model assumptions, we build a single-pore copy number likelihood by
first computing the variance in estimates of the mean copy number from the
individual Nup kinetics data and marginalize that prior against the Gaussian
likelihood of single pore copy numbers with that spread fixed to a given
arbitrary value.

Different localizations of the Y-complex are enumerated,
while the localizations of other subcomplexes are assumed to be insignificant, and are assigned randomly.
The connecting complexes due not have FCS data - therefore, we assume they are added to the assembling pore
at the same time as the first nuclear ring on the same side of the NPC.

Listing of the nup subcomplex groupings:

YC:
    - Nup133
    - Nup107
    - Nup96
    - Sec13
    - Seh1
    - Nup85
    - Nup43
    - Nup160
    - Nup37

core1,3:
    - Nup93
    - Nup205
    - Nup155

core2,4:
    - Nup93
    - Nup188
    - Nup155

chan:
    - p54
    - p58_p45
    - p62

conn:
    - Nup155

available tag data:
    - Nup107
    - Nup153
    - Nup93
    - Nup62
    - Nup214
    - Nup358
    - TPR
    - Seh1
    - Pom121
    - Nup205
"""

from itertools import combinations
import numpy as np
import itertools
import sys
import os
import math
import pandas as pd

# params
outdir='/.../included_nups/'
data_dir = "/.../data/qfluor_data/"
nmodels = 4

times = ["5min", "6min", "8min", "10min", "15min"]

nups = ["Nup107", "Nup62", "Nup93", "Nup205"]

# scaling factors for each nup
cp = {"Nup107" : 31.3454,
      "Nup62" : 46.4844,
      "Nup205" : 16.4733,
      "Nup93" : 47.2745}


# enumerate all possible compositions (independent of time)
yc = range(5) # nup107 tag
core13 = range(3) # nup205 tag
core24 = range(3) # nup93 tag
chan = range(5) # nup62 tag
conn = range(3) # no tag

# generate all possible combinations filtering states where core1,3 is greater or equalt to core2,4
all_library = [l for l in itertools.product(yc, core13, core24, chan) if l[1] >= l[2]]
# remove empty state
all_library.pop(all_library.index((0,0,0,0)))

# Function to convert a series of time based intensities to mean and standard deviation as a function of time
def convert_qflour(fluor_data_dir,nups,final_cn,times,output_dir):
    # read in available csv files
    nd = {nup: pd.read_csv(fluor_data_dir + "total_" + nup + "_homoZ.csv") for nup in nups}
    sigmas = {key: [] for key in nd.keys()}
    means = {key: [] for key in nd.keys()}
    # Read in data. Use 30 min for mature
    for t in [5, 6, 8, 10, 15, 30]:
        # for each protein
        for key in nd.keys():
            m5_data = nd[key].loc[nd[key]["Time (min)"] == t]
            for m5_key in m5_data.keys():
                for value in m5_data.loc[:, m5_key]:
                    # replace negative intensity with 0
                    if value < 0:
                        m5_data = m5_data.replace(value, 0.0)
            # Read in mean and std_dev from the data
            traces = nd[key].columns[2:]  # from the CSV, the row index is stored as well
            cn = m5_data[traces] * final_cn[key]
            means[key].append(float(cn.mean(1)))
            sigmas[key].append(float(cn.std(1)))
    # convert means and std dictionaries into pd dataframes and save to csv
    for key in nd.keys():
        nup_cn=pd.DataFrame(index=times,columns=['Time','mean','std'])
        for i in range(len(times)):
            nup_cn['Time'][times[i]]=times[i]
            nup_cn['mean'][times[i]]=means[key][i]
            nup_cn['std'][times[i]] = sigmas[key][i]
        # save to csv
        output_fn=output_dir+'exp_comp_'+key+'.csv'
        nup_cn.to_csv(output_fn)
    return

def composition_likelihood_function(mean, std, t, prots, model_cn):
    """Function that calculates the likelihood of an individual node, used by
    calc_likelihood().

    @param mean: dictionary of dictionaries where the first key is the protein,
           the second key is the time, and the expected mean copy number
           from experiment is returned.
    @param std: dictionary of dictionaries where the first key is the protein,
           the second key is the time, and the expected standard deviation
           of protein copy number from experiment is returned.
    @param t: time to find the mean / std
    @param prots: list of proteins or subcomplexes which will be scored
           according to this likelihood function
    @param model_cn: dictionary with the copy numbers of each protein in the forward model
    @return w: float, the weight of the graphNode according to the composition
            likelihood function.
    """
    w = 0
    for prot in prots:
        # x counts the number of proteins of a given type in the node
        x = model_cn[prot]
        # check std is greater than 0
        if std[prot][t] > 0:
            pass
        else:
            warnings.warn(
                'WARNING!!! Standard deviation of protein ' + prot
                + ' 0 or less at time ' + t
                + '. May lead to illogical results.')
        w += (0.5 * ((x - mean[prot][t]) / std[prot][t])**2 + np.log(std[prot][t] * np.sqrt(2 * np.pi)))
    return w

def calc_likelihood(exp_comp_map, t, state):
    """
    Function that adds a score for the compositional likelihood for all
    states represented as nodes in the graph. The composition likelihood
    assumes a Gaussian distribution for copy number of each protein or
    subcomplex with means and standard deviatiations derived from experiment.
    Returns the nodes, with the new weights added.

    @param exp_comp_map: dictionary, which describes protein stoicheometery.
           The key describes the protein, which should correspond to names
           within the expected_subcomplexes. Only copy numbers for proteins
           or subcomplexes included in this dictionary will be scored. For
           each of these proteins, a csv file should be provided with protein
           copy number data. The csv file should have 3 columns,
           1) "Time", which matches up to the possible times in the graph,
           2) "mean", the average protein copy number at that time point
           from experiment, and 3) "std", the standard deviation of that
           protein copy number from experiment.
    @param nodes: list of graphNode objects, which have been already been
           initiated with static scores
    @return nodes: editted list of graphNode objects, which now have static
            and composition scores
    """
    import pandas as pd
    prots=exp_comp_map.keys()
    # Data is stored as a dictionary of dictionaries. The first dictionary
    # references which protein you are refering to.
    # the 2nd dictionary references which time you are refering to. The return
    # is the mean or standard deviation of the protein copy number
    mean = {}
    std = {}
    # import csv file as pandas data frame
    for prot in prots:
        prot_dict_mean = {}
        prot_dict_std = {}
        if os.path.exists(exp_comp_map[prot]):
            exp = pd.read_csv(exp_comp_map[prot])
        else:
            raise Exception(
                "Error!!! Check exp_comp_map. Unable to find composition "
                "file: " + exp_comp_map[prot] + '\nClosing...')
        for i in range(len(exp)):
            prot_dict_mean[exp['Time'][i]] = exp['mean'][i]
            prot_dict_std[exp['Time'][i]] = exp['std'][i]
        mean[prot] = prot_dict_mean
        std[prot] = prot_dict_std
    # compute the compositional likelihood of the nodes. Convert state into dictionary of copy numbers
    state_cn={'Nup107':8.0*state[0],'Nup93':8.0*(state[1]+state[2]),'Nup205':8.0*(state[1]),'Nup62':8.0*state[3]}
    weight = composition_likelihood_function(mean, std, t, prots, state_cn)
    return weight


convert_qflour(data_dir,nups,cp,times,outdir)

os.chdir(outdir)
# compare the complete composition library against the composition model 
# at each time and extract top scoring compositions
for time in times:
    exp_comp = {'Nup107': 'exp_comp_Nup107.csv', 'Nup205': 'exp_comp_Nup205.csv', 'Nup62': 'exp_comp_Nup62.csv',
                'Nup93': 'exp_comp_Nup93.csv'}
    unnormalized_weights =[]
    for state in all_library:
        unnormalized_weights.append(calc_likelihood(exp_comp,time,state))
        #print(state)
        #print(calc_likelihood(exp_comp,time,state))

    unw = np.array(unnormalized_weights)
    # get top scoring nmodels
    mindx = np.argsort(unw)[0:nmodels]
    print(time)

    # write out library
    olist = []
    nuplist = []
    for m in mindx:
        state = all_library[m]
        print('State:')
        print(state)
        print('Weight:')
        print(unw[m])

        # convert state counts to nup list
        # enumerate of yc locations
        for j in combinations(["yc_inner_cr", "yc_inner_nr", "yc_outer_cr", "yc_outer_nr"], state[0]):
            # append to state log file for each combination
            olist.append([state[0], state[1], state[2], state[3]])

            npl = []
            npl += list(j)
            npl += ["ir_core_1", "ir_core_3"][0:state[1]]
            npl += ["ir_core_2", "ir_core_4"][0:state[2]]
            npl += ["ir_chan_1", "ir_chan_2", "ir_chan_3", "ir_chan_4"][0:state[3]]
            # Add conn complex with the yc of the same side. conn_1 - cr, conn_2 - nr
            if "yc_inner_cr" in npl or "yc_outer_cr" in npl:
                npl += ["conn_1"]
            if "yc_inner_nr" in npl or "yc_outer_nr" in npl:
                npl += ["conn_2"]

            nuplist.append(npl)


    # write top "scoring" compositions to file
    oary = np.array(olist, dtype=int)
    np.savetxt(outdir + "/" + time + ".txt", oary)

    # write nup config library
    for indx,npl in enumerate(nuplist):
        with open(outdir + "/" + str(indx+1) + "_" + time + ".config", "w") as fh:
            for sc in npl:
                fh.write(sc +"\n")

# finally, write out the fully assembled state for the mature simulation
np.savetxt(outdir + "/mature.txt", np.array([[4, 2, 2, 4,]]))

# write nup config library
with open(outdir + "/1_mature.config", "w") as fh:
    npl = ["yc_inner_cr", "yc_inner_nr", "yc_outer_cr", "yc_outer_nr",
           "ir_core_1", "ir_core_3",
           "ir_core_2", "ir_core_4",
           "ir_chan_1", "ir_chan_2", "ir_chan_3", "ir_chan_4",
           "conn_1", "conn_2"]

    for sc in npl:
        fh.write(sc +"\n")
