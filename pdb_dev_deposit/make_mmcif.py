"""
Script to generate an mmCIF file suitable for deposition in PDB-dev.

Entities correspond to unique protein sequences.

Asymmetric units correspond to nucleoporins.

multiple AS can point to the same entity.

Each snapshot is represented as a single assembly object, from which a model is generated.
1 model (the centroid structure from the most representative cluster) is chosen per time point.
Therefore, each state and each model group consist of a single model
All models are then included in the state_group,
and the ordered process of NPC assembly is represented as a series of model groups
"""

import RMF
import IMP
import IMP.rmf
import ihm
import ihm.representation
import ihm.protocol
import ihm.dumper
import ihm.citations
import ihm.reference
import ihm.analysis
import ihm.startmodel
from Bio import SeqIO


main_dir='../'

times=['5min','6min','8min','10min','15min','mature']
best_states={'5min':'2_5min','6min':'10_6min','8min':'14_8min',
             '10min':'4_10min','15min':'8_15min','mature':'1_mature'}
# Time dependent parameters for NPC assembly. Includes results and input information
CCC={'5min':0.6661717242187478,'6min':0.6034028217725574,'8min':0.7963839913096493,
     '10min':0.7757774850030276,'15min':0.7618849068526939,'mature':0.7542426645010205}
nstructures={'5min':1121,'6min':881,'8min':801,'10min':801,'15min':1121,'mature':801}
structural_precision={'5min':232.552,'6min':144.083,'8min':104.579,'10min':99.893,'15min':100.826,'mature':104.258}
membrane_D={'5min':515.3,'6min':583.9,'8min':727.4,'10min':845.9,'15min':798.3,'mature':870.0}
membrane_H={'5min':427.2,'6min':424.9,'8min':429.9,'10min':405.4,'15min':342.5,'mature':300.0}
# Dictionary with the starting residue number for each protein modeled
start_res = {'Nup133': 518, 'Nup107': 150, 'Nup96': 277, 'SEC13': 14,
             'SEH1': 1, 'Nup85': 9, 'Nup43': 4, 'Nup160': 41,
             'Nup37': 9, 'Nup93': 1, 'Nup205': 9, 'Nup155': 20,
             'Nup188': 1, 'p54': 128, 'p58': 248, 'p62': 334}
# Dictionary with the ending residue number for each protein modeled
end_res = {'Nup133': 1156, 'Nup107': 924, 'Nup96': 751, 'SEC13': 304,
           'SEH1': 322, 'Nup85': 475, 'Nup43': 380, 'Nup160': 1201,
           'Nup37': 326, 'Nup93': 815, 'Nup205': 1692, 'Nup155': 1375,
           'Nup188': 1564, 'p54': 493, 'p58': 418, 'p62': 502}

nstates=len(times)

title = ("Integrative spatiotemporal modeling of biomolecular processes: "
         "application to the assembly of the Nuclear Pore Complex")
system = ihm.System(title=title)
# start IMP model
IMP_m = IMP.Model()

# Function to load hierarchy from RMF. By default, loads the last frame
def load_hc(best_state):
    """Load the first frame of the centroid model for a given state (best state)
    """
    # path to cluster centroid
    rmf_filename=(main_dir+'simulations_round2/Refined_energies_1model_460/filtered_noNup188/total/'
                  +best_state+'/cluster.0/cluster_center_model.rmf3')
    rmf_fh = RMF.open_rmf_file_read_only(rmf_filename)
    h = IMP.rmf.create_hierarchies(rmf_fh, IMP_m)[0]
    lowest_frame_id = rmf_fh.get_number_of_frames() - 1
    IMP.rmf.load_frame(rmf_fh, RMF.FrameID(lowest_frame_id))
    return h

print('Additing software and citations...')
# Include citation for NPC assembly paper
system.citations.append(ihm.Citation(
          pmid='XXX', title=title,
          journal="XXX", volume=1, page_range=(1,1),
          year=2024,
          authors=['Latham, A.P.', 'Tempkin, J.O.B.', 'Otsuka, S.',
                   'Zhang, W.', 'Ellenberg, J.', 'Sali, A.'],
          doi='XXX'))

# Modeling was performed in IMP
imp = ihm.Software(
          name="Integrative Modeling Platform (IMP)",
          version="2.20",
          classification="Integrative model building.",
          description="Used for representation, scoring, search process, and assessment of the model.",
          location='https://integrativemodeling.org')
if hasattr(ihm, 'citations'):
    imp.citation = ihm.citations.imp
system.software.append(imp)

# Fitting to GMM was performed with
imp = ihm.Software(
          name="gmconvert",
          version="2022/05/09",
          classification="Processing of input data to modeling.",
          description="Program used to fit GMMs to 3DEM maps at each snapshot as "
                      "well as fitting of forward model GMMs used in the restraints.",
          location='https://pdbj.org/gmfit/doc_gmconvert/README_gmconvert.html')
system.software.append(imp)
print('Done.')

# Read in experimental datasets. Make list with different data at each time point,
print('Adding experimental data...')
# as ET data changes as a function of time
exp_data=[]
processed_EM_list=[]
# Original EM
em_database_loc = ihm.location.EMPIARLocation("EMD-3820")
raw_em = ihm.dataset.EMDensityDataset(em_database_loc,details='Original electron tomography dataset.')
# Original FCS
FCS_database_loc=ihm.location.DatabaseLocation("idr0115",details=
    'Raw FCS data available on the Image Data Resource: https://idr.openmicroscopy.org/')
raw_FCS = ihm.dataset.Dataset(FCS_database_loc,details='Original fluorescence correlation spectroscopy dataset.')
# Processed FCS
FCS_excel_loc=ihm.location.InputFileLocation(main_dir+"data/qfluor_data/Exp_data.txt",repo=None)
processed_FCS = ihm.dataset.Dataset(FCS_excel_loc,details='Processed fluorescence correlation spectroscopy dataset.')
processed_FCS.parents.append(raw_FCS)
# PDB
pdb1_loc=ihm.location.PDBLocation("5a9q")
pdb_structure1=ihm.dataset.PDBDataset(pdb1_loc,details='PDB for the Y-complex and connecting complex')
pdb2_loc=ihm.location.PDBLocation("5ijo")
pdb_structure2=ihm.dataset.PDBDataset(pdb2_loc,details='PDB for the inner ring')
# Processed EM
for i in range(0,len(times)):
    em_processed_loc=ihm.location.InputFileLocation(main_dir+
        'data/fit_etdata_with_gmm/Andrew_run/data/'+times[i]+'_150.mrc',repo=None)
    processed_EM=ihm.dataset.EMDensityDataset(em_processed_loc,details='Processed EM data at '+times[i])
    processed_EM.parents.append(raw_em)
    processed_EM_list.append(processed_EM)
    # append all experimental data
    exp_data.append(ihm.dataset.DatasetGroup((raw_em,raw_FCS,processed_FCS,pdb_structure1,pdb_structure2,processed_EM),
                                              name='Snapshot model data for '+times[i]))
print('Done.')

# Include entities for each protein sequence
print('Adding entities...')
def build_entity_template(hc_tmpl,start_res,end_res):
    """Return an entity object for a Nup domain.
    """
    pdb_path='../data/cg_models/'
    # Dictionary converting protein names to Uniprot entries
    Uniprot_dict={'Nup133':'Q8WUM0', 'Nup107':'P57740', 'Nup96':'P52948', 'SEC13':'P55735',
                  'SEH1':'Q96EE3', 'Nup85':'Q9BW27', 'Nup43':'Q8NFH3', 'Nup160':'Q12769',
                  'Nup37':'Q8NFH4', 'Nup93':'Q8N1F7', 'Nup205':'Q92621', 'Nup155':'O75694',
                  'Nup188':'Q5SRE5', 'p54':'Q7Z3B4', 'p58':'Q9BVL2', 'p62':'P37198'}
    # Dictionary with pdbID and chainID for each protein
    pdb_label = {'Nup133': ['5a9q','L'], 'Nup107': ['5a9q','M'], 'Nup96': ['5a9q','N'], 'SEC13': ['5a9q','O'],
                    'SEH1': ['5a9q','P'], 'Nup85': ['5a9q','Q'], 'Nup43': ['5a9q','R'], 'Nup160': ['5a9q','S'],
                    'Nup37': ['5a9q','T'], 'Nup93': ['5ijo','C'], 'Nup205': ['5ijo','D'], 'Nup155': ['5ijo','E'],
                    'Nup188': ['5ijo','J'], 'p54': ['5ijo','R'], 'p58': ['5ijo','S'], 'p62': ['5ijo','T']}
    # Loop over all proteins
    entities_dict = {}
    for subcomplex in hc_tmpl.get_children():
        for template in subcomplex.get_children():
            name = template.get_name().split("_")[0]
            if name not in entities_dict.keys():

                # Get sequence for correct uniprot entity depending on the name
                ref = ihm.reference.UniProtSequence.from_accession(Uniprot_dict[name])
                #ref.alignments.append(ihm.reference.Alignment(db_begin=start_res[name],entity_begin=start_res[name],db_end=end_res[name],entity_end=end_res[name]))
                # UniprotID is for the Nup98-Nup96 fusion protein. Offset the alignment by the length of Nup98
                if name=='Nup96':
                    ref.alignments.append(ihm.reference.Alignment(db_begin=start_res[name]+880,db_end=end_res[name]+880))
                else:
                    ref.alignments.append(ihm.reference.Alignment(db_begin=start_res[name], db_end=end_res[name]))
                query_seqres = SeqIO.PdbIO.PdbSeqresIterator(pdb_path + pdb_label[name][0] + '.pdb')
                count = 0
                for chain in query_seqres:
                    if chain.id[-1] == pdb_label[name][1]:
                        sequence = chain.seq[start_res[name]-1:end_res[name]]
                        count += 1
                if count == 1:
                    entity = ihm.Entity(sequence, description="_".join([name, subcomplex.get_name()]), references=[ref])
                else:
                    print('Error! Check pdb file, ' + str(count) + ' sequences were found for ' + name)

                entities_dict[name] = entity
                system.entities.append(entity)
    return entities_dict
# Load mature hierarchy
hc_mature=load_hc(best_states[times[-1]])
# Set entities from mature hierarchy
possible_entities=build_entity_template(hc_mature,start_res,end_res)
print('Done.')

# Define asymeteric units for each state
print('Building asymmetric units...')
def build_new_assembly_from_entities(hc_sc, edict, syst, chain_offset, asym_unit_map=None):
    """Return an asymmetric units and representation for each subcomplex
    """
    # Dictionary that connects each protein to its location in the original PDB
    # (structure1 - 5a9q, structure2 - 5ijo)
    start_model_label = {'Nup133_yc_inner_cr': [pdb_structure1, '3'], 'Nup107_yc_inner_cr': [pdb_structure1, '4'],
                         'Nup96_yc_inner_cr': [pdb_structure1, '5'], 'SEC13_yc_inner_cr': [pdb_structure1, '6'],
                         'SEH1_yc_inner_cr': [pdb_structure1, '7'], 'Nup85_yc_inner_cr': [pdb_structure1, '8'],
                         'Nup43_yc_inner_cr': [pdb_structure1, '9'], 'Nup160_yc_inner_cr': [pdb_structure1, 'a'],
                         'Nup37_yc_inner_cr': [pdb_structure1, 'b'], 'Nup133_yc_outer_cr': [pdb_structure1, 'U'],
                         'Nup107_yc_outer_cr': [pdb_structure1, 'V'], 'Nup96_yc_outer_cr': [pdb_structure1, 'W'],
                         'SEC13_yc_outer_cr': [pdb_structure1, 'X'], 'SEH1_yc_outer_cr': [pdb_structure1, 'Y'],
                         'Nup85_yc_outer_cr': [pdb_structure1, 'Z'], 'Nup43_yc_outer_cr': [pdb_structure1, '0'],
                         'Nup160_yc_outer_cr': [pdb_structure1, '1'], 'Nup37_yc_outer_cr': [pdb_structure1, '2'],
                         'Nup133_yc_inner_nr': [pdb_structure1, 'L'], 'Nup107_yc_inner_nr': [pdb_structure1, 'M'],
                         'Nup96_yc_inner_nr': [pdb_structure1, 'N'], 'SEC13_yc_inner_nr': [pdb_structure1, 'O'],
                         'SEH1_yc_inner_nr': [pdb_structure1, 'P'], 'Nup85_yc_inner_nr': [pdb_structure1, 'Q'],
                         'Nup43_yc_inner_nr': [pdb_structure1, 'R'], 'Nup160_yc_inner_nr': [pdb_structure1, 'S'],
                         'Nup37_yc_inner_nr': [pdb_structure1, 'T'], 'Nup133_yc_outer_nr': [pdb_structure1, 'C'],
                         'Nup107_yc_outer_nr': [pdb_structure1, 'D'], 'Nup96_yc_outer_nr': [pdb_structure1, 'E'],
                         'SEC13_yc_outer_nr': [pdb_structure1, 'F'], 'SEH1_yc_outer_nr': [pdb_structure1, 'G'],
                         'Nup85_yc_outer_nr': [pdb_structure1, 'H'], 'Nup43_yc_outer_nr': [pdb_structure1, 'I'],
                         'Nup160_yc_outer_nr': [pdb_structure1, 'J'], 'Nup37_yc_outer_nr': [pdb_structure1, 'K'],
                         'Nup93_ir_core_1': [pdb_structure2, 'C'], 'Nup205_ir_core_1': [pdb_structure2, 'D'],
                         'Nup155_ir_core_1': [pdb_structure2, 'E'],
                         'Nup93_ir_core_2': [pdb_structure2, 'I'], 'Nup188_ir_core_2': [pdb_structure2, 'J'],
                         'Nup155_ir_core_2': [pdb_structure2, 'K'],
                         'Nup93_ir_core_3': [pdb_structure2, 'O'], 'Nup205_ir_core_3': [pdb_structure2, 'P'],
                         'Nup155_ir_core_3': [pdb_structure2, 'Q'],
                         'Nup93_ir_core_4': [pdb_structure2, 'U'], 'Nup188_ir_core_4': [pdb_structure2, 'V'],
                         'Nup155_ir_core_4': [pdb_structure2, 'W'],
                         'p54_ir_chan_1': [pdb_structure2, 'F'], 'p58_ir_chan_1': [pdb_structure2, 'G'],
                         'p62_ir_chan_1': [pdb_structure2, 'H'],
                         'p54_ir_chan_2': [pdb_structure2, 'L'], 'p58_ir_chan_2': [pdb_structure2, 'M'],
                         'p62_ir_chan_2': [pdb_structure2, 'N'],
                         'p54_ir_chan_3': [pdb_structure2, 'R'], 'p58_ir_chan_3': [pdb_structure2, 'S'],
                         'p62_ir_chan_3': [pdb_structure2, 'T'],
                         'p54_ir_chan_4': [pdb_structure2, 'X'], 'p58_ir_chan_4': [pdb_structure2, 'Y'],
                         'p62_ir_chan_4': [pdb_structure2, 'Z'],
                         'Nup155_conn_1': [pdb_structure1, 'B'], 'Nup155_conn_2': [pdb_structure1, 'A']}

    nups = [child for child in hc_sc.get_children() if "Density" not in child.get_name()]

    asym_units = []
    sc_representation = []
    for nup in nups:
        name = nup.get_name().split("_")[0]
        unique_name = "_".join([name, hc_sc.get_name()])
        asu = ihm.AsymUnit(edict[name], auth_seq_id_map=chain_offset[name]-1, details=unique_name)
        nres = len(nup.get_children())
        spoke_name = unique_name.replace('spoke_1_', '').replace('spoke_2_', '').replace('spoke_3_', '').replace(
            'spoke_4_', '').replace('spoke_5_', '').replace('spoke_6_', '').replace('spoke_7_', '').replace('spoke_8_','')
        sm = ihm.startmodel.StartingModel(asu, start_model_label[spoke_name][0], start_model_label[spoke_name][1],
                                          offset=chain_offset[name]-1, description=spoke_name)
        rep = ihm.representation.FeatureSegment(asu,
                                                rigid=True,
                                                primitive="sphere",
                                                count=nres,
                                                starting_model = sm)
        sc_representation.append(rep)
        asym_units.append(asu)

        if asym_unit_map is not None:
            asym_unit_map[unique_name] = asu

    syst.asym_units.extend(asym_units)

    return sc_representation, asym_units

# Determine the overall representation and asymmetric units from the mature state
# Variables to save
assemblies_list=[]
asym_unit_map={}
# variables for output
representations=[]
asym=[]
# For each subcomplex
for sc in hc_mature.get_children():
    sc_rep, sc_asym = build_new_assembly_from_entities(sc, possible_entities, system, start_res, asym_unit_map=asym_unit_map)
    # Add each subcomplex to the appropriate variables
    representations.extend(sc_rep)
    asym.extend(sc_asym)
representation=ihm.representation.Representation(representations,name='Representation of '
                                                                      'Nups in NPC assembly snapshot models')

# Determine the list of assemblies at each timepoint
# Loop over all times
for i in range(nstates):
    asym_list=[]
    # Loop over all subcomplexes at each timepoint
    hc_temp=load_hc(best_states[times[i]])
    for sc in hc_temp.get_children():
        # For each protein in the subcomplex
        for nup in sc.get_children():
            # Use unique name to find the corresponding asymmetric unit
            name = nup.get_name().split("_")[0]
            unique_name = "_".join([name, sc.get_name()])
            asym_list.append(asym_unit_map[unique_name])
    # Create assembly for that timepoint
    assemblies_list.append(ihm.Assembly(asym_list,name='Assembly at '+times[i]))
print('Done')


# Define sampling of single snapshot model
print('Adding protocol...')
def snapshot_model_protocol(t,exp,assembly):
    """Return protocol object used to model each snapshot
    """
    # Define main_dir
    main_dir='../'
    # Describe the modeling protocol
    prot = ihm.protocol.Protocol(name='Snapshot modeling at '+t)
    # Step 1 of Monte Carlo, with filtering after:
    prot.steps.append(ihm.protocol.Step(
        assembly=assembly, dataset_group=exp,
        method='MC sampling',
        name='Monte Carlo sampling step 1. Repeated for each snapshot model.',
        num_models_begin=200, num_models_end=1600000,
        script_file=ihm.location.WorkflowFileLocation(main_dir + 'start_sim/main.py'),
        ordered=True, software=imp))
    prot.steps.append(ihm.protocol.Step(
        assembly=assembly, dataset_group=exp,
        method='Filtering',
        name='Choose the lowest scoring structure from each simulation.',
        num_models_begin=1600000, num_models_end=200,
        script_file=ihm.location.WorkflowFileLocation(main_dir + 'start_sim/extract_lowestE.py'),
        software=imp))
    # Step 2 of Monte Carlo, with filtering after:
    prot.steps.append(ihm.protocol.Step(
        assembly=assembly, dataset_group=exp,
        method='MC sampling',
        name='Monte Carlo sampling step 2. Repeated for each snapshot model.',
        num_models_begin=200, num_models_end=320000,
        script_file=ihm.location.WorkflowFileLocation(main_dir + 'refine_sim/main_restart1.py'),
        ordered=True, software=imp))
    prot.steps.append(ihm.protocol.Step(
        assembly=assembly, dataset_group=exp,
        method='Filtering',
        name='Choose the lowest scoring structure from each simulation.',
        num_models_begin=320000, num_models_end=200,
        script_file=ihm.location.WorkflowFileLocation(main_dir + 'refine_sim/extract_lowestE_step1.py'),
        software=imp))
    # Step 3 of Monte Carlo, with filtering after:
    prot.steps.append(ihm.protocol.Step(
        assembly=assembly, dataset_group=exp,
        method='MC sampling',
        name='Monte Carlo sampling step 3. Repeated for each snapshot model.',
        num_models_begin=200, num_models_end=320000,
        script_file=ihm.location.WorkflowFileLocation(main_dir + 'refine_sim/main_restart2.py'),
        ordered=True, software=imp))
    prot.steps.append(ihm.protocol.Step(
        assembly=assembly, dataset_group=exp,
        method='Filtering',
        name='Choose the lowest scoring structure from each '
                'replica, and further filter by the median score.',
        num_models_begin=320000, num_models_end=801,
        script_file=ihm.location.WorkflowFileLocation(main_dir + 'score_graph/prepare_filtered_noNup188.py'),
        ordered=True, software=imp))
    prot.steps.append(ihm.protocol.Step(
        assembly=assembly, dataset_group=exp,
        method='Trajectory construction',
        name='Score trajectories based on the scores of snapshot models and the transitions between them.',
        script_file=ihm.location.WorkflowFileLocation(main_dir + 'score_graph/noNup188/score_trj_noNup188.py'),
        ordered=True, software=imp, multi_state=True))
    return prot

# Create list for protocol of each snapshot
protocol_list=[]
# Loop over each state and determine protocol for producing that state
for i in range(nstates):
    protocol_list.append(snapshot_model_protocol(times[i],exp_data[i],assemblies_list[i]))
print('Done.')

# Define class to extract spheres
class MyModel(ihm.model.Model):
    """Create a unique class for model objects to read in coordinates from hierarchy
    """
    # override the get_spheres subclass to pull coordinates from rmf file on disk
    def get_spheres(self):
        # Load hierarchy
        hc = load_hc(self.name)
        # for each subcomplex
        for nups in hc.get_children():
            # for each nup
            for nup in nups.get_children():
                name = nup.get_name().split("_")[0]
                unique_name = "_".join([name, nups.get_name()])
                # Get the asymmetric unit for that Nup
                _asym = asym_unit_map[unique_name]
                # Read in leaves
                for leaf in IMP.atom.get_leaves(nup):
                    # parse info for sphere
                    x, y, z = IMP.core.XYZR(leaf).get_coordinates()
                    resids = IMP.atom.Fragment(leaf).get_residue_indexes()-(start_res[name]-1)
                    R = IMP.core.XYZR(leaf).get_radius()
                    yield ihm.model.Sphere(asym_unit=_asym,
                                           seq_id_range = (resids[0],resids[-1]),
                                           x=x,
                                           y=y,
                                           z=z,
                                           radius=R)

# function to find position restraint at a give time
def find_pos_restraint(_assemblies, chain_offset):
    """Function to calculate position restraints
    """
    # List of all assemblies in the current representation
    _asym_units_list=[]
    # Convert list of asymmetric unit objects to list of names
    for _assembly in _assemblies:
        _asym_units_list.append(_assembly.details)
    _pos_restraint_list=[]
    # Load the template rmf of the mature NPC
    template_rmf=main_dir+'data/cg_models/10/npc_cg.rmf'
    template_rmf_fh = RMF.open_rmf_file_read_only(template_rmf)
    template_h = IMP.rmf.create_hierarchies(template_rmf_fh, IMP_m)[0]
    lowest_frame_id = template_rmf_fh.get_number_of_frames() - 1
    IMP.rmf.load_frame(template_rmf_fh, RMF.FrameID(lowest_frame_id))
    # for each subcomplex
    for subcomplex in template_h.get_children():
        # Exclude EM particles
        if subcomplex.get_name() not in ["Data_density", "Sigma"]:
            # for each nup
            for nup in subcomplex.get_children():
                # Exclude EM particles
                if nup.get_name() =='Density':
                    continue
                name = nup.get_name().split("_")[0]
                unique_name = "_".join([name, subcomplex.get_name()])
                # Only include assemblies in the current NPC representation
                if unique_name in _asym_units_list:
                    # print(name)
                    # Get the asymmetric unit for that Nup
                    _asym = asym_unit_map[unique_name]
                    # Read in leaves
                    for leaf in IMP.atom.get_leaves(nup):
                        x, y, z = IMP.core.XYZR(leaf).get_coordinates()
                        resids = IMP.atom.Fragment(leaf).get_residue_indexes()
                        # feature in new system
                        # print('Bead- '+str(resids[0]-(chain_offset[name]-1))+':'+str(resids[-1]-(chain_offset[name]-1)))
                        new_feature=ihm.restraint.ResidueFeature([_asym(resids[0]-(chain_offset[name]-1),
                                                                        resids[-1]-(chain_offset[name]-1))])
                        # feature in old system
                        old_feature=ihm.restraint.PseudoSiteFeature(ihm.restraint.PseudoSite(x,y,z))
                        # Subcomplexes derived from 5ijo
                        subcomplexes_5ijo=['ir_core_1', 'ir_core_2', 'ir_core_3', 'ir_core_4',
                                           'ir_chan_1', 'ir_chan_2', 'ir_chan_3', 'ir_chan_4']
                        # Subcomplexes derived from 5a9q
                        subcomplexes_5a9q=['yc_inner_cr','yc_inner_nr','yc_outer_cr','yc_outer_nr','conn_1','conn_2']
                        count=0
                        # Find if subcomplex is in 5a9q
                        for subcomplex_5a9q in subcomplexes_5a9q:
                            if subcomplex_5a9q in unique_name:
                                dataset_choice=pdb_structure1
                                count+=1
                        # Find if subcomplex is in 5ijo
                        for subcomplex_5ijo in subcomplexes_5ijo:
                            if subcomplex_5ijo in unique_name:
                                dataset_choice=pdb_structure2
                                count += 1
                        # Check that exactly 1 subcomplex was found for each unique_name
                        if count!=1:
                            print('Warning. Incorrect subcomplex number found.')
                            print(str(count)+' numbers of subcomplex '+unique_name+' were found')
                        _pos_restraint_list.append(ihm.restraint.DerivedDistanceRestraint(dataset_choice,new_feature,
                            old_feature,ihm.restraint.HarmonicDistanceRestraint(0.0)))
    _pos_restraints=ihm.restraint.RestraintGroup(_pos_restraint_list)
    return _pos_restraints

# Function to find membrane bound portions of Nup155 and Nup160
def find_MBM_residues(_asymm_units,chain_offset):
    """Function to calculate membrane binding restraints
    """
    _objects=[]
    for _asymm_unit in _asymm_units:
        if "Nup155" in _asymm_unit.details:
            obj1=_asymm_unit(262-(chain_offset[name]-1),271-(chain_offset[name]-1))
            _objects.append(obj1)
        if "Nup160" in _asymm_unit.details:
            obj2=_asymm_unit(260-(chain_offset[name]-1),268-(chain_offset[name]-1))
            _objects.append(obj2)
    return ihm.restraint.ResidueFeature(_objects,details='Membrane bound residues on Nup155 and Nup160.')


# Output the final models as a time ordered sequence
# Includes only the centroid structure of each most
# likely state at each time point
# Add states to model
print('Adding ordered models to system...')
# empty restraint vector
restraint_list=[]
pos_restraint_list=[]
state_list=[]
for i in range(nstates):
    # generate model from hc
    m = MyModel(assembly=assemblies_list[i], protocol=protocol_list[i], representation=representation, name=best_states[times[i]])
    # Add restraints
    # EM restraint
    em_restraint=ihm.restraint.EM3DRestraint(processed_EM_list[i],assemblies_list[i],number_of_gaussians=150)
    em_restraint.fits[m]=ihm.restraint.EM3DRestraintFit(cross_correlation_coefficient=CCC[times[i]])
    restraint_list.append(em_restraint)
    # Membrane EV
    NE=ihm.geometry.HalfTorus(ihm.geometry.Center(0,0,0),(membrane_D[times[i]]+membrane_H[times[i]])/2,
        membrane_H[times[i]]/2,0,name='Nuclear envelope at time '+times[i]+'.')
    membrane_EV_restraint=ihm.restraint.OuterSurfaceGeometricRestraint(processed_EM_list[i],NE,
            ihm.restraint.ResidueFeature(assemblies_list[i]),
            ihm.restraint.HarmonicDistanceRestraint(0.0),harmonic_force_constant=0.01,restrain_all=True)
    restraint_list.append(membrane_EV_restraint)
    # Membrane attraction
    find_data=find_MBM_residues(assemblies_list[i],start_res)
    membrane_attraction_restraint=ihm.restraint.InnerSurfaceGeometricRestraint(processed_EM_list[i],NE,
        find_data,
        ihm.restraint.HarmonicDistanceRestraint(0.0),harmonic_force_constant=0.001,restrain_all=True)
    restraint_list.append(membrane_attraction_restraint)
    # Position restraint
    pos_restraint_list.append(find_pos_restraint(assemblies_list[i], start_res))
    # add model to model group
    model_group = ihm.model.ModelGroup([m], name="Centroid model at time "+times[i]+'.')
    # add ensemble
    ensemble = ihm.model.Ensemble(model_group,num_models=nstructures[times[i]],
        precision=structural_precision[times[i]],file=ihm.location.InputFileLocation(
        main_dir+'simulations_round2/Refined_energies_1model_460/filtered_noNup188/total/'
        +best_states[times[i]]+'/'+best_states[times[i]]+'_ensemble_v2.rmf',repo=None))
    # add state
    state = ihm.model.State([model_group])
    state_list.append(state)
    # add to system
    system.ensembles.append(ensemble)
system.state_groups.extend([ihm.model.StateGroup(state_list)])
# Add restraints to system
system.restraints.extend(restraint_list)
system.restraint_groups.extend(pos_restraint_list)
# Add ordered process
g = ihm.model.OrderedProcess('Time steps in NPC assembly.',
                            description='Time steps based on enumerating and scoring all '
                            'sufficiently good scoring possible trajectories of NPC assembly. '
                            'Step 1: 5min, Step 2: 6min, Step 3: 8min, Step 4: 10min '
                            'Step 5: 15min, Step 6: mature. '
                            'Representative structures are the centroid structure from '
                            'the most populated cluster, from RMSD clustering.')
# add connections between adjacent snapshots
step = ihm.model.ProcessStep(description='Assembly trajectory of the NPC.')
for i in range(len(system.state_groups[0])-1):
    step.append(ihm.model.ProcessEdge(system.state_groups[0][i][0],system.state_groups[0][i+1][0]))
g.steps.append(step)
system.ordered_processes.append(g)
print('Done.')

# Update local paths
print('Updating paths...')
repos = []
repos.append(ihm.location.Repository(
             doi="10.5281/zenodo.11129725", root="../",
             url="https://zenodo.org/record/11129725/files/PM_NPC_Assembly.zip",
             top_directory="PM_NPC_Assembly"))
system.update_locations_in_repositories(repos)
print('Done.')

# Write out in mmCIF
print('Writing mmCIF...')
with open("NPC_assembly.cif", "w") as fh:
    ihm.dumper.write(fh, [system])
print('Done.')
