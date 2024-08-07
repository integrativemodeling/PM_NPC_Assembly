"""
Script to generate an mmCIF file suitable for deposition in PDB-dev.

Entities correspond to unique protein sequences.

Asymmetric units correspond to CG beads (in our model).

multiple AS can point to the same entity.

Each NPC subcomplex is represented by grouping it's constituent beads into an
ihm.Assembly grouping and the whole pore site is represented as ihm.Assembly of
subcomplex assemblies (whose relationship is denoted by the parent trait of
assemblies.
"""

import RMF
import IMP
import IMP.rmf
import sys
import ihm
import ihm.representation
import ihm.protocol
import ihm.dumper
import ihm.citations
import ihm.reference
import ihm.analysis
import random
import os
import re
from Bio import SeqIO
from io import StringIO
import requests

main_dir='../'

"""times=['5min','6min','8min','10min','15min','mature']
best_states={'5min':'2_5min','6min':'10_6min','8min':'14_8min',
             '10min':'4_10min','15min':'8_15min','mature':'1_mature'}"""
times=['mature']
best_states={'mature':'1_mature'}
nstates=len(times)

title = ("Integrative spatiotemporal modeling of biomolecular processes: "
         "application to the assembly of the Nuclear Pore Complex")
system = ihm.System(title=title)
# start IMP model
IMP_m = IMP.Model()

# Function to load hierarchy from RMF. By default, loads the last frame
def load_hc(best_state):
    # path to cluster centroid
    rmf_filename=main_dir+'simulations_round2/Refined_energies_1model_460/filtered_noNup188/total/'+best_state+'/cluster.0/cluster_center_model.rmf3'
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
          description="Program used to fit GMMs to 3DEM maps at each snapshot as well as fitting of forward model GMMs used in the restraints.",
          location='https://pdbj.org/gmfit/doc_gmconvert/README_gmconvert.html')
system.software.append(imp)
print('Done.')

# Include entities for each protein sequence
print('Adding entities...')
def build_entity_template(hc_tmpl, system):
    """Return an entity object for a Nup domain.
    """
    # Dictionary converting protein names to Uniprot entries
    Uniprot_dict={'Nup133':'Q8WUM0', 'Nup107':'P57740', 'Nup96':'P52948', 'SEC13':'P55735', 'SEH1':'Q96EE3', 'Nup85':'Q9BW27',
                  'Nup43':'Q8NFH3', 'Nup160':'Q12769', 'Nup37':'Q8NFH4', 'Nup93':'Q8N1F7', 'Nup205':'Q92621',
                  'Nup155':'O75694', 'Nup188':'Q5SRE5', 'p54':'Q7Z3B4', 'p58':'Q9BVL2', 'p62':'P37198'}
    # Loop over all proteins
    entities_dict = {}
    for subcomplex in hc_tmpl.get_children():
        for template in subcomplex.get_children():
            name = template.get_name().split("_")[0]
            if name not in entities_dict.keys():

                # Get sequence for correct uniprot entity depending on the name
                """
                #ref = ihm.reference.UniProtSequence.from_accession(Uniprot_dict[name])
                URL="http://www.uniprot.org/uniprot/"+Uniprot_dict[name]+".fasta"
                response = requests.post(URL)
                cData = ''.join(response.text)
                Seq = StringIO(cData)
                seq_dat = list(SeqIO.parse(Seq, 'fasta'))
                sequence=seq_dat[0].seq
                entity = ihm.Entity(sequence, description="_".join([name, subcomplex.get_name()]),references=[ref])
                """

                # Add empty sequence to save time. For debugging
                AAs=['R','H','K','D','E','S','T','N','Q','C','G','P','A','V','I','L','M','F','Y','W']
                sequence=''
                for i in range(10000):
                    AA_choice=random.randint(0, len(AAs)-1)
                    sequence += AAs[AA_choice]
                entity = ihm.Entity(sequence, description=name)

                entities_dict[name] = entity
                system.entities.append(entity)

    return entities_dict
# Load mature hierarchy
hc_mature=load_hc(best_states[times[-1]])
# Set entities from mature hierarchy
possible_entities=build_entity_template(hc_mature,system)
print('Done.')

# Define asymeteric units for each state
print('Building asymeteric units...')
def build_new_assembly_from_entities(hc_sc, edict, syst, asym_unit_map=None):
    """Return an assembly instance pointing to the provided entity dictionary template.
    """

    nups = [child for child in hc_sc.get_children() if "Density" not in child.get_name()]

    asym_units = []
    sc_representation = []
    for nup in nups:
        name = nup.get_name().split("_")[0]
        unique_name = "_".join([name, hc_sc.get_name()])
        asu = ihm.AsymUnit(edict[name], details=unique_name)
        nres = len(nup.get_children())
        rep = ihm.representation.FeatureSegment(asu,
                                                rigid=True,
                                                primitive="sphere",
                                                count=nres)
        sc_representation.append(rep)
        asym_units.append(asu)

        if asym_unit_map is not None:
            asym_unit_map[unique_name] = asu

    syst.asym_units.extend(asym_units)
    sc_assembly = ihm.Assembly(asym_units, description=hc_sc.get_name())

    return sc_assembly, sc_representation, asym_units

# Variables to save
representation_list=[]
assemblies_list=[]
asym_unit_map={}
# Loop over all times
for i in range(nstates):
    hc_temp=load_hc(best_states[times[i]])
    # variables for output
    assemblies=[]
    representations=[]
    asym=[]
    # For each subcomplex
    for sc in hc_temp.get_children():
        sc_assemble, sc_rep, sc_asym = build_new_assembly_from_entities(sc, possible_entities, system, asym_unit_map=asym_unit_map)
        # Add each subcomplex to the appropriate variables
        assemblies.extend(sc_assemble)
        representations.extend(sc_rep)
        asym.extend(sc_asym)
    # Add variables to their respective lists
    assemblies_list.append(ihm.Assembly(asym,name='Assembly at '+times[i]))
    representation_list.append(ihm.representation.Representation(representations,name='Representation at '+times[i]))
print('Done')



# Read in experimental datasets. Make list with different data at each time point,
print('Adding experimental data...')
# as ET data changes as a function of time
exp_data=[]
processed_EM_list=[]
# Original EM
em_database_loc = ihm.location.EMPIARLocation("EMD-3820")
raw_em = ihm.dataset.EMDensityDataset(em_database_loc,details='Original electron tomography dataset.')
# Original FCS
FCS_database_loc=ihm.location.DatabaseLocation("idr0115",details='Raw FCS data available on the Image Data Resource: https://idr.openmicroscopy.org/')
raw_FCS = ihm.dataset.Dataset(FCS_database_loc,details='Original fluorescence correlation spectroscopy dataset.')
# Processed FCS
FCS_excel_loc=ihm.location.InputFileLocation(main_dir+"data/qfluor_data/Exp_data.txt",repo=None)
processed_FCS = ihm.dataset.Dataset(FCS_excel_loc,details='Processed fluorescence correlation spectroscopy dataset.')
# PDB
pdb1_loc=ihm.location.PDBLocation("5a9q")
pdb_structure1=ihm.dataset.PDBDataset(pdb1_loc,details='PDB for the Y-complex and connecting complex')
pdb2_loc=ihm.location.PDBLocation("5ijo")
pdb_structure2=ihm.dataset.PDBDataset(pdb1_loc,details='PDB for the inner ring')
# Processed EM
for i in range(0,len(times)):
    em_processed_loc=ihm.location.InputFileLocation(main_dir+'data/fit_etdata_with_gmm/Andrew_run/data/'+times[i]+'_150.mrc',repo=None)
    processed_EM=ihm.dataset.EMDensityDataset(em_processed_loc,details='Processed EM data at '+times[i])
    processed_EM_list.append(processed_EM)
    # append all experimental data
    exp_data.append(ihm.dataset.DatasetGroup((raw_em,raw_FCS,processed_FCS,pdb_structure1,pdb_structure2,processed_EM),
                                              name='Snapshot model data for '+times[i]))
print('Done.')

# Add restraints to the model
print('Adding restraints...')
em_restraint_list=[]
for i in range(nstates):
    em_restraint=ihm.restraint.EM3DRestraint(processed_EM_list[i],assemblies_list[i],number_of_gaussians=150)
    em_restraint_list.append(em_restraint)
system.restraints.extend(em_restraint_list)
print('Done.')


# Define sampling of single snapshot model
print('Adding protocol...')
def snapshot_model_protocol(t,exp,assembly):
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
    return prot

# Create list for protocol of each snapshot
protocol_list=[]
# Loop over each state and determine protocol for producing that state
for i in range(nstates):
    protocol_list.append(snapshot_model_protocol(times[i],exp_data[i],assemblies_list[i]))
# Skip protocol for scoring graphs. Included in output models
"""# Protocols for graph modeling
graph_prot=ihm.protocol.Protocol('Spatiotemporal modeling')
# Build graphs of snapshot models
graph_prot.steps.append(ihm.protocol.Step(
                 assembly=assembly, dataset_group=All_data,
                 method='Graph modeling',
                 name='Enumerate and score graphs of trajectories of NPC assembly.',
                 num_models_begin=1782000, num_models_end=5184,
                 script_file=main_dir+'score_graph/score_trj_noNup188.py',
                 ordered=True,software=imp))
# Validation of graph
graph_prot.steps.append(ihm.analysis.ValidationStep(
                 feature='other', dataset_group=FCS_data_validate,
                 details='Validate model against FCS data for Nup188 and Seh1,'
                         ' which are not included in model scoring.',
                 script_file=main_dir+'analysis/analyze_copy_number.py',
                 software=imp))"""
print('Done.')

# Define class to extract spheres
class MyModel(ihm.model.Model):
    # override the get_spheres subclass to pull coordinates from rmf file on disk
    def get_spheres(self):
        # for each subcomplex
        for nups in hc.get_children():
            # for each nup
            for nup in nups.get_children():
                name = nup.get_name().split("_")[0]
                unique_name = "_".join([name, nups.get_name()])
                # Get the asymmeteric unit for that Nup
                _asym = asym_unit_map[unique_name]
                # Read in leaves
                for leaf in IMP.atom.get_leaves(nup):
                    # parse info for sphere
                    x, y, z = IMP.core.XYZR(leaf).get_coordinates()
                    resids = IMP.atom.Fragment(leaf).get_residue_indexes()
                    R = IMP.core.XYZR(leaf).get_radius()
                    yield ihm.model.Sphere(asym_unit=_asym,
                                           seq_id_range = (resids[0],resids[-1]),
                                           x=x,
                                           y=y,
                                           z=z,
                                           radius=R)

# Output the final models as a time ordered sequence
# Includes only the centroid structure of each most
# likely state at each time point
# Add states to model
print('Adding ordered models to system...')
"""for i in range(nstates):
    hc=load_hc(best_states[times[i]])
    m = MyModel(assembly=assemblies_list[i], protocol=protocol_list[i], representation=representation_list[i])
    model_group = ihm.model.ModelGroup([m], name="Centroid model at time "+times[i]+'.')
    state = ihm.model.State([model_group])
    system.state_groups.append(ihm.model.StateGroup([state]))
g = ihm.model.OrderedProcess('Time steps in NPC assembly.',
                            description='Time steps based on enumerating and scoring all '
                            'sufficiently good scoring possible trajectories of NPC assembly. '
                            'Step 1: 5min, Step 2: 6min, Step 3: 8min, Step 4: 10min '
                            'Step 5: 15min, Step 6: mature. '
                            'Representative structures are the centroid structure from '
                            'the most populated cluster, from RMSD clustering.')
system.ordered_processes.append(g)
step = ihm.model.ProcessStep(description='Assembly trajectory of the NPC.')
for i in range(len(system.state_groups[0])-1):
    step.append(ihm.model.ProcessEdge(system.state_groups[0][i],system.state_groups[0][i+1]))
g.steps.append(step)"""
i=0
print(system.restraints)
hc=load_hc(best_states[times[i]])
m = MyModel(assembly=assemblies_list[i], protocol=protocol_list[i], representation=representation_list[i])
model_group = ihm.model.ModelGroup([m], name="Centroid model at time "+times[0]+'.')
ensemble = ihm.model.Ensemble(model_group,num_models=801,precision=104.258,file=ihm.location.InputFileLocation(
    main_dir+'simulations_round2/Refined_energies_1model_460/filtered_noNup188/total/'
    +best_states[times[i]]+'/'+best_states[times[i]]+'_ensemble_v2.rmf',repo=None))
state = ihm.model.State([model_group])
system.ensembles.append(ensemble)
system.state_groups.append(ihm.model.StateGroup([state]))
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
