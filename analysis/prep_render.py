"""Code to render a more visually appealing version of our NPC model. Renders the last rmf frame. Takes up to 2 inputs: #1. an RMF file
2nd: if True, will write an RMF for both individual subomplexes and the entire model
global imports"""
import IMP
import RMF
import IMP.rmf
import IMP.em
import IMP.core
import IMP.atom
import glob
import sys
sys.path.insert(0, "/wynton/home/sali/aplatham/NPC_Assembly/Tempkin_rerun/pm_assembly_code/src/analysis/sampling_analysis")


def write_rmf_for_chimera(model, hc, geom, write_subcomplex):
    """Write color coded set of RMF files for reading into Chimera.

    Creates a display hierarchy with just the protein particles containing all
    components as well as one RMF per color code.
    """

    sc_colors = {"yc_outer_cr": "green",
                 "yc_outer_nr": "green",
                 "yc_inner_cr": "green",
                 "yc_inner_nr": "green",
                 "ir_core_1": "cyan",  # contains 205, 93
                 "ir_core_2": "blue",  # contains 93
                 "ir_core_3": "cyan",  # contains 205, 93
                 "ir_core_4": "blue",  # contains 93
                 "ir_chan_1": "purple",
                 "ir_chan_2": "purple",
                 "ir_chan_3": "purple",
                 "ir_chan_4": "purple",
                 "conn_1": "orange",
                 "conn_2": "orange"}
    rgb_codes = {"orange": (0.133, 0.545, 0.133),
                 "green": (0, 0.447, 0.741),
                 "cyan": (0.85, 0.325, 0.098),
                 "blue": (0.494, 0.184, 0.556),
                 "purple": (0.635, 0.078, 0.184),
                 "yellow": (0.929, 0.694, 0.125)}


    # create a display structure to write to RMF
    hc_display_all = IMP.atom.Hierarchy(model, IMP.Particle(model))
    hc_display_all.set_name("npc")

    for sc in hc.get_children():
        # skip data densities
        if "Data_density" in sc.get_name():
            continue
        if "Sigma" in sc.get_name():
            continue
        if "bond" in sc.get_name():
            continue

        sc_display = IMP.atom.Hierarchy(model, IMP.Particle(model))
        sc_display.set_name(sc.get_name())
        hc_display_all.add_child(sc_display)

        # add only non-density nups
        for nup in sc.get_children():

            if "Density" in nup.get_name():
                continue

            nup_display = IMP.atom.Hierarchy(model, IMP.Particle(model))
            # remove density particles
            for p in nup.get_children():
                if IMP.core.Gaussian.get_is_setup(p):
                    continue
                nup_display.add_child(p)


            sc_display.add_child(nup_display)

    # collect hc_display for each colored subcomplex
    hc_display_sc = {}

    # color the representative display structure
    for sc in hc_display_all.get_children():
        # figure out sc label from name, ignore density
        k = None
        for key in sc_colors.keys():
            if key in sc.get_name():
                k = key
                break
        if k is None:
            break

        color = sc_colors[k]

        # if the color is in hc_display_sc, add it, if not create new entry
        if color not in hc_display_sc.keys():
            hc_sc = IMP.atom.Hierarchy(model, IMP.Particle(model))
            hc_sc.set_name(color + "_npc_subcomplexes")
            hc_display_sc[color] = hc_sc

        # add child to colored hierarchy
        hc_display_sc[color].add_child(sc)

        # look up rgb code
        rgb = rgb_codes[color]

        # color each particle
        for leaf in IMP.atom.get_leaves(sc):
            IMP.display.Colored.setup_particle(model, leaf, IMP.display.Color(rgb[0],rgb[1],rgb[2]))

    # write file wth all subcomplexes to disk
    rmf_out_fh = RMF.create_rmf_file("display_apl2.rmf")
    IMP.rmf.add_hierarchy(rmf_out_fh, hc_display_all)
    IMP.rmf.add_geometries(rmf_out_fh, geom)
    IMP.rmf.save_frame(rmf_out_fh, "0")

    # write file with just colored subcomplexes to disk
    if write_subcomplex:
        for key in hc_display_sc.keys():
            hc_sc = hc_display_sc[key]

            rmf_out_fh = RMF.create_rmf_file("display_" + str(key) + "_apl2.rmf")
            IMP.rmf.add_hierarchy(rmf_out_fh, hc_sc)
            IMP.rmf.add_geometries(rmf_out_fh, geom)
            IMP.rmf.save_frame(rmf_out_fh, "0")

    return

def process_traj(rmf_file,write_subcomplex=False):
    """Reads RMF file from disks to select lowest scoring state and write RMF and MRC files to disk.

    Parameters:
    ---------------------------
    oheader - directory to write output
    rmf_file - input rmf file to parse
    ref_mrc_file - MRC density to use as template for raster
    nframes - use just the last nframes
    rasterize_mrc - boolean for whether to rasterize subcomplex colored densities
    parallel_raster - select whether to call parallel raster routine
    """

    print('Running...')

    m = IMP.Model()

    # import the target file
    rmf_fh = RMF.open_rmf_file_read_only(rmf_file)

    # open hierarchy and remake scoring function
    npc = IMP.rmf.create_hierarchies(rmf_fh, m)[0]
    rs = IMP.rmf.create_restraints(rmf_fh, m)
    geom = IMP.rmf.create_geometries(rmf_fh)

    # get the last frame and write to a new file
    lowest_frame_id = rmf_fh.get_number_of_frames()-1
    IMP.rmf.load_frame(rmf_fh, RMF.FrameID(lowest_frame_id))

    # write new RMF file with lowest energy structure to disk
    write_rmf_for_chimera(m, npc, geom,write_subcomplex)
    print('Done.')

    return



RMF_filename=sys.argv[1]
if len(sys.argv)>2:
    if sys.argv[2]=='True' or sys.argv[2]=='TRUE' or sys.argv[2]=='true':
        process_traj(RMF_filename,True)
    else:
        process_traj(RMF_filename)
else:
    process_traj(RMF_filename)

