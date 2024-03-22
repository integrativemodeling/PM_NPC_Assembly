"""
Implementations of score terms used in the structural sampling of NPC intermediates.
"""


from gmm_util import decorate_gmm_from_text
import IMP.core
import IMP.algebra
import IMP.container
import IMP.atom
import IMP.isd
import IMP.npc
import sys
import numpy as np


def Overal_pos(m, ps, rs, sigma, tol, name="overall_pos_restraints"):
    """
    Function to restrain Nups to their original position
    m- IMP Model
    ps- heirarchy with only protein atoms
    rs- scoring function, includes all scores
    sigma- strength of OverallPositionRestraint term. Note: larger sigma, weaker restraint
    tol- tolerance of restraint. Restraint is 0 within this window.
    """

    # Get number of atoms and put their indicies in a list
    N=len(IMP.atom.get_leaves(ps))
    scnt = IMP.container.ListSingletonContainer(m, IMP.atom.get_leaves(ps))
    beads=scnt.get_contents()


    count_xyz=0
    # Loop over all atoms and restrain them by their position
    for i in range(0,N):
        i1 = beads[i]
        p1 = m.get_particle(i1)

        if (IMP.core.Gaussian.get_is_setup(p1)):
            continue

        sc = IMP.container.ListSingletonContainer(m, [p1], "atom")
        xyz = IMP.core.XYZ(p1).get_coordinates()
        pos = IMP.npc.OverallPositionRestraint(m, sc, xyz[0], xyz[1], xyz[2], tol, False, sigma)
        rs.append(pos)
        count_xyz=count_xyz+1

    print("created ", count_xyz, " position restraints.")

    return rs

def Z_restraint(m, ps, Z_lower, Z_upper, sigma, name="Z_pos_restraint"):
    """
    Function to restrain all beads to be below a given Z-position
    m- IMP Model
    ps- heirarchy with only protein atoms
    rs- scoring function, includes all scores
    sigma- strength of OverallPositionRestraint term. Note: larger sigma, weaker restraint
    tol- tolerance of restraint. Restraint is 0 within this window.
    """

    # Get number of atoms and put their indicies in a list
    scnt = IMP.container.ListSingletonContainer(m, IMP.atom.get_leaves(ps))
    pos = IMP.npc.ZAxialPositionRestraint(m, scnt, Z_lower, Z_upper, False, sigma)

    return pos


def XYRadial_restraint(m, ps, rs, sigma, tol, name="overall_pos_restraints"):
    """
    Function to restrain Nups to based on their original position
    This function specifically restrains the nuclear side Y-complex to be the same distance from the Z-axis
    m- IMP Model
    ps- heirarchy with only protein atoms
    rs- scoring function, includes all scores
    sigma- strength of OverallPositionRestraint term. Note: larger sigma, weaker restraint
    tol- tolerance of restraint. Restraint is 0 within this window.
    """

    # Get number of atoms and put their indicies in a list
    N = len(IMP.atom.get_leaves(ps))
    scnt = IMP.container.ListSingletonContainer(m, IMP.atom.get_leaves(ps))
    beads = scnt.get_contents()

    count_xyz = 0
    # Loop over all atoms and restrain them by their position
    for i in range(0, N):
        i1 = beads[i]
        p1 = m.get_particle(i1)
        if "yc_outer_nr" in IMP.atom.Hierarchy(p1).get_parent().get_parent().get_name() or "yc_inner_nr" in IMP.atom.Hierarchy(p1).get_parent().get_parent().get_name():
            if (IMP.core.Gaussian.get_is_setup(p1)):
                continue

            sc = IMP.container.ListSingletonContainer(m, [p1], "atom")
            xyz = IMP.core.XYZ(p1).get_coordinates()
            xy=np.sqrt(xyz[0]*xyz[0]+xyz[1]*xyz[1])
            pos = IMP.npc.XYRadialPositionRestraint(m, sc, xy-tol/2, xy+tol/2, False, sigma)

            rs.append(pos)
            count_xyz = count_xyz + 1

    print("created ", count_xyz, " position restraints.")

    return rs

def XYAxial_restraint(m, ps, rs, sigma, tol, name="overall_pos_restraints"):
    """
    Function to restrain Nups to their original position
    m- IMP Model
    ps- heirarchy with only protein atoms
    rs- scoring function, includes all scores
    sigma- strength of OverallPositionRestraint term. Note: larger sigma, weaker restraint
    tol- tolerance of restraint. Restraint is 0 within this window.
    """

    # Get number of atoms and put their indicies in a list
    N = len(IMP.atom.get_leaves(ps))
    scnt = IMP.container.ListSingletonContainer(m, IMP.atom.get_leaves(ps))
    beads = scnt.get_contents()

    count_xyz = 0
    # Loop over all atoms and restrain them by their position
    for i in range(0, N):
        i1 = beads[i]
        p1 = m.get_particle(i1)
        if "yc_outer_nr" in IMP.atom.Hierarchy(p1).get_parent().get_parent().get_name() or "yc_inner_nr" in IMP.atom.Hierarchy(p1).get_parent().get_parent().get_name():
            if (IMP.core.Gaussian.get_is_setup(p1)):
                continue

            sc = IMP.container.ListSingletonContainer(m, [p1], "atom")
            xyz = IMP.core.XYZ(p1).get_coordinates()
            x_pos = IMP.npc.XAxialPositionRestraint(m, sc, xyz[0]-tol/2, xyz[0]+tol/2, False, sigma)
            y_pos = IMP.npc.YAxialPositionRestraint(m, sc, xyz[1]-tol/2, xyz[1]+tol/2, False, sigma)

            rs.append(x_pos)
            rs.append(y_pos)
            count_xyz = count_xyz + 1

    print("created ", count_xyz, " position restraints.")

    return rs

def go_harmonic_bonds_2(m, ps, k, cutoff, name="go_restraints"):
    """
    Apply go terms by instantiating IMP:atom::Bonds and a BongSingletonScore.

    Ignores non-protein atoms.
    """
    bonds = []
    nbonds = 0

    scnt = IMP.container.ListSingletonContainer(m, IMP.atom.get_leaves(ps))
#    cpc = IMP.container.ClosePairContainer(scnt, cutoff, IMP.core.XYZClosePairsFinder(), 0.0)
    cpc = IMP.container.ClosePairContainer(scnt, cutoff, 0.0)

    m.update()

    pairs = cpc.get_contents()
    # add pairs within the cutoff.
    # not using closepairs container since we want the XYZ interpretation of distance
    for pair in pairs:
        i1, i2 = pair
        p1 = m.get_particle(i1)
        p2 = m.get_particle(i2)
        if IMP.atom.Hierarchy(p1).get_parent().get_parent().get_name() == IMP.atom.Hierarchy(p2).get_parent().get_parent().get_name():
            continue
        if (IMP.core.Gaussian.get_is_setup(p1)) or (IMP.core.Gaussian.get_is_setup(p2)):
            continue

        dist = IMP.core.get_distance(IMP.core.XYZ(p1), IMP.core.XYZ(p2))
        # only include distances that are less than the cutoff. Note that this is the direct particle-particle distance,
        # while ClosePairContainer includes excluded volume, cutting off less atoms
        if dist>cutoff:
            continue
        # decorate each as Bonded
        if not IMP.atom.Bonded.get_is_setup(p1):
            b1 = IMP.atom.Bonded.setup_particle(p1)
        if not IMP.atom.Bonded.get_is_setup(p2):
            b2 = IMP.atom.Bonded.setup_particle(p2)
        # create Bond particle
        bond = IMP.atom.create_custom_bond(IMP.atom.Bonded(p1),
                                           IMP.atom.Bonded(p2),
                                           dist,
                                           k)
        bonds.append(bond)
        nbonds += 1

    # create restraint
    bnds = IMP.container.ListSingletonContainer(m, bonds, "bonds")
    uf = IMP.core.Harmonic(0, 1)
    scf = IMP.atom.BondSingletonScore(uf)
    rs = IMP.container.SingletonsRestraint(scf, bnds, name)

    print("created ", nbonds, " bonds.")

    return bnds, rs



def go_truncated_2(m, ps, k, cutoff, cut2, name="go_restraints"):
    """
    Apply go terms by instantiating IMP:atom::Bonds and a BongSingletonScore.

    Uses a Gaussian to retrain atomic positions.
    m: model
    ps: portion of hierarchy to restrain (in this case, the Nup subcomplexes)
    k: depth of Gaussian well. Positive values are attractive, negative values are repulsive (so use positive k)
    cut2: where the truncated Gaussian begins to go to level off
    cutoff: cutoff for "bonded" Go pairs

    Ignores non-protein atoms. (those not included in ps)
    """
    bonds = []
    nbonds = 0

    scnt = IMP.container.ListSingletonContainer(m, IMP.atom.get_leaves(ps))
#    cpc = IMP.container.ClosePairContainer(scnt, cutoff, IMP.core.XYZClosePairsFinder(), 0.0)
    cpc = IMP.container.ClosePairContainer(scnt, cutoff, 0.0)

    m.update()

    pairs = cpc.get_contents()
    # add pairs within the cutoff.
    # not using closepairs container since we want the XYZ interpretation of distance
    for pair in pairs:
        i1, i2 = pair
        p1 = m.get_particle(i1)
        p2 = m.get_particle(i2)
        if IMP.atom.Hierarchy(p1).get_parent().get_parent().get_name() == IMP.atom.Hierarchy(p2).get_parent().get_parent().get_name():
            continue
        if (IMP.core.Gaussian.get_is_setup(p1)) or (IMP.core.Gaussian.get_is_setup(p2)):
            continue

        dist = IMP.core.get_distance(IMP.core.XYZ(p1), IMP.core.XYZ(p2))
        # only include distances that are less than the cutoff. Note that this is the direct particle-particle distance,
        # while ClosePairContainer includes excluded volume, cutting off less atoms
        if dist > cutoff:
            continue
        # decorate each as Bonded
        if not IMP.atom.Bonded.get_is_setup(p1):
            b1 = IMP.atom.Bonded.setup_particle(p1)
        if not IMP.atom.Bonded.get_is_setup(p2):
            b2 = IMP.atom.Bonded.setup_particle(p2)
        # create Bond particle
        bond = IMP.atom.create_custom_bond(IMP.atom.Bonded(p1),
                                           IMP.atom.Bonded(p2),
                                           dist)
        bonds.append(bond)
        nbonds += 1

    # create restraint
    bnds = IMP.container.ListSingletonContainer(m, bonds, "bonds")
    uf = IMP.core.TruncatedHarmonicBound(0, k, cut2)
    scf = IMP.atom.BondSingletonScore(uf)
    rs = IMP.container.SingletonsRestraint(scf, bnds, name)

    print("created ", nbonds, " bonds.")

    return bnds, rs

def go_gaussian_2(m, ps, k, sigma, cutoff, name="go_restraints"):
    """
    Apply go terms by instantiating IMP:atom::Bonds and a BongSingletonScore.

    Uses a Gaussian to retrain atomic positions.
    m: model
    ps: portion of hierarchy to restrain (in this case, the Nup subcomplexes)
    k: depth of Gaussian well. Positive values are attractive, negative values are repulsive (so use positive k)
    sigma: width of Gaussian well
    cutoff: cutoff for "bonded" Go pairs

    Ignores non-protein atoms. (those not included in ps)
    """
    import IMP.spb
    bonds = []
    nbonds = 0

    scnt = IMP.container.ListSingletonContainer(m, IMP.atom.get_leaves(ps))
#    cpc = IMP.container.ClosePairContainer(scnt, cutoff, IMP.core.XYZClosePairsFinder(), 0.0)
    cpc = IMP.container.ClosePairContainer(scnt, cutoff, IMP.core.QuadraticClosePairsFinder(), 0.0)

    m.update()

    pairs = cpc.get_contents()
    # add pairs within the cutoff.
    # not using closepairs container since we want the XYZ interpretation of distance
    for pair in pairs:
        i1, i2 = pair
        p1 = m.get_particle(i1)
        p2 = m.get_particle(i2)
        if IMP.atom.Hierarchy(p1).get_parent().get_parent().get_name() == IMP.atom.Hierarchy(p2).get_parent().get_parent().get_name():
            continue
        if (IMP.core.Gaussian.get_is_setup(p1)) or (IMP.core.Gaussian.get_is_setup(p2)):
            continue

        dist = IMP.core.get_distance(IMP.core.XYZ(p1), IMP.core.XYZ(p2))
        # only include distances that are less than the cutoff. Note that this is the direct particle-particle distance,
        # while ClosePairContainer includes excluded volume, cutting off less atoms
        if dist > cutoff:
            continue
        # decorate each as Bonded
        if not IMP.atom.Bonded.get_is_setup(p1):
            b1 = IMP.atom.Bonded.setup_particle(p1)
        if not IMP.atom.Bonded.get_is_setup(p2):
            b2 = IMP.atom.Bonded.setup_particle(p2)
        # create Bond particle
        bond = IMP.atom.create_custom_bond(IMP.atom.Bonded(p1),
                                           IMP.atom.Bonded(p2),
                                           dist)
        bonds.append(bond)
        nbonds += 1

    # create restraint
    bnds = IMP.container.ListSingletonContainer(m, bonds, "bonds")
    uf = IMP.spb.Gaussian(-k, 0.0, sigma)
    scf = IMP.atom.BondSingletonScore(uf)
    rs = IMP.container.SingletonsRestraint(scf, bnds, name)

    print("created ", nbonds, " bonds.")

    return bnds, rs


def sterics(ps, k, slack, name="excluded_volume"):
    """
    Instantiate sterics.
    """

    # add excluded volume terms between proteins
    ev = IMP.core.ExcludedVolumeRestraint(ps, k, slack, name)

    return ev


def membrane_go_restraints(m, ps, mem, k, zo, name="membrane go restraint"):
    """
    Instantiate membrane restraints using spherical cap representation.

    Parameters:
    --------------
    m - model
    ps - list of components to appyly restraints to
    k - force parameter for restraint
    zo - equilibrium depth of attachement
    """

    sps = IMP.npc.SlabWithSphericalIndentPairScore(zo, k)
    mbps = []
    for sc in ps:
        if "yc" in sc.get_name():

            nups = sc.get_children()
            for nup in nups:
                if "Nup160" in nup.get_name():
                    mbps.append(nup.get_child(2))  # nup160
                if "Nup133" in nup.get_name():
                    mbps.append(nup.get_child(0))  # nup133

    mlsc = IMP.container.ListSingletonContainer(m, [mem])
    plsc = IMP.container.ListSingletonContainer(m, mbps)
    bpc = IMP.container.AllBipartitePairContainer(mlsc, plsc)
    # NTS: CONSIDER IMP::container::create_restraint
    #prs = imp.container.PairsRestraint(sps, bpc, name)
    prs = IMP.container.create_restraint(sps, bpc)

    return prs


def membrane_indent_ev_restraints(m, mem, ps, k, name="membrane indent ev score"):
    """
    Instantiate EV restraints for membrane-protein interactions.
    """
    # add excluded volume restraints to membrane
    mlsc = IMP.container.ListSingletonContainer(m, [mem])
    sdc = IMP.npc.SphericalIndentSurfaceDepthPairScore(k)
    alist = IMP.container.ListSingletonContainer(m, IMP.atom.get_leaves(ps))
    bpc = IMP.container.AllBipartitePairContainer(mlsc, alist)
    # NTS: CONSIDER IMP::container::create_restraint
    sdpr = IMP.container.create_restraint(sdc, bpc)

    return sdpr


def em_scores(m, hc, target_density_file, sig_mass = 1.0, scaling_factor = 1.0, model_cutoff = 1000.0, data_cutoff = 1000.0, slope = 10.0, parallel=False):
    """
    Build EM scores.
    """
    # here we can use the IMP.isd.gmm_tools.decorate_gmm_from_text routine for the decoration of the 
    # particles. First we'll need to get from the quasi-pdb format to the ISD format (done)
    tar_ps = []
    data_density = IMP.atom.Hierarchy(m, IMP.Particle(m))
    data_density.set_name("Data_density")
    decorate_gmm_from_text(target_density_file, tar_ps, m)  # copied from gmm_tools 
    gmm_data = [IMP.core.Gaussian(p) for p in tar_ps]

    for indx,p in enumerate(gmm_data):
        p.set_name("Data_density_" + str(indx))
        p.set_coordinates_are_optimized(False)
        data_density.add_child(p)

    # fit each subcomplex with a GMM 
    gmm_model = []

    # get model density particles from structure
    for sc in hc.get_children():
        nups = sc.get_children()
        for nup in nups:
            if "Density" in nup.get_name():
                gmm_model += nup.get_children()

    # apply scaling factor to variances of model gmm
    for g in gmm_model:
        IMP.core.Gaussian(g).set_variances(IMP.core.Gaussian(g).get_variances()*scaling_factor)

    # set global sigma parameter
    sigma = IMP.core.XYZR.setup_particle(m, IMP.Particle(m))
    sigma.set_name("Sigma")
    sigma.set_radius(1.0)
    sigma.set_coordinates(IMP.algebra.Vector3D(0.0, 0.0, 0.0))
    msig = IMP.atom.Mass.setup_particle(m, sigma, sig_mass)

    if parallel:
        gmm_rst = IMP.isd.GaussianEMRestraintParallel(m,
                                                      gmm_data,
                                                      gmm_model,
                                                      msig,
                                                      model_cutoff,
                                                      data_cutoff,
                                                      slope)

    else:
        gmm_rst = IMP.isd.GaussianEMRestraint(m,
                                              gmm_data,
                                              gmm_model,
                                              msig,
                                              model_cutoff,
                                              data_cutoff,
                                              slope)

    gmm_rst.compute_initial_scores()

    return gmm_rst, data_density, sigma


def membrane_toroid_go_restraints(m, ps, mem, k, zo, name="membrane go restraint"):
    """
    Instantiate membrane restraints using spherical cap representation.

    Parameters:
    --------------
    m - model
    ps - list of components to appyly restraints to
    k - force parameter for restraint
    zo - equilibrium depth of attachement
    """

    sps = IMP.npc.SlabWithToroidalPoreGoPairScore(zo, k)

    mlsc = IMP.container.ListSingletonContainer(m, [mem])
    plsc = IMP.container.ListSingletonContainer(m, ps)
    bpc = IMP.container.AllBipartitePairContainer(mlsc, plsc)
    prs = IMP.container.PairsRestraint(sps, bpc, name)

    return prs


def membrane_toroid_mbm_restraints(m, ps, mem, k, zo_upper, zo_lower, name="membrane mbm restraint"):
    """
    Instantiate membrane restraints using spherical cap representation.

    Parameters:
    --------------
    m - model
    ps - list of components to appyly restraints to
    k - force parameter for restraint
    zo_upper - equilibrium of upper harmonic wall
    zo_lower - equilibrium of lower harmonic wall
    """

    sps = IMP.npc.SlabWithToroidalPoreMBMScore(zo_upper, zo_lower, k)

    mlsc = IMP.container.ListSingletonContainer(m, [mem])
    plsc = IMP.container.ListSingletonContainer(m, ps)
    bpc = IMP.container.AllBipartitePairContainer(mlsc, plsc)
    prs = IMP.container.PairsRestraint(sps, bpc, name)

    return prs


def membrane_indent_mbm_restraints(m, ps, mem, k, zo_upper, name="membrane mbm restraints"):
    """
    Instantiate membrane restraints for spherical cap representation.
    """
    sps = IMP.npc.SlabWithSphericalIndentMBMScore(zo_upper, k)

    mlsc = IMP.container.ListSingletonContainer(m, [mem])
    plsc = IMP.container.ListSingletonContainer(m, ps)
    bpc = IMP.container.AllBipartitePairContainer(mlsc, plsc)
    prs = IMP.container.PairsRestraint(sps, bpc, name)

    return prs


def membrane_toroid_ev_restraints(m, mem, ps, k, name="membrane ev score"):
    """
    Instantiate EV restraints for membrane-protein interactions.
    """
    # add excluded volume restraints to membrane
    mlsc = IMP.container.ListSingletonContainer(m, [mem])
    sdc = IMP.npc.ToroidalPoreSurfaceDepthPairScore(k)
    alist = IMP.container.ListSingletonContainer(m, IMP.atom.get_leaves(ps))
    bpc = IMP.container.AllBipartitePairContainer(mlsc, alist)
    sdpr = IMP.container.PairsRestraint(sdc, bpc, name)

    return sdpr

def initialize_restraints(m, hc, boundary, prm):
    """Initializes and returns a list of restraint functions.
    Reads parameters from input dictionary.
    """

    # handle for model components
    mdlc = prm["model_components"]

    # list of restraint terms
    rs = []


    subcomplexes = hc.get_children()

    # prepare a hierarchy without model density particles
    # this representation simplifies some 
    # restraint initlization calls but isn't required otherwise
    prot_sc = IMP.atom.Hierarchy(m, IMP.Particle(m))
    for sc in subcomplexes:
        phc = IMP.atom.Hierarchy(m, IMP.Particle(m))
        phc.set_name(sc.get_name())

        for nup in sc.get_children():
            #if not IMP.core.Gaussian.get_is_setup(nup):
            if "Density" not in nup.get_name():
                phc.add_child(nup)
        prot_sc.add_child(phc)

    # prepare a hierarchy with only yc_nr
    yc_nr = IMP.atom.Hierarchy(m, IMP.Particle(m))
    for sc in subcomplexes:
        phc = IMP.atom.Hierarchy(m, IMP.Particle(m))
        if "yc_outer_nr" in sc.get_name() or "yc_inner_nr" in sc.get_name():
            phc.set_name(sc.get_name())
            for nup in sc.get_children():
            # if not IMP.core.Gaussian.get_is_setup(nup):
                if "Density" not in nup.get_name():
                    phc.add_child(nup)
            yc_nr.add_child(phc)

    # prepare a hierarchy with only yc_cr
    yc_cr = IMP.atom.Hierarchy(m, IMP.Particle(m))
    for sc in subcomplexes:
        phc = IMP.atom.Hierarchy(m, IMP.Particle(m))
        if "yc_outer_cr" in sc.get_name() or "yc_inner_cr" in sc.get_name():
            phc.set_name(sc.get_name())
            for nup in sc.get_children():
                # if not IMP.core.Gaussian.get_is_setup(nup):
                if "Density" not in nup.get_name():
                    phc.add_child(nup)
            yc_cr.add_child(phc)


    # build boundary and score
    print("Building boundary box and boundary score...")

    bndp = prm["boundary"]
    bscr = IMP.core.BoundingBox3DSingletonScore(IMP.core.HarmonicUpperBound(bndp["boundary_zo"], bndp["boundary_k"]),
                                                boundary)
    brst = IMP.container.SingletonsRestraint(bscr,
                                            IMP.container.ListSingletonContainer(m, subcomplexes),
                                            "Boundary score")
    rs.append(brst)
    print(" Done.")

    # add go model terms
    if mdlc["go_model"]:
        gomp = prm["go_model"]
        print("Building Go-like scores...")

        if mdlc["useTruncated"]:
            bonds, rs_go = go_truncated_2(m, prot_sc, gomp["go_k"], gomp["go_cutoff"], gomp["trunc_cut"])

            # set the weight of the go-parameters to be 1 / number_of_bonds
            if gomp["weight"] == "recip":
                bondsweight = float(1.0 / len(bonds.get_contents()))
            else:
                bondsweight = float(gomp["weight"])
            print("setting native contacts restraint weights to", bondsweight)
            rs_go.set_weight(bondsweight)

            # add bonded term to restraints
            rs.append(rs_go)

        elif mdlc["useGaussian"]:
            bonds, rs_go = go_gaussian_2(m, prot_sc, gomp["go_k"], gomp["go_sigma"], gomp["go_cutoff"])

            # set the weight of the go-parameters to be 1 / number_of_bonds
            if gomp["weight"] == "recip":
                bondsweight = float(1.0 / len(bonds.get_contents()))
            else:
                bondsweight = float(gomp["weight"])
            print("setting native contacts restraint weights to", bondsweight)
            rs_go.set_weight(bondsweight)

            # add bonded term to restraints
            rs.append(rs_go)

        elif mdlc["useBonds"]:
            bonds, rs_go = go_harmonic_bonds_2(m, prot_sc, gomp["go_k"], gomp["go_cutoff"])
            #bonds, rs_go = scores.go_harmonic_bonds(m, subcomplexes, gomp["go_k"], gomp["go_cutoff"])

            # set the weight of the go-parameters to be 1 / number_of_bonds 
            if gomp["weight"] == "recip":
                bondsweight = float(1.0 / len(bonds.get_contents()))
            else:
                bondsweight = float(gomp["weight"])
            print("setting native contacts restraint weights to", bondsweight)
            rs_go.set_weight(bondsweight)

            # add bonded term to restraints
            rs.append(rs_go)

        else:
            print("Note: No go bonds used in this model")
            pass
            #rs_go = go_harmonic(m, subcomplexes, gomp["go_k"], gomp["go_cutoff"])

            # add go terms to restraints
            #for r in rs_go:
            #    rs.append(r)

        print(" Done.")


    # structure based membrane restraints
    if mdlc["membrane"]:
        print("Constructing membrane and membrane restraints... ")
        memp = prm["membrane"]
        """
        y-complex (5a9q)
        Nup155 -> betaProp: 58-520
                  MBM: 241-275
        Nup160 -> betaProp: 41-438

        ir (5ijo)
        Nup155 -> betaProp: 58-520

        From ignacia:
        hNup160 -> 260-268
        hNup133 -> 245-267
        hNup155 -> 262-271
        """

        assert memp["mem_type"] in ["toroid", "indent"], "Unknown membrane model requested."

        pmbm = []
        for sc in subcomplexes:
            for nup in sc.get_children():
                if "Nup155" in nup.get_name():
                    for frg in nup.get_children():
                        if bool(set(IMP.atom.Fragment(frg).get_residue_indexes()).intersection(range(262,271))):
                            pmbm.append(frg)

                if "Nup160" in nup.get_name():
                    for frg in nup.get_children():
                        if bool(set(IMP.atom.Fragment(frg).get_residue_indexes()).intersection(range(260,268))):
                            pmbm.append(frg)

        if memp["mem_type"] == "toroid":
            major_r = (memp["R"] + memp["h"])/2.0
            membrane = IMP.npc.SlabWithToroidalPore.setup_particle(m, IMP.Particle(m), memp["h"], major_r)
            #mem_geom = IMP.npc.SlabWithToroidalPoreWireGeometry(memp["h"], major_r, bndp["L"])

            rs_mbm_mem = membrane_toroid_mbm_restraints(m, pmbm, membrane, memp["k_mem"], memp["zo_upper"], memp["zo_lower"])

            # rs_mbm_mem.set_weight(float(1.0 / len(pmbm)))
            # set the weight of the membrane binding harmonics term to 1.0
            rs_mbm_mem.set_weight(float(1.0))

            rs_ev_mem = membrane_toroid_ev_restraints(m, membrane, subcomplexes, memp["ev_mem"])

            rs.append(rs_mbm_mem)
            rs.append(rs_ev_mem)

        elif memp["mem_type"] == "indent":
            membrane = IMP.npc.SlabWithSphericalIndent.setup_particle(m, IMP.Particle(m), memp["R"], memp["h"])
            #mem_geom = IMP.npc.SlabWithSphericalIndentGeometry(memp["R"], bndp["L"], memp["h"])

            rs_mbm_mem = membrane_indent_mbm_restraints(m, pmbm, membrane, memp["k_mem"], memp["zo_mem"])

            # rs_mbm_mem.set_weight(float(1.0 / len(pmbm)))
            # set the weight of the membrane binding harmonics term to 1.0
            rs_mbm_mem.set_weight(float(1.0))
            rs_ev_mem = membrane_indent_ev_restraints(m, membrane, subcomplexes, memp["ev_mem"])

            rs.append(rs_mbm_mem)
            rs.append(rs_ev_mem)

        else:
            print("Error: Not recognized membrane type. Membrane not initialized")

        print(" Done.")

    # structure based position restraints --------------------------------------------------------------------------
    print('Building position based scores...')

    pos_prm = prm["pos_score"]
    # OverallPos - restrains Nups to their original position
    if mdlc["OverallPos"]:
        rs=Overal_pos(m, prot_sc, rs, pos_prm["pos_sigma"], pos_prm["pos_tol"], name="overall_pos_restraints")

    # XYRadial - restrains selected Nups (nr y-complex) to the radial distance from the Z-axis
    elif mdlc["XYRadial"]:
        rs=XYRadial_restraint(m, prot_sc, rs, pos_prm["pos_sigma"], pos_prm["pos_tol"], name="overall_pos_restraints")

    # XYAxial - adds seperate restrains to selected Nups (nr y-complex) on the X position and Y position
    elif mdlc["XYAxial"]:
        rs=XYAxial_restraint(m, prot_sc, rs, pos_prm["pos_sigma"], pos_prm["pos_tol"], name="overall_pos_restraints")

    else:
        print('Note: no position restraints used')
    print('Done')

    # Z-position restraint -----------------------------------------------------------------------------------------
    if mdlc["Zscore"]:
        Z_prm = prm["Zscore"]
        Z_rs1 = Z_restraint(m, yc_nr, Z_prm["Z_lower"], Z_prm["Z_mid"], Z_prm["Z_sigma"], name="Zmax_restraint_nr")
        Z_rs2 = Z_restraint(m, yc_cr, Z_prm["Z_mid"], Z_prm["Z_upper"], Z_prm["Z_sigma"], name="Zmax_restraint_cr")
        rs.append(Z_rs1)
        rs.append(Z_rs2)


    # ------------ EM DENSITY RESTRAINTS --------------------
    if mdlc["gmm_restraint"]:
        print("Initializing Gaussian EM Restraint...")
        gmmp = prm["em_score"]
        gmm_rest, data_density, msig = em_scores(m,
                                                 hc,
                                                 gmmp["gmm_target_density"],
                                                 gmmp["gmm_sig_mass"],
                                                 gmmp["gmm_scaling_factor"],
                                                 gmmp["gmm_model_cutoff"],
                                                 gmmp["gmm_data_cutoff"],
                                                 gmmp["gmm_slope"],
                                                 gmmp["omp_parallel"])
        gmm_rest.set_weight(gmmp["gmm_weight"])
        gmm_rest.compute_initial_scores()

        rs.append(gmm_rest)

        hc.add_child(data_density)
        hc.add_child(msig)
        print(" Done.")

    # ---- EXCLUDED VOLUME ------
    # add excluded volume terms between proteins
    if mdlc["sterics"]:
        gomp = prm["go_model"]
        print("Building sterics...")
        ev = sterics(prot_sc.get_children(), gomp["sterics_k"], gomp["sterics_slack"])
        rs.append(ev)
        print(" Done.")

    return rs

