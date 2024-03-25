"""
Script that converts quasi-PDB file from gmconvert to gmm file used by IMP.isd
"""

import sys


def convert_pdb_to_isd(pdbfile, isdfile):
    """Write PDB file with GMM information (see gmconvert PDB output format)
    into ISD format in give file.
    """

    # scale parameter used for units conversion
    scale = 1.0  # convert from nm to ang

    with open(pdbfile, "r") as f_handle:
        lines = f_handle.readlines()
        lnum = []
        for num, line in enumerate(lines):
            if "HETATM" in line:
                lnum.append(num)

        natoms = len(lnum)

        # write out lines 
        with open(isdfile, "w") as f_out:

            f_out.write("#|num|weight|mean|covariance matrix|\n")

            for indx,l in enumerate(lnum):
                assert lines[l][0:6] == "HETATM"

                weight = float(lines[l+1].split()[-1])*scale
                mean = (str(float(lines[l+4].split()[-3])*scale),
                        str(float(lines[l+4].split()[-2])*scale),
                        str(float(lines[l+4].split()[-1])*scale))

                # extract covariance terms
                c_d = {}
                lcov1 = lines[l+5].split()
                lcov2 = lines[l+6].split()
                c_d["xx"] = float(lcov1[lcov1.index("xx") + 1])*scale**2
                c_d["xy"] = float(lcov1[lcov1.index("xy") + 1])*scale**2
                c_d["xz"] = float(lcov1[lcov1.index("xz") + 1])*scale**2
                c_d["yy"] = float(lcov2[lcov2.index("yy") + 1])*scale**2
                c_d["yz"] = float(lcov2[lcov2.index("yz") + 1])*scale**2
                c_d["zz"] = float(lcov2[lcov2.index("zz") + 1])*scale**2

                cov = (str(c_d["xx"]), # xx
                       str(c_d["xy"]), # xy
                       str(c_d["xz"]), # xz
                       str(c_d["xy"]), # yx
                       str(c_d["yy"]), # yy
                       str(c_d["yz"]), # yz
                       str(c_d["xz"]), # zx
                       str(c_d["yz"]), # zy
                       str(c_d["zz"])) # zz

                # write outfile
                f_out.write("|{0}|{1}|{2}|{3}|\n".format(indx,
                                                         weight,
                                                         " ".join(mean),
                                                         " ".join(cov)))


# convert output ISD file to map
"""
m = IMP.Model()
ps = []
gmm_util.decorate_gmm_from_text(sys.argv[2], ps, m)

gmm_util.write_gmm_to_map(ps, "50_scaled.mrc", 7.5, fast=True)
"""

if __name__ == "__main__":
    convert_pdb_to_isd(sys.argv[1], sys.argv[2])
