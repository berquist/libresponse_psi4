#
# @BEGIN LICENSE
#
# libresponse_psi4 by Eric Berquist, a plugin to:
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2019 The Psi4 Developers.
#
# The copyrights for code used from other parties are included in
# the corresponding files.
#
# This file is part of Psi4.
#
# Psi4 is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, version 3.
#
# Psi4 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License along
# with Psi4; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
# @END LICENSE
#

import psi4
import psi4.driver.p4util as p4util
from psi4.driver.procrouting import proc_util


def disable_symmetry(molecule):
    """Return a molecule with symmetry completely disabled."""
    molecule.update_geometry()
    if molecule.schoenflies_symbol() != "c1":
        psi4.core.print_out(
            """  A requested method does not make use of molecular symmetry: """
            """further calculations in C1 point group.\n"""
        )
        molecule = molecule.clone()
        molecule.reset_point_group("c1")
        # TODO the orientation and absolute position (COM translation)
        # has already been messed with at this point! Need to disable
        # in the input file.
        molecule.fix_orientation(True)
        molecule.fix_com(True)
        molecule.update_geometry()
        # psi4_string = molecule.create_psi4_string_from_molecule()
        # print(psi4_string)
    return molecule


def run_libresponse_psi4(name, **kwargs):
    r"""Function encoding sequence of PSI module and plugin calls so that
    libresponse_psi4 can be called via :py:func:`~driver.energy`. For post-scf
    plugins.

    >>> energy('libresponse_psi4')
    """

    kwargs = p4util.kwargs_lower(kwargs)

    # Your plugin's psi4 run sequence goes here
    psi4.core.set_local_option("MYPLUGIN", "PRINT", 1)

    # The response density is asymmetric; need to build generalized
    # J/K.
    proc_util.check_non_symmetric_jk_density("libresponse")

    # Disable symmetry, even for passed-in wavefunctions.
    molecule = kwargs.get("molecule", None)
    if molecule is None:
        molecule = psi4.core.get_active_molecule()
    molecule = disable_symmetry(molecule)
    kwargs["molecule"] = molecule

    # Compute a SCF reference, a wavefunction is return which holds the
    # molecule used, orbitals Fock matrices, and more
    ref_wfn = kwargs.get("ref_wfn", None)
    if ref_wfn is None:
        ref_wfn = psi4.driver.scf_helper(name, **kwargs)

    # Ensure IWL files have been written when not using DF/CD
    proc_util.check_iwl_file_from_scf_type(
        psi4.core.get_option("SCF", "SCF_TYPE"), ref_wfn
    )

    # Call the Psi4 plugin
    # Please note that setting the reference wavefunction in this way is ONLY for plugins
    libresponse_psi4_wfn = psi4.core.plugin("libresponse_psi4.so", ref_wfn)

    return libresponse_psi4_wfn


# Integration with driver routines
psi4.driver.procedures["energy"]["libresponse_psi4"] = run_libresponse_psi4
