/*
 * @BEGIN LICENSE
 *
 * libresponse_psi4 by Eric Berquist, a plugin to:
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2019 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This file is part of Psi4.
 *
 * Psi4 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * Psi4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with Psi4; if not, write to the Free Software Foundation, Inc., 51
 * Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

#include <armadillo>

#include <libresponse/libresponse.h>
#include <libresponse/linear/interface.h>

#include "psi4/liboptions/liboptions.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/psi4-dec.h"

#include "psi4/libmints/basisset.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/vector.h"
#include "psi4/libmints/wavefunction.h"
#include "psi4/libpsi4util/PsiOutStream.h"

#include "matvec_psi4.h"
// #include "operatordatatype.h"
#include "parse_operators.h"
#include "wrappers.h"

namespace psi {
namespace libresponse_psi4 {

extern "C" PSI_API int read_options(const std::string &name, Options &options) {
    if (name == "LIBRESPONSE_PSI4" || options.read_globals()) {
        options.add_str("order", "linear");
        options.add_str("solver", "linear");
        options.add_str("hamiltonian", "rpa");
        options.add_str("spin", "singlet");
        options.add_int("maxiter", 60);
        // void add(std::string key, DataType* data);
        // OperatorDataType etc.
        options.add_array("OPERATOR_DIPOLE");
        options.add_array("OPERATOR_QUADRUPOLE");
        // options.add_array("OPERATOR_MULTIPOLE");
        options.add_array("OPERATOR_NABLA");
        options.add_array("OPERATOR_ANGMOM");
    }

    return true;
}

extern "C" PSI_API SharedWavefunction libresponse_psi4(SharedWavefunction ref_wfn,
                                                       Options &options) {
    // This is the Psi4 output file. Replace the cout stream buffer
    // with the one from Psi4.
    std::ostream *ofs = outfile->stream();
    std::streambuf *cout_original_buff = std::cout.rdbuf();
    std::cout.rdbuf(ofs->rdbuf());

    options.print();

    const size_t NOrb = ref_wfn->nmo();
    const size_t NOa = ref_wfn->nalpha();
    const size_t NOb = ref_wfn->nbeta();
    const size_t NVa = NOrb - NOa;
    const size_t NVb = NOrb - NOb;

    const size_t NDen = ref_wfn->same_a_b_dens() ? 1 : 2;

    const std::shared_ptr<BasisSet> bset = ref_wfn->basisset();
    const size_t NBasis = bset->nbf();

    arma::cube C(NBasis, NOrb, NDen);
    arma::mat moene(NOrb, NDen);

    // Psi4 is C++; storage is row major!
    const SharedMatrix mCa = ref_wfn->Ca()->transpose();
    const SharedMatrix mCb = ref_wfn->Cb()->transpose();
    const SharedVector vEa = ref_wfn->epsilon_a();
    const SharedVector vEb = ref_wfn->epsilon_b();
    const arma::mat Ca(mCa->get_pointer(), mCa->rowdim(), mCa->coldim(), true);
    const arma::vec Ea(vEa->pointer(), vEa->dim(), true);
    C.slice(0) = Ca;
    moene.col(0) = Ea;
    if (NDen == 2) {
        const arma::mat Cb(mCb->get_pointer(), mCb->rowdim(), mCb->coldim(), true);
        const arma::vec Eb(vEb->pointer(), vEb->dim(), true);
        C.slice(1) = Cb;
        moene.col(1) = Eb;
    }

    arma::uvec occupations(4);
    occupations(0) = NOa;
    occupations(1) = NVa;
    occupations(2) = NOb;
    occupations(3) = NVb;

    // Configuration.
    libresponse::configurable libresponse_options;
    set_defaults(libresponse_options);
    // TODO input option
    libresponse_options.cfg<int>("print_level", 2);
    // Psi4 must take C_left/C_right, not Dg!
    libresponse_options.cfg<bool>("_do_compute_generalized_density", false);

    // Frequencies.
    std::vector<double> omega;
    // Force the static case for now.
    omega.push_back(0.0);

    // 1-electron integrals for operators.
    std::vector<libresponse::operator_spec> operators;
    parse_operators(ref_wfn, options, operators);

    // 2-electron integral engine.
    MatVec_Psi4 *matvec = new MatVec_Psi4(ref_wfn, options);

    libresponse::SolverIterator_linear *solver_iterator = new libresponse::SolverIterator_linear();

    arma::cube results;

    libresponse::solve_linear_response(results,
                                       matvec,
                                       solver_iterator,
                                       C,
                                       moene,
                                       occupations,
                                       omega,
                                       operators,
                                       libresponse_options);

    delete matvec;
    delete solver_iterator;

    // Final printing.
    std::ostringstream os;
    os.precision(6);
    os.setf(std::ios_base::fixed, std::ios_base::floatfield);

    // os << " " << dashes << std::endl;
    // print_polarizability(os, results.slice(0));
    // os << " " << rcarats << std::endl;

    std::cout << os.str();

    // Put the original cout stream buffer back.
    std::cout.rdbuf(cout_original_buff);

    // Typically you would build a new wavefunction and populate it with data
    return ref_wfn;
}

} // namespace libresponse_psi4
} // namespace psi
