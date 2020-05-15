#include "matvec_psi4.h"

#include "psi4/libmints/basisset.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libqt/qt.h"
#include "psi4/libscf_solver/hf.h"

#include <libresponse/utils.h>

MatVec_Psi4::MatVec_Psi4(SharedWavefunction wfn, Options &options) : options_(options) {

    wfn_ = wfn;
    std::shared_ptr<BasisSet> basisset_ = wfn_->basisset();

    // TODO memory is in doubles, not MB...
    const size_t memory_ = 100000;

    // see psi4/libfock/apps.cc/RBase::preiterations()
    if (!jk_) {
        if (options_.get_bool("SAVE_JK")) {
            jk_ = (static_cast<psi::scf::HF *>(wfn_.get()))->jk();
            outfile->Printf("    Reusing JK object from SCF.\n\n");
        } else {
            if (options_.get_str("SCF_TYPE") == "DF") {
                jk_ = JK::build_JK(basisset_, wfn->get_basisset("DF_BASIS_SCF"), options_);
            } else {
                jk_ = JK::build_JK(basisset_, BasisSet::zero_ao_basis_set(), options_);
            }
            size_t effective_memory =
                (size_t)(0.125 * options_.get_double("CPHF_MEM_SAFETY_FACTOR") * memory_);
            jk_->set_memory(effective_memory);
            jk_->initialize();
        }
    }
}

MatVec_Psi4::~MatVec_Psi4() { jk_->finalize(); }

void MatVec_Psi4::compute(arma::cube &J,
                          arma::cube &K,
                          const std::vector<arma::mat> &L,
                          const std::vector<arma::mat> &R) {

    if (J.n_rows != J.n_cols)
        throw std::runtime_error("J.n_rows != J.n_cols");
    if (K.n_rows != K.n_cols)
        throw std::runtime_error("K.n_rows != K.n_cols");
    if (J.n_slices != K.n_slices)
        throw std::runtime_error("J.n_slices != K.n_slices");
    if (L.size() != R.size())
        throw std::runtime_error("L.size() != R.size()");
    if (L.size() != K.n_slices)
        throw std::runtime_error("L.size() != K.n_slices");

    const size_t nden = L.size();

    if (R[0].n_rows != L[0].n_rows)
        throw std::runtime_error("R[0].n_rows != L[0].n_rows");
    if (R[0].n_cols != L[0].n_cols)
        throw std::runtime_error("R[0].n_cols != L[0].n_cols");
    if (nden == 2) {
        if (R[1].n_rows != L[1].n_rows)
            throw std::runtime_error("R[1].n_rows != L[1].n_rows");
        if (R[1].n_cols != L[1].n_cols)
            throw std::runtime_error("R[1].n_cols != L[1].n_cols");
    }

    const size_t nbasis = L[0].n_rows;

    std::vector<SharedMatrix> &vpL = jk_->C_left();
    std::vector<SharedMatrix> &vpR = jk_->C_right();
    vpL.clear();
    vpR.clear();

    const size_t nvirt_alph = L[0].n_cols;
    SharedMatrix pL_alph(new Matrix("L (alpha)", nbasis, nvirt_alph));
    SharedMatrix pR_alph(new Matrix("R (alpha)", nbasis, nvirt_alph));
    // TODO transpose without copy?
    arma::mat tmpL = L[0].t();
    arma::mat tmpR = R[0].t();
    C_DCOPY(nbasis * nvirt_alph, tmpL.memptr(), 1, pL_alph->pointer()[0], 1);
    C_DCOPY(nbasis * nvirt_alph, tmpR.memptr(), 1, pR_alph->pointer()[0], 1);
    vpL.push_back(pL_alph);
    vpR.push_back(pR_alph);
    if (nden == 2) {
        const size_t nvirt_beta = L[1].n_cols;
        SharedMatrix pL_beta(new Matrix("L (beta)", nbasis, nvirt_beta));
        SharedMatrix pR_beta(new Matrix("R (beta)", nbasis, nvirt_beta));
        tmpL = L[1].t();
        tmpR = R[1].t();
        C_DCOPY(nbasis * nvirt_beta, tmpL.memptr(), 1, pL_beta->pointer()[0], 1);
        C_DCOPY(nbasis * nvirt_beta, tmpR.memptr(), 1, pR_beta->pointer()[0], 1);
        vpL.push_back(pL_beta);
        vpR.push_back(pR_beta);
    }

    jk_->compute();

    std::vector<SharedMatrix> vpJ = jk_->J();
    std::vector<SharedMatrix> vpK = jk_->K();

    if (vpJ.size() != nden)
        throw std::runtime_error("vpJ.size() != nden");
    if (vpK.size() != nden)
        throw std::runtime_error("vpK.size() != nden");

    for (size_t s = 0; s < nden; s++) {
        arma::mat mJ(vpJ[s]->get_pointer(), vpJ[s]->rowdim(), vpJ[s]->coldim(), false, true);
        arma::mat mK(vpK[s]->get_pointer(), vpK[s]->rowdim(), vpK[s]->coldim(), false, true);
        // Place things back in Fortran order.
        J.slice(s) = mJ.t();
        K.slice(s) = mK.t();
    }

    return;
}
