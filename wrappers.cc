#include <armadillo>

#include "psi4/libmints/basisset.h"
#include "psi4/libmints/integral.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/mintshelper.h"
#include "psi4/libmints/vector3.h"

using namespace psi;


void AO_dipole_length(arma::cube &M_AO, const arma::vec &origin, SharedWavefunction wfn, Options &options) {

    const Vector3 sm_origin(origin.memptr());

    const size_t nbf = wfn->basisset()->nbf();

    std::shared_ptr<IntegralFactory> intfactory = wfn->integral();

    std::vector<SharedMatrix> v_dipole;

    v_dipole.push_back(SharedMatrix(new Matrix("AO Mux", nbf, nbf)));
    v_dipole.push_back(SharedMatrix(new Matrix("AO Muy", nbf, nbf)));
    v_dipole.push_back(SharedMatrix(new Matrix("AO Muz", nbf, nbf)));

    std::shared_ptr<OneBodyAOInt> ints(intfactory->ao_dipole());
    ints->set_origin(sm_origin);
    ints->compute(v_dipole);

    M_AO.set_size(nbf, nbf, 3);

    for (size_t c = 0; c < 3; c++) {
        SharedMatrix msm_dipole = v_dipole[c];
        arma::mat ma_dipole(msm_dipole->get_pointer(), msm_dipole->rowdim(), msm_dipole->coldim(), false, true);
        M_AO.slice(c) = ma_dipole;
    }

    return;

}

void AO_quadrupole_length(arma::cube &M_AO, const arma::vec &origin, SharedWavefunction wfn, Options &options) {

    return;

}

void AO_multipole(arma::cube &M_AO, const arma::uvec &order, const arma::vec &origin, SharedWavefunction wfn, Options &options) {

    return;

}

void AO_dipole_velocity(arma::cube &D_AO, SharedWavefunction wfn, Options &options) {

    return;

}

void AO_angular_momentum(arma::cube &L_AO, const arma::vec &origin, SharedWavefunction wfn, Options &options) {

    return;

}
