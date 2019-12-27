#include <armadillo>

#include "psi4/libmints/basisset.h"
#include "psi4/libmints/integral.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/mintshelper.h"
#include "psi4/libmints/vector3.h"

using namespace psi;

void AO_dipole_length(arma::cube &M_AO,
                      const arma::vec &origin,
                      SharedWavefunction wfn,
                      Options &options) {

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

    M_AO.set_size(nbf, nbf, v_dipole.size());

    for (size_t c = 0; c < v_dipole.size(); c++) {
        SharedMatrix msm_dipole = v_dipole[c];
        arma::mat ma_dipole(msm_dipole->get_pointer(),
                            msm_dipole->rowdim(),
                            msm_dipole->coldim(),
                            false,
                            true);
        M_AO.slice(c) = ma_dipole;
    }

    return;
}

void AO_quadrupole_length(arma::cube &M_AO,
                          const arma::vec &origin,
                          SharedWavefunction wfn,
                          Options &options) {

    const Vector3 sm_origin(origin.memptr());

    const size_t nbf = wfn->basisset()->nbf();

    std::shared_ptr<IntegralFactory> intfactory = wfn->integral();

    std::vector<SharedMatrix> v_quadrupole;

    v_quadrupole.push_back(SharedMatrix(new Matrix("AO Quadrupole XX", nbf, nbf)));
    v_quadrupole.push_back(SharedMatrix(new Matrix("AO Quadrupole XY", nbf, nbf)));
    v_quadrupole.push_back(SharedMatrix(new Matrix("AO Quadrupole XZ", nbf, nbf)));
    v_quadrupole.push_back(SharedMatrix(new Matrix("AO Quadrupole YY", nbf, nbf)));
    v_quadrupole.push_back(SharedMatrix(new Matrix("AO Quadrupole YZ", nbf, nbf)));
    v_quadrupole.push_back(SharedMatrix(new Matrix("AO Quadrupole ZZ", nbf, nbf)));

    std::shared_ptr<OneBodyAOInt> ints(intfactory->ao_quadrupole());
    ints->set_origin(sm_origin);
    ints->compute(v_quadrupole);

    M_AO.set_size(nbf, nbf, v_quadrupole.size());

    for (size_t c = 0; c < v_quadrupole.size(); c++) {
        SharedMatrix msm_quadrupole = v_quadrupole[c];
        arma::mat ma_quadrupole(msm_quadrupole->get_pointer(),
                                msm_quadrupole->rowdim(),
                                msm_quadrupole->coldim(),
                                false,
                                true);
        M_AO.slice(c) = ma_quadrupole;
    }

    return;
}

void AO_multipole(arma::cube &M_AO,
                  const arma::uvec &order,
                  const arma::vec &origin,
                  SharedWavefunction wfn,
                  Options &options) {

    // TODO
    return;
}

void AO_dipole_velocity(arma::cube &D_AO, SharedWavefunction wfn, Options &options) {

    const size_t nbf = wfn->basisset()->nbf();

    std::shared_ptr<IntegralFactory> intfactory = wfn->integral();

    std::vector<SharedMatrix> v_nabla;

    v_nabla.push_back(SharedMatrix(new Matrix("AO Px", nbf, nbf)));
    v_nabla.push_back(SharedMatrix(new Matrix("AO Py", nbf, nbf)));
    v_nabla.push_back(SharedMatrix(new Matrix("AO Pz", nbf, nbf)));

    std::shared_ptr<OneBodyAOInt> ints(intfactory->ao_nabla());
    ints->compute(v_nabla);

    D_AO.set_size(nbf, nbf, v_nabla.size());

    for (size_t c = 0; c < v_nabla.size(); c++) {
        SharedMatrix msm_nabla = v_nabla[c];
        arma::mat ma_nabla(msm_nabla->get_pointer(),
                           msm_nabla->rowdim(),
                           msm_nabla->coldim(),
                           false,
                           true);
        D_AO.slice(c) = ma_nabla;
        // TODO check for antisymmetry due to imaginary operator
    }

    return;
}

void AO_angular_momentum(arma::cube &L_AO,
                         const arma::vec &origin,
                         SharedWavefunction wfn,
                         Options &options) {

    const Vector3 sm_origin(origin.memptr());

    const size_t nbf = wfn->basisset()->nbf();

    std::shared_ptr<IntegralFactory> intfactory = wfn->integral();

    std::vector<SharedMatrix> v_angmom;

    v_angmom.push_back(SharedMatrix(new Matrix("AO Lx", nbf, nbf)));
    v_angmom.push_back(SharedMatrix(new Matrix("AO Ly", nbf, nbf)));
    v_angmom.push_back(SharedMatrix(new Matrix("AO Lz", nbf, nbf)));

    std::shared_ptr<OneBodyAOInt> ints(intfactory->ao_angular_momentum());
    ints->set_origin(sm_origin);
    ints->compute(v_angmom);

    L_AO.set_size(nbf, nbf, v_angmom.size());

    for (size_t c = 0; c < v_angmom.size(); c++) {
        SharedMatrix msm_angmom = v_angmom[c];
        arma::mat ma_angmom(msm_angmom->get_pointer(),
                            msm_angmom->rowdim(),
                            msm_angmom->coldim(),
                            false,
                            true);
        L_AO.slice(c) = ma_angmom;
        // TODO check for antisymmetry due to imaginary operator
    }

    return;
}
