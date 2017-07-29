#ifndef WRAPPERS_H_
#define WRAPPERS_H_

#include <armadillo>

#include "psi4/libmints/wavefunction.h"

using namespace psi;


void AO_dipole_length(arma::cube &M_AO, const arma::vec &origin, SharedWavefunction wfn, Options &options);
void AO_quadrupole_length(arma::cube &M_AO, const arma::vec &origin, SharedWavefunction wfn, Options &options);
void AO_multipole(arma::cube &M_AO, const arma::uvec &order, const arma::vec &origin, SharedWavefunction wfn, Options &options);
void AO_multipole(arma::cube &M_AO, const arma::uvec &order, const arma::vec &origin, SharedWavefunction wfn, Options &options);
void AO_dipole_velocity(arma::cube &D_AO, SharedWavefunction wfn, Options &options);
void AO_angular_momentum(arma::cube &L_AO, const arma::vec &origin, SharedWavefunction wfn, Options &options);

#endif // WRAPPERS_H_
