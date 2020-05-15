#include <libresponse/operator_spec.h>

#include "psi4/libmints/wavefunction.h"
#include "psi4/liboptions/liboptions.h"

#include "wrappers.h"

using namespace libresponse;

void parse_operators(psi::SharedWavefunction &ref_wfn,
                     psi::Options &options,
                     std::vector<operator_spec> &operators) {

    arma::cube integrals;
    arma::vec origin(3, arma::fill::zeros);
    std::string origin_label;

    const bool do_response = true;

    if (options["OPERATOR_DIPOLE"].has_changed()) {
        std::cout << " Adding operator: dipole (length gauge)" << std::endl;
        const std::string operator_label = "dipole";
        AO_dipole_length(integrals, origin, ref_wfn, options);
        operator_metadata metadata(operator_label, origin_label, -1, false, false);
        operator_spec operator_(metadata, integrals, origin, do_response);
        operators.push_back(operator_);
    }

    if (options["OPERATOR_QUADRUPOLE"].has_changed()) {
        std::cout << " Adding operator: quadrupole" << std::endl;
        const std::string operator_label = "quadrupole";
        AO_quadrupole_length(integrals, origin, ref_wfn, options);
        operator_metadata metadata(operator_label, origin_label, -1, false, false);
        operator_spec operator_(metadata, integrals, origin, do_response);
        operators.push_back(operator_);
    }

    // if (options["OPERATOR_MULTIPOLE"].has_changed()) {
    // }

    if (options["OPERATOR_NABLA"].has_changed()) {
        std::cout << " Adding operator: dipole (velocity gauge) / linear momentum" << std::endl;
        const std::string operator_label = "dipvel";
        AO_dipole_velocity(integrals, ref_wfn, options);
        operator_metadata metadata(operator_label, origin_label, -1, true, false);
        operator_spec operator_(metadata, integrals, origin, do_response);
        operators.push_back(operator_);
    }

    if (options["OPERATOR_ANGMOM"].has_changed()) {
        std::cout << " Adding operator: angular momentum" << std::endl;
        const std::string operator_label = "angmom";
        AO_angular_momentum(integrals, origin, ref_wfn, options);
        operator_metadata metadata(operator_label, origin_label, -1, true, false);
        operator_spec operator_(metadata, integrals, origin, do_response);
        operators.push_back(operator_);
    }

    return;
}
