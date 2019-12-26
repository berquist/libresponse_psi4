#ifndef PARSE_OPERATORS_H_
#define PARSE_OPERATORS_H_

#include "psi4/libmints/wavefunction.h"
#include "psi4/liboptions/liboptions.h"
#include "libresponse/operator_spec.h"

void parse_operators(
    psi::SharedWavefunction &ref_wfn,
    psi::Options &options,
    std::vector<libresponse::operator_spec> &operators);

#endif // PARSE_OPERATORS_H_
