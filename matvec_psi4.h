#ifndef MATVEC_PSI4_H_
#define MATVEC_PSI4_H_

#include <libresponse/matvec_i.h>

#include "psi4/libfock/jk.h"
#include "psi4/libmints/wavefunction.h"
#include "psi4/liboptions/liboptions.h"

using namespace psi;

class MatVec_Psi4 : public MatVec_i {

public:
    MatVec_Psi4(SharedWavefunction wfn, Options &options);
    ~MatVec_Psi4();

    void compute(arma::cube &J,
                 arma::cube &K,
                 const std::vector<arma::mat> &L,
                 const std::vector<arma::mat> &R);

private:
    SharedWavefunction wfn_;
    Options &options_;
    std::shared_ptr<JK> jk_;
};

#endif // MATVEC_PSI4_H_
