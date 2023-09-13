#ifndef PSI_INITIALIZER_ATOMIC_RANDOM_H
#define PSI_INITIALIZER_ATOMIC_RANDOM_H
#include "psi_initializer_atomic.h"

template <typename FPTYPE>
class psi_initializer_atomic_random : public psi_initializer_atomic<FPTYPE>
{
    public:
        psi_initializer_atomic_random(Structure_Factor* sf_in, ModulePW::PW_Basis_K* pw_wfc_in);
        ~psi_initializer_atomic_random();
        void initialize(psi::Psi<std::complex<double>>& psi, int ik);
    private:

};
#endif