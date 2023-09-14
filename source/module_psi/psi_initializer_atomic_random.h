#ifndef PSI_INITIALIZER_ATOMIC_RANDOM_H
#define PSI_INITIALIZER_ATOMIC_RANDOM_H
#include "psi_initializer_atomic.h"

class psi_initializer_atomic_random : public psi_initializer_atomic
{
    public:
        psi_initializer_atomic_random(Structure_Factor* sf_in, ModulePW::PW_Basis_K* pw_wfc_in);
        ~psi_initializer_atomic_random();
        psi::Psi<std::complex<double>>* cal_psig(int ik) override;
    private:

};
#endif