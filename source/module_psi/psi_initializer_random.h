#ifndef PSI_INITIALIZER_RANDOM_H
#define PSI_INITIALIZER_RANDOM_H

#include "psi_initializer.h"

class psi_initializer_random : public psi_initializer
{
    public:
        psi_initializer_random(Structure_Factor* sf_in, ModulePW::PW_Basis_K* pw_wfc_in);
        ~psi_initializer_random();

        // methods
        void random(std::complex<double>* psi,
                    const int iw_start,
                    const int iw_end,
                    const int ik,
                    const ModulePW::PW_Basis_K* wfc_basis);
        psi::Psi<std::complex<double>>* cal_psig(int ik) override;
};
#endif