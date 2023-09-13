#ifndef PSI_INITIALIZER_NAO_RANDOM_H
#define PSI_INITIALIZER_NAO_RANDOM_H
#include "psi_initializer_nao.h"

template <typename FPTYPE>
class psi_initializer_nao_random : public psi_initializer_nao<FPTYPE>
{
    public:
        psi_initializer_nao_random(Structure_Factor* sf_in, ModulePW::PW_Basis_K* pw_wfc_in);
        ~psi_initializer_nao_random();
        psi::Psi<std::complex<FPTYPE>>* cal_psig(int ik) override;
    private:
};
#endif