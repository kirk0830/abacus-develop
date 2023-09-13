#ifndef PSI_INITIALIZER_NAO_RANDOM_H
#define PSI_INITIALIZER_NAO_RANDOM_H
#include "psi_initializer_nao.h"

class psi_initializer_nao_random : public psi_initializer_nao
{
    public:
        psi_initializer_nao_random(Structure_Factor* sf_in, ModulePW::PW_Basis_K* pw_wfc_in);
        ~psi_initializer_nao_random();
        void initialize(psi::Psi<std::complex<double>>& psi, int ik);
    private:
};
#endif