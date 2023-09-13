#include "psi_initializer_random.h"
#include <mpi.h>
#include "module_base/parallel_global.h"
#include "module_cell/parallel_kpoints.h"

template <typename FPTYPE>
psi_initializer_random<FPTYPE>::psi_initializer_random(Structure_Factor* sf_in, ModulePW::PW_Basis_K* pw_wfc_in) : psi_initializer<FPTYPE>(sf_in, pw_wfc_in)
{
    this->set_method("random");
}

template <typename FPTYPE>
psi_initializer_random<FPTYPE>::~psi_initializer_random()
{
}

template <typename FPTYPE>
void psi_initializer_random<FPTYPE>::random(std::complex<double>* psi,
                                    const int iw_start,
                                    const int iw_end,
                                    const int ik,
                                    const ModulePW::PW_Basis_K* wfc_basis)
{
    ModuleBase::timer::tick("psi_initializer_random", "random");
    this->random_t(psi, iw_start, iw_end, ik, wfc_basis);
    ModuleBase::timer::tick("psi_initializer_random", "random");
}

template <typename FPTYPE>
void psi_initializer_random<FPTYPE>::random(std::complex<float>* psi,
                                    const int iw_start,
                                    const int iw_end,
                                    const int ik,
                                    const ModulePW::PW_Basis_K* wfc_basis)
{
    ModuleBase::timer::tick("psi_initializer_random", "random");
    this->random_t(psi, iw_start, iw_end, ik, wfc_basis);
    ModuleBase::timer::tick("psi_initializer_random", "random");
}

template <typename FPTYPE>
psi::Psi<std::complex<FPTYPE>>* psi_initializer_random<FPTYPE>::cal_psig(int ik)
{
    ModuleBase::timer::tick("psi_initializer_random", "initialize");
    //this->print_status(psi);
    this->psig->fix_k(ik);
    this->random(this->psig->get_pointer(), 0, this->psig->get_nbands(), ik, this->pw_wfc);
    // we still need to diagonalize the obtained psi from hsolver::DiagoIterAssist::diagH_subspace
    // will do it in HSolver function...
    ModuleBase::timer::tick("psi_initializer_random", "initialize");
    return this->psig;
}