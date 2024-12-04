#include "psi_initializer_nao_random.h"
#include "module_base/timer.h"

#ifdef __MPI
template <typename T, typename Device>
void PsiInitializerNAORandom<T, Device>::initialize(Structure_Factor* sf,
                                                       ModulePW::PW_Basis_K* pw_wfc,
                                                       UnitCell* p_ucell,
                                                       Parallel_Kpoints* p_parakpts,
                                                       const int& random_seed,
                                                       pseudopot_cell_vnl* p_pspot_nl,
                                                       const int& rank)
{
    PsiInitializerNAO<T, Device>::initialize(sf, pw_wfc, p_ucell, p_parakpts, random_seed, p_pspot_nl, rank);
}
#else
template <typename T, typename Device>
void PsiInitializerNAORandom<T, Device>::initialize(Structure_Factor* sf,
                                                       ModulePW::PW_Basis_K* pw_wfc,
                                                       UnitCell* p_ucell,
                                                       const int& random_seed,
                                                       pseudopot_cell_vnl* p_pspot_nl)
{
    PsiInitializerNAO<T, Device>::initialize(sf, pw_wfc, p_ucell, random_seed, p_pspot_nl);
}
#endif

template <typename T, typename Device>
void PsiInitializerNAORandom<T, Device>::proj_ao_onkG(const int ik)
{
    ModuleBase::timer::tick("PsiInitializerNAORandom", "proj_ao_onkG");
    double rm = this->random_mix();
    const int ik_psig = (this->d_psig_->get_nk() == 1) ? 0 : ik;
    this->d_psig_->fix_k(ik_psig);
    PsiInitializerNAO<T, Device>::proj_ao_onkG(ik);
    psi::Psi<T, Device> psi_random(1, this->d_psig_->get_nbands(), this->d_psig_->get_nbasis(), nullptr);
    psi_random.fix_k(0);
    this->random_t(psi_random.get_pointer(), 0, psi_random.get_nbands(), ik);
    for(int iband = 0; iband < this->d_psig_->get_nbands(); iband++)
    {
        for(int ibasis = 0; ibasis < this->d_psig_->get_nbasis(); ibasis++)
        {
            (*(this->d_psig_))(iband, ibasis) = ((Real)(1-rm))*(*(this->d_psig_))(iband, ibasis) + ((Real)rm)*psi_random(iband, ibasis);
        }
    }
    ModuleBase::timer::tick("PsiInitializerNAORandom", "proj_ao_onkG");
}

template class PsiInitializerNAORandom<std::complex<double>, base_device::DEVICE_CPU>;
template class PsiInitializerNAORandom<std::complex<float>, base_device::DEVICE_CPU>;
// gamma point calculation
template class PsiInitializerNAORandom<double, base_device::DEVICE_CPU>;
template class PsiInitializerNAORandom<float, base_device::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
template class PsiInitializerNAORandom<std::complex<double>, base_device::DEVICE_GPU>;
template class PsiInitializerNAORandom<std::complex<float>, base_device::DEVICE_GPU>;
// gamma point calculation
template class PsiInitializerNAORandom<double, base_device::DEVICE_GPU>;
template class PsiInitializerNAORandom<float, base_device::DEVICE_GPU>;
#endif