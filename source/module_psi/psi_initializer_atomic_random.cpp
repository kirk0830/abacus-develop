#include "psi_initializer_atomic_random.h"

#ifdef __MPI
template <typename T, typename Device>
void PsiInitializerAtomicRandom<T, Device>::initialize(Structure_Factor* sf,                                //< structure factor
                                                          ModulePW::PW_Basis_K* pw_wfc,                        //< planewave basis
                                                          UnitCell* p_ucell,                                   //< unit cell
                                                          Parallel_Kpoints* p_parakpts,                        //< parallel kpoints
                                                          const int& random_seed,                          //< random seed
                                                          pseudopot_cell_vnl* p_pspot_nl,
                                                          const int& rank)
{
    PsiInitializerAtomic<T, Device>::initialize(sf, pw_wfc, p_ucell, p_parakpts, random_seed, p_pspot_nl, rank);
}
#else
template <typename T, typename Device>
void PsiInitializerAtomicRandom<T, Device>::initialize(Structure_Factor* sf,                                //< structure factor
                                                          ModulePW::PW_Basis_K* pw_wfc,                        //< planewave basis
                                                          UnitCell* p_ucell,                                   //< unit cell
                                                          const int& random_seed,                          //< random seed
                                                          pseudopot_cell_vnl* p_pspot_nl)
{
    PsiInitializerAtomic<T, Device>::initialize(sf, pw_wfc, p_ucell, random_seed, p_pspot_nl);
}
#endif

template <typename T, typename Device>
void PsiInitializerAtomicRandom<T, Device>::proj_ao_onkG(const int ik)
{
    double rm = this->random_mix();
    const int ik_psig = (this->d_psig_->get_nk() == 1) ? 0 : ik;
    this->d_psig_->fix_k(ik_psig);
    PsiInitializerAtomic<T, Device>::proj_ao_onkG(ik);
    psi::Psi<T, Device> psi_random(1, this->d_psig_->get_nbands(), this->d_psig_->get_nbasis(), nullptr);
    psi_random.fix_k(0);
    this->random_t(psi_random.get_pointer(), 0, psi_random.get_nbands(), ik);
    
    #pragma omp parallel for collapse(2)
    for(int iband = 0; iband < this->d_psig_->get_nbands(); iband++)
    {
        for(int ibasis = 0; ibasis < this->d_psig_->get_nbasis(); ibasis++)
        {
            (*(this->d_psig_))(iband, ibasis) = ((Real)(1-rm))*(*(this->d_psig_))(iband, ibasis) + ((Real)rm)*psi_random(iband, ibasis);
        }
    }
}

template class PsiInitializerAtomicRandom<std::complex<double>, base_device::DEVICE_CPU>;
template class PsiInitializerAtomicRandom<std::complex<float>, base_device::DEVICE_CPU>;
// gamma point calculation
template class PsiInitializerAtomicRandom<double, base_device::DEVICE_CPU>;
template class PsiInitializerAtomicRandom<float, base_device::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
template class PsiInitializerAtomicRandom<std::complex<double>, base_device::DEVICE_GPU>;
template class PsiInitializerAtomicRandom<std::complex<float>, base_device::DEVICE_GPU>;
// gamma point calculation
template class PsiInitializerAtomicRandom<double, base_device::DEVICE_GPU>;
template class PsiInitializerAtomicRandom<float, base_device::DEVICE_GPU>;
#endif
