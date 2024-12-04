#ifndef PSI_INITIALIZER_ATOMIC_RANDOM_H
#define PSI_INITIALIZER_ATOMIC_RANDOM_H
#include "module_psi/psi_initializer_atomic.h"
#include "module_hamilt_pw/hamilt_pwdft/VNL_in_pw.h"
#include "module_cell/parallel_kpoints.h"

/*
Psi (planewave based wavefunction) initializer: atomic+random
*/
template <typename T, typename Device>
class PsiInitializerAtomicRandom : public PsiInitializerAtomic<T, Device>
{
    private:
        using Real = typename GetTypeReal<T>::type;
    public:
        PsiInitializerAtomicRandom() {this->set_method("atomic+random"); this->set_random_mix(0.05);}
        ~PsiInitializerAtomicRandom() {};

        #ifdef __MPI // MPI additional implementation
        /// @brief initialize the PsiInitializer with external data and methods
        virtual void initialize(Structure_Factor*,                      //< structure factor
                                ModulePW::PW_Basis_K*,                  //< planewave basis
                                UnitCell*,                              //< unit cell
                                Parallel_Kpoints*,                      //< parallel kpoints
                                const int& = 1,                         //< random seed
                                pseudopot_cell_vnl* = nullptr,          //< nonlocal pseudopotential
                                const int& = 0) override;               //< MPI rank
        #else
        /// @brief serial version of initialize function, link PsiInitializer with external data and methods
        virtual void initialize(Structure_Factor*,                      //< structure factor
                                ModulePW::PW_Basis_K*,                  //< planewave basis
                                UnitCell*,                              //< unit cell
                                const int& = 1,                         //< random seed
                                pseudopot_cell_vnl* = nullptr) override;//< nonlocal pseudopotential
        #endif

        virtual void proj_ao_onkG(const int ik) override;
        virtual void tabulate() override {PsiInitializerAtomic<T, Device>::tabulate();};
        
    private:
};
#endif