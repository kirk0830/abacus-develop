#ifndef PSI_INITIALIZER_ATOMIC_H
#define PSI_INITIALIZER_ATOMIC_H
#include "psi_initializer.h"
#include "module_base/realarray.h"

/*
Psi (planewave based wavefunction) initializer: atomic
*/
template <typename T, typename Device>
class PsiInitializerAtomic : public PsiInitializer<T, Device>
{
    private:
        using Real = typename GetTypeReal<T>::type;
    public:
        PsiInitializerAtomic() {this->set_method("atomic");}
        ~PsiInitializerAtomic() {};

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
        virtual void allocate_table() override;
        virtual void tabulate() override;
        virtual void proj_ao_onkG(const int ik) override;
        // additional getter
        std::vector<std::string> pseudopot_files() const { return pseudopot_files_; }

    private:
        std::vector<std::string> pseudopot_files_;
        ModuleBase::realArray ovlp_pswfcjlq_;
};
#endif