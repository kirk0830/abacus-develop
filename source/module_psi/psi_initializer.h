// three global variables definition
#include "module_base/global_variable.h"
#include "module_base/global_function.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
// data structure support
#include "module_psi/psi.h" // for psi data structure
#include "module_hamilt_pw/hamilt_pwdft/structure_factor.h"
#include "module_basis/module_pw/pw_basis_k.h" // for kpoint related data structure
// basic math support
#include "module_base/complexmatrix.h"
#include "module_base/realarray.h"
#include "module_base/vector3.h"
// numerical algorithm support
#include "module_base/math_polyint.h" // for polynomial interpolation
#include "module_base/math_ylmreal.h" // for real spherical harmonics
#include "module_base/spherical_bessel_transformer.h" // for spherical bessel transform
#include "module_base/math_integral.h" // for numerical integration
/*
Psi (planewave based wavefunction) initializer
Auther: Kirk0830
Institute: AI for Science Institute, BEIJING

This class is used to allocate memory and give initial guess for Psi.
Following methods are available:
    1. random: use random number to initialize psi
               implemented in psi_initializer_random.h
    2. atomic: use pseudo-wavefunction in pseudopotential file to initialize psi
               implemented in psi_initializer_atomic.h
    3. atomic+random: mix 'atomic' with some random numbers to initialize psi
                      not implemented yet
    4. nao: use numerical orbitals to initialize psi
            implemented in psi_initializer_nao.h
    5. nao+random: mix 'nao' with some random numbers to initialize psi
                   not implemented yet
*/
class psi_initializer
{
    public:
        psi_initializer(Structure_Factor* sf_in, ModulePW::PW_Basis_K* pw_wfc_in) : sf(sf_in), pw_wfc(pw_wfc_in) { };
        ~psi_initializer() { };

        /// @brief calculate number of wavefunctions to be initialized
        /// @return if GlobalV::init_wfc has valid value, return the number of wavefunctions to be initialized, otherwise throw an error
        int get_starting_nw() const;

        /// @brief allocate memory for psi
        /// @return pointer to psi, memory allocated
        psi::Psi<std::complex<double>>* allocate();

        /// @brief initialize planewave represented psi
        /// @param psi psi
        /// @param ik index of kpoint
        virtual void initialize(psi::Psi<std::complex<double>>& psi, int ik) = 0;

        /// @brief initialize planewave represented psi after change of cell volume
        /// @note due to wanf2, the Wannier function in pw representation, is not used anymore, this function is not needed anymore
        void init_after_vc(); // in the past this function will refresh the wanf2, the Wannier function in pw representation. However for now it is not needed, because lcao_in_pw now is designed for initial guess

        Structure_Factor* sf;
        ModulePW::PW_Basis_K* pw_wfc; // I dont think it should appear here. It should, only be modified by ESolver object.

        ModuleBase::SphericalBesselTransformer sbt; // useful for atomic-like methods
    private:
        int mem_saver = 0; // will deprecated this variable soon
};