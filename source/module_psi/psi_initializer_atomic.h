#ifndef PSI_INITIALIZER_ATOMIC_H
#define PSI_INITIALIZER_ATOMIC_H
#include "psi_initializer.h"
#include "module_base/complexmatrix.h"
#include "module_base/realarray.h"

class psi_initializer_atomic : public psi_initializer
{
    public:
        psi_initializer_atomic(Structure_Factor* sf_in, ModulePW::PW_Basis_K* pw_wfc_in);
        ~psi_initializer_atomic();

        // methods
        void initialize(psi::Psi<std::complex<double>>& psi, int ik) override;

        // setters
        void set_pseudopot_files(std::string* pseudopot_files);
        // I wont write a function to set ovlp_pswfcjlq, it is totally useless
        void normalize_pswfc(int n_rgrid, double* pswfc, double* rgrid);
        /// @brief simple unitary phase factor
        /// @param arg the argument of the phase factor
        /// @param mode +1 for real part, -1 for imaginary part, 0 for the whole
        /// @return the phase factor
        std::complex<double> phase_factor(double arg, int mode = 0);
        /// @brief calculate the overlap between pseudo atomic wavefunctions and planewave basis
        void cal_ovlp_pswfcjlq();

        // historically left functions
        // getters
        std::vector<std::string> get_pseudopot_files() const { return pseudopot_files; }
        ModuleBase::realArray get_ovlp_pswfcjlq() const { return ovlp_pswfcjlq; }
    private:
        std::vector<std::string> pseudopot_files;
        ModuleBase::realArray ovlp_pswfcjlq;
};
#endif