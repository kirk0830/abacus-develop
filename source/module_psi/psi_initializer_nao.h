#ifndef PSI_INITIALIZER_NAO_H
#define PSI_INITIALIZER_NAO_H
#include "psi_initializer.h"

template<typename FPTYPE>
class psi_initializer_nao : public psi_initializer<FPTYPE>
{
    public:
        psi_initializer_nao(Structure_Factor* sf_in, ModulePW::PW_Basis_K* pw_wfc_in);
        ~psi_initializer_nao();

        // methods
        psi::Psi<std::complex<FPTYPE>>* cal_psig(int ik) override;

        // setters
        void set_orbital_files(std::string* orbital_files);
        // I wont write a function to set ovlp_flzjlq, it is totally useless
        
        /// @brief calculate overlap integral between f_{l\\zeta} the radial numerical orbital and spherical Bessel function
        void cal_ovlp_flzjlq();

        // historically left functions
        void integral(const int meshr, const double *psir, const double *r,
            const double *rab, const int &l, double* table);
        double cosine_interpolation(const double &energy_x, const double &ecut,
            const double &beta);
        
        // getters
        std::vector<std::string> get_orbital_files() const { return orbital_files; }
        ModuleBase::realArray get_ovlp_flzjlq() const { return ovlp_flzjlq; }
    private:
        std::vector<std::string> orbital_files;
        ModuleBase::realArray ovlp_flzjlq;
};
#endif