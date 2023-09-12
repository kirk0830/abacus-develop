#include "psi_initializer.h"

class psi_initializer_random : public psi_initializer
{
    public:
        psi_initializer_random(Structure_Factor* sf_in, ModulePW::PW_Basis_K* pw_wfc_in);
        ~psi_initializer_random();

        // methods
        void random(std::complex<double>* psi,
                    const int iw_start,
                    const int iw_end,
                    const int ik,
                    const ModulePW::PW_Basis_K* wfc_basis);
        void random(std::complex<float>* psi,
                    const int iw_start,
                    const int iw_end,
                    const int ik,
                    const ModulePW::PW_Basis_K* wfc_basis);
        template <typename FPTYPE>
        void random_t(std::complex<FPTYPE>* psi,
                      const int iw_start,
                      const int iw_end,
                      const int ik,
                      const ModulePW::PW_Basis_K* wfc_basis);
        void initialize(psi::Psi<std::complex<double>>& psi, int ik) override;
#ifdef __MPI
        void stick_to_pool(float* stick, const int& ir, float* out, const ModulePW::PW_Basis_K* wfc_basis) const;
        void stick_to_pool(double* stick, const int& ir, double* out, const ModulePW::PW_Basis_K* wfc_basis) const;
#endif
    private:
        int* ixy2is;
};