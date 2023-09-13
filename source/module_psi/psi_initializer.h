#ifndef PSI_INITIALIZER_H
#define PSI_INITIALIZER_H
// three global variables definition
#include "module_base/global_variable.h"
#include "module_base/global_function.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
// timer and memory support
#include "module_base/timer.h"
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
#ifdef __MPI
// MPI support
#include <mpi.h>
#include "module_base/parallel_global.h"
#include "module_cell/parallel_kpoints.h"
#endif
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
        psi_initializer() : sf(nullptr), pw_wfc(nullptr) { };
        psi_initializer(Structure_Factor* sf_in, ModulePW::PW_Basis_K* pw_wfc_in);
        ~psi_initializer() { delete[] this->ixy2is; };

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
        
        /// @brief get method of initializing psi
        /// @return the method
        std::string get_method() const { return this->method; }
        /// @brief get complement number of bands
        /// @return nbands_complem
        int get_nbands_complem() const { return this->nbands_complem; }

        void print_status(psi::Psi<std::complex<double>>& psi) const;

        /// @brief set method manually
        /// @param method_in initialization method
        void set_method(std::string method_in) { this->method = method_in; }
        /// @brief set number of complementary bands
        /// @param nbands_in nbands_complem
        void set_nbands_complem(int nbands_in) { this->nbands_complem = nbands_in; }

        // virtual functions, will be implemented in derived classes
        // random
        virtual void random(std::complex<double>* psi,
                            const int iw_start,
                            const int iw_end,
                            const int ik,
                            const ModulePW::PW_Basis_K* wfc_basis) { ModuleBase::WARNING_QUIT("psi_initializer::random", "Polymorphism error"); }
        virtual void random(std::complex<float>* psi,
                            const int iw_start,
                            const int iw_end,
                            const int ik,
                            const ModulePW::PW_Basis_K* wfc_basis) { ModuleBase::WARNING_QUIT("psi_initializer::random", "Polymorphism error"); }
        #ifdef __MPI
        void stick_to_pool(float* stick, const int& ir, float* out, const ModulePW::PW_Basis_K* wfc_basis) const;
        void stick_to_pool(double* stick, const int& ir, double* out, const ModulePW::PW_Basis_K* wfc_basis) const;
        #endif
        template <typename FPTYPE>
        void random_t(std::complex<FPTYPE>* psi, const int iw_start, const int iw_end, const int ik, const ModulePW::PW_Basis_K* wfc_basis)
        {
            ModuleBase::timer::tick("psi_initializer", "random_t");
            assert(iw_start >= 0);
            const int ng = wfc_basis->npwk[ik];
        #ifdef __MPI
            if (INPUT.pw_seed > 0) // qianrui add 2021-8-13
            {
                srand(unsigned(INPUT.pw_seed + GlobalC::Pkpoints.startk_pool[GlobalV::MY_POOL] + ik));
                const int nxy = wfc_basis->fftnxy;
                const int nz = wfc_basis->nz;
                const int nstnz = wfc_basis->nst*nz;

                FPTYPE *stickrr = new FPTYPE[nz];
                FPTYPE *stickarg = new FPTYPE[nz];
                FPTYPE *tmprr = new FPTYPE[nstnz];
                FPTYPE *tmparg = new FPTYPE[nstnz];
                for (int iw = iw_start ;iw < iw_end;iw++)
                {   
                    // get the starting memory address of iw band
                    std::complex<FPTYPE>* psi_slice = &(psi[iw * this->pw_wfc->npwk_max * GlobalV::NPOL]);
                    int startig = 0;
                    for(int ipol = 0 ; ipol < GlobalV::NPOL ; ++ipol)
                    {
                    
                        for(int ir=0; ir < nxy; ir++)
                        {
                            if(wfc_basis->fftixy2ip[ir] < 0) continue;
                            if(GlobalV::RANK_IN_POOL==0)
                            {
                                for(int iz=0; iz<nz; iz++)
                                {
                                    stickrr[ iz ] = std::rand()/FPTYPE(RAND_MAX);
                                    stickarg[ iz ] = std::rand()/FPTYPE(RAND_MAX);
                                }
                            }
                            stick_to_pool(stickrr, ir, tmprr, wfc_basis);
                            stick_to_pool(stickarg, ir, tmparg, wfc_basis);
                        }

                        for (int ig = 0;ig < ng;ig++)
                        {
                            const FPTYPE rr = tmprr[wfc_basis->getigl2isz(ik,ig)];
                            const FPTYPE arg= ModuleBase::TWO_PI * tmparg[wfc_basis->getigl2isz(ik,ig)];
                            const FPTYPE gk2 = wfc_basis->getgk2(ik,ig);
                            psi_slice[ig+startig] = std::complex<FPTYPE>(rr * cos(arg), rr * sin(arg)) / FPTYPE(gk2 + 1.0);
                        }
                        startig += this->pw_wfc->npwk_max;
                    }
                }
                delete[] stickrr;
                delete[] stickarg;
                delete[] tmprr;
                delete[] tmparg;
            }
            else
            {
        #else  // !__MPI
            if (INPUT.pw_seed > 0) // qianrui add 2021-8-13
            {
                srand(unsigned(INPUT.pw_seed + ik));
            }
        #endif
                for (int iw = iw_start ;iw < iw_end;iw++)
                {
                    std::complex<FPTYPE>* psi_slice = &(psi[iw * this->pw_wfc->npwk_max * GlobalV::NPOL]);
                    for (int ig = 0; ig < ng; ig++)
                    {
                        const FPTYPE rr = std::rand()/FPTYPE(RAND_MAX); //qianrui add RAND_MAX
                        const FPTYPE arg= ModuleBase::TWO_PI * std::rand()/FPTYPE(RAND_MAX);
                        const FPTYPE gk2 = wfc_basis->getgk2(ik,ig);
                        psi_slice[ig] = std::complex<FPTYPE>(rr * cos(arg), rr * sin(arg)) / FPTYPE(gk2 + 1.0);
                    }
                    if(GlobalV::NPOL==2)
                    {
                        for (int ig = this->pw_wfc->npwk_max; ig < this->pw_wfc->npwk_max + ng; ig++)
                        {
                            const FPTYPE rr = std::rand()/FPTYPE(RAND_MAX);
                            const FPTYPE arg= ModuleBase::TWO_PI * std::rand()/FPTYPE(RAND_MAX);
                            const FPTYPE gk2 = wfc_basis->getgk2(ik,ig-this->pw_wfc->npwk_max);
                            psi_slice[ig] = std::complex<FPTYPE>(rr * cos(arg), rr * sin(arg)) / FPTYPE(gk2 + 1.0);
                        }
                    }

                }
        #ifdef __MPI
            }
        #endif
            ModuleBase::timer::tick("psi_initializer_random", "random_t");
        }
        
        // atomic
        virtual void set_pseudopot_files(std::string* pseudopot_files) { ModuleBase::WARNING_QUIT("psi_initializer::set_pseudopot_files", "Polymorphism error"); }
        virtual void normalize_pswfc(int n_rgrid, double* pswfc, double* rgrid) { ModuleBase::WARNING_QUIT("psi_initializer::normalize_pswfc", "Polymorphism error"); }
        virtual std::complex<double> phase_factor(double arg, int mode = 0) { ModuleBase::WARNING_QUIT("psi_initializer::phase_factor", "Polymorphism error"); return std::complex<double>(0.0,0.0);}
        virtual void cal_ovlp_pswfcjlq() { ModuleBase::WARNING_QUIT("psi_initializer::calc_ovlp_pswfcjlq", "Polymorphism error"); }
        // nao
        virtual void set_orbital_files(std::string* orbital_files) { ModuleBase::WARNING_QUIT("psi_initializer::set_orbital_files", "Polymorphism error"); }
        virtual void cal_ovlp_flzjlq() { ModuleBase::WARNING_QUIT("psi_initializer::cal_ovlp_flzjlq", "Polymorphism error"); }
        // atomic+random
        // nao+random
        Structure_Factor* sf;
        ModulePW::PW_Basis_K* pw_wfc; // I dont think it should appear here. It should, only be modified by ESolver object.

        ModuleBase::SphericalBesselTransformer sbt; // useful for atomic-like methods
    private:
        int mem_saver = 0; // will deprecated this variable soon
        std::string method = "none";

        int nbands_complem = 0;
        int* ixy2is;
};
#endif