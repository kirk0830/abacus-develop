#ifndef PSI_INITIALIZER_H
#define PSI_INITIALIZER_H
// three global variables definition
#include "module_base/global_variable.h"
#include "module_base/global_function.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
// basic functions support
#include "module_base/tool_quit.h"
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

This class is used to allocate memory and give initial guess for psi (not kspw_psi the FPTYPE, Device template one)
therefore only double datatype is needed to be supported.
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
        /// @brief default constructor of psi initializer
        psi_initializer() : sf(nullptr), pw_wfc(nullptr) { };
        /// @brief parameterized constructor of psi initializer
        /// @param sf_in structure factor pointer
        /// @param pw_wfc_in pw_basis_k pointer
        psi_initializer(Structure_Factor* sf_in, ModulePW::PW_Basis_K* pw_wfc_in);
        /// @brief destructor
        ~psi_initializer();

        /// @brief calculate number of wavefunctions to be initialized
        /// @return if GlobalV::init_wfc has valid value, return the number of wavefunctions to be initialized, otherwise throw an error
        int get_starting_nw() const;

        /// @brief allocate memory for psi
        /// @return pointer to psi, memory allocated
        psi::Psi<std::complex<double>>* allocate();

        /// @brief calculate psi in planewave representation
        /// @param psi psi
        /// @param ik index of kpoint
        virtual psi::Psi<std::complex<double>>* cal_psig(int ik) = 0;

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
        /// @brief output psig to file, given number of kpoints, bands and basis, diagonalization method. Note: because it is complex number, therefore every number has format (real, imag)
        void write_psig() const;

        // virtual functions, will be implemented in derived classes

        // random to complement bands not initialized by pswfc or nao, therefore it is a basic function, or psi_initializer_random will be inherented by all other methods.
        /// @brief kernel to generate and assign random number for psi
        /// @param psi for psi::Psi<FPTYPE, Device>, first use psi.fix(ik), then psi.get_pointer() to pass to parameter list
        /// @param iw_start index of wavefunction (start), the index of first band to initialize
        /// @param iw_end index of wavefunction (end)
        /// @param ik index of kpoint
        /// @param wfc_basis because do not want to change contents in it, this parameter is reserved. should give this->pw_wfc
        void random_t(std::complex<double>* psi, const int iw_start, const int iw_end, const int ik, const ModulePW::PW_Basis_K* wfc_basis)
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

                double *stickrr = new double[nz];
                double *stickarg = new double[nz];
                double *tmprr = new double[nstnz];
                double *tmparg = new double[nstnz];
                for (int iw = iw_start; iw < iw_end; iw++)
                {   
                    // get the starting memory address of iw band
                    std::complex<double>* psi_slice = &(psi[iw * this->pw_wfc->npwk_max * GlobalV::NPOL]);
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
                                    stickrr[ iz ] = std::rand()/double(RAND_MAX);
                                    stickarg[ iz ] = std::rand()/double(RAND_MAX);
                                }
                            }
                            stick_to_pool(stickrr, ir, tmprr, wfc_basis);
                            stick_to_pool(stickarg, ir, tmparg, wfc_basis);
                        }

                        for (int ig = 0;ig < ng;ig++)
                        {
                            const double rr = tmprr[wfc_basis->getigl2isz(ik,ig)];
                            const double arg= ModuleBase::TWO_PI * tmparg[wfc_basis->getigl2isz(ik,ig)];
                            const double gk2 = wfc_basis->getgk2(ik,ig);
                            psi_slice[ig+startig] = std::complex<double>(rr * cos(arg), rr * sin(arg)) / double(gk2 + 1.0);
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
                for (int iw = iw_start ;iw < iw_end; iw++)
                {
                    std::complex<double>* psi_slice = &(psi[iw * this->pw_wfc->npwk_max * GlobalV::NPOL]);
                    for (int ig = 0; ig < ng; ig++)
                    {
                        const double rr = std::rand()/double(RAND_MAX); //qianrui add RAND_MAX
                        const double arg= ModuleBase::TWO_PI * std::rand()/double(RAND_MAX);
                        const double gk2 = wfc_basis->getgk2(ik,ig);
                        psi_slice[ig] = std::complex<double>(rr * cos(arg), rr * sin(arg)) / double(gk2 + 1.0);
                    }
                    if(GlobalV::NPOL==2)
                    {
                        for (int ig = this->pw_wfc->npwk_max; ig < this->pw_wfc->npwk_max + ng; ig++)
                        {
                            const double rr = std::rand()/double(RAND_MAX);
                            const double arg= ModuleBase::TWO_PI * std::rand()/double(RAND_MAX);
                            const double gk2 = wfc_basis->getgk2(ik,ig-this->pw_wfc->npwk_max);
                            psi_slice[ig] = std::complex<double>(rr * cos(arg), rr * sin(arg)) / double(gk2 + 1.0);
                        }
                    }

                }
        #ifdef __MPI
            }
        #endif
            ModuleBase::timer::tick("psi_initializer_random", "random_t");
        }
        // random
        /// @brief wrapper of random_t
        /// @param psi for psi::Psi<FPTYPE, Device>, first use psi.fix(ik), then psi.get_pointer() to pass to parameter list
        /// @param iw_start index of wavefunction (start), the index of first band to initialize
        /// @param iw_end index of wavefunction (end)
        /// @param ik index of kpoint
        /// @param wfc_basis because do not want to change contents in it, this parameter is reserved. should give this->pw_wfc
        virtual void random(std::complex<double>* psi, const int iw_start, const int iw_end,
                            const int ik, const ModulePW::PW_Basis_K* wfc_basis) { ModuleBase::WARNING_QUIT("psi_initializer::random", "Polymorphism error"); }
        #ifdef __MPI
        /// @brief (about planewaves distribution) from stick mapping to pool
        /// @param stick 
        /// @param ir 
        /// @param out 
        /// @param wfc_basis because do not want to change contents in it, this parameter is reserved. should give this->pw_wfc
        void stick_to_pool(double* stick, const int& ir, double* out, const ModulePW::PW_Basis_K* wfc_basis) const;
        #endif
        // atomic
        /// @brief setter of pseudopotential files, useful when init_wfc = atomic
        virtual void set_pseudopot_files(std::string* pseudopot_files) { ModuleBase::WARNING_QUIT("psi_initializer::set_pseudopot_files", "Polymorphism error"); }
        /// @brief normalize pseudo wavefunction
        /// @param n_rgrid level of realspace grid
        /// @param pswfc pseudowavefunction read from pseudopotential file
        /// @param rgrid realspace grid read from pseudopotential file
        virtual void normalize_pswfc(int n_rgrid, double* pswfc, double* rgrid) { ModuleBase::WARNING_QUIT("psi_initializer::normalize_pswfc", "Polymorphism error"); }
        /// @brief calculate cos(arg)+isin(arg)
        /// @param arg argument
        /// @param mode if 1, return cos(arg), 0, return cos(arg)+isin(arg), -1, return sin(arg)
        /// @return it depends
        virtual std::complex<double> phase_factor(double arg, int mode = 0) { ModuleBase::WARNING_QUIT("psi_initializer::phase_factor", "Polymorphism error"); return std::complex<double>(0.0,0.0);}
        /// @brief calculate overlap table between pseudowavefunction and spherical bessel function
        virtual void cal_ovlp_pswfcjlq() { ModuleBase::WARNING_QUIT("psi_initializer::calc_ovlp_pswfcjlq", "Polymorphism error"); }
        // nao
        /// @brief setter of numerical orbital files, useful when init_wfc = nao
        virtual void set_orbital_files(std::string* orbital_files) { ModuleBase::WARNING_QUIT("psi_initializer::set_orbital_files", "Polymorphism error"); }
        /// @brief calculate overlap between numerical orbital and spherical bessel function
        virtual void cal_ovlp_flzjlq() { ModuleBase::WARNING_QUIT("psi_initializer::cal_ovlp_flzjlq", "Polymorphism error"); }
        // atomic+random
        // nao+random
        /// @brief setter of random_mix
        /// @param random_mix_in new value of random_mix
        void set_random_mix(const double random_mix_in) { this->random_mix = random_mix_in; }
        /// @brief getter of random_mix
        /// @return this->random_mix
        double get_random_mix() const { return this->random_mix; }
        /// @brief getter of ixy2is, the mapping from fftixy to stick index
        /// @return this->ixy2is
        int* get_ixy2is() const { return this->ixy2is; }
        /// @brief setter of ixy2is, the mapping from fftixy to stick index
        void set_ixy2is(int* ixy2is_in) { this->ixy2is = ixy2is_in; }
        // member variables
        /// @brief interface to the psi::Psi data structure class
        psi::Psi<std::complex<double>>* psig;
        /// @brief interface to the Structure_Factor method class
        Structure_Factor* sf;
        /// @brief interface to the PW_Basis_K data structure class
        ModulePW::PW_Basis_K* pw_wfc;
        /// @brief method of Spherical Bessel Transformation
        ModuleBase::SphericalBesselTransformer sbt; // useful for atomic-like methods

    private:
        // basic properties
        int mem_saver = 0; // will deprecated this variable soon
        std::string method = "none";
        // non-random case
        int nbands_complem = 0;
        // random
        int* ixy2is;

        // atomic+random or nao+random
        double random_mix = 0;
};
#endif