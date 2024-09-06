#include <vector>

#include "charge.h"
#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_base/libm/libm.h"
#include "module_base/math_integral.h"
#include "module_base/math_sphbes.h"
#include "module_base/memory.h"
#include "module_base/parallel_reduce.h"
#include "module_base/timer.h"
#include "module_base/tool_threading.h"
#include "module_elecstate/magnetism.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_hamilt_pw/hamilt_pwdft/parallel_grid.h"
#include "module_io/rho_io.h"
#include "module_io/rhog_io.h"
#ifdef USE_PAW
#include "module_cell/module_paw/paw_cell.h"
#endif
#include <numeric>
#include <algorithm>
void Charge::init_rho(elecstate::efermi& eferm_iout, const ModuleBase::ComplexMatrix& strucFac, const int& nbz, const int& bz)
{
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "init_chg", GlobalV::init_chg);

    std::cout << " START CHARGE      : " << GlobalV::init_chg << std::endl;
    bool read_error = false;
    if (GlobalV::init_chg == "file" || GlobalV::init_chg == "auto")
    {
        GlobalV::ofs_running << " try to read charge from file : " << std::endl;

        // try to read charge from binary file first, which is the same as QE
        // liuyu 2023-12-05
        std::stringstream binary;
        binary << GlobalV::global_readin_dir << PARAM.inp.suffix + "-CHARGE-DENSITY.restart";
        if (ModuleIO::read_rhog(binary.str(), rhopw, rhog))
        {
            GlobalV::ofs_running << " Read in the charge density: " << binary.str() << std::endl;
            // for (int is = 0; is < GlobalV::NSPIN; ++is)
            // {
            //     rhopw->recip2real(rhog[is], rho[is]);
            // }
            std::cout << " Nomenclature:\n read-in rho in realspace: rho(0), in recip.space: rhog(0)." << std::endl;
            std::cout << " Read rhog from " << binary.str() << std::endl;
        }
        if(true) // else
        {
            for (int is = 0; is < GlobalV::NSPIN; ++is)
            {
                std::stringstream ssc;
                ssc << GlobalV::global_readin_dir << "SPIN" << is + 1 << "_CHG.cube";
                std::cout << " Read rho from " << ssc.str() << std::endl;
                double& ef_tmp = eferm_iout.get_ef(is);
                if (ModuleIO::read_rho(
#ifdef __MPI
                        &(GlobalC::Pgrid),
#endif
                        GlobalV::MY_RANK,
                        GlobalV::ESOLVER_TYPE,
                        GlobalV::RANK_IN_STOGROUP,
                        is,
                        GlobalV::ofs_running,
                        GlobalV::NSPIN,
                        ssc.str(),
                        this->rho[is],
                        this->rhopw->nx,
                        this->rhopw->ny,
                        this->rhopw->nz,
                        ef_tmp,
                        &(GlobalC::ucell),
                        this->prenspin))
                {
                    GlobalV::ofs_running << " Read in the charge density: " << ssc.str() << std::endl;
                    double err_r = 0;
                    double err_g = 0;
                    std::vector<std::complex<double>> rhog_i(this->rhopw->npw);
                    std::vector<double> rho_i(this->rhopw->nrxx); // nrxx is the dimension of rho in present processor
                    std::vector<double> rho_temp(this->rhopw->nrxx);
                    std::vector<std::complex<double>> rhog_temp(this->rhopw->npw);
                    std::vector<double> rho1(this->rhopw->nrxx);
                    double* rho0 = this->rho[is];
                    std::complex<double>* rhog0 = this->rhog[is];

                    // charge integrity check...
                    double charge = std::accumulate(rho0, rho0 + this->rhopw->nrxx, 0.0) * GlobalC::ucell.omega / this->rhopw->nxyz;
                    #ifdef __MPI
                        Parallel_Reduce::reduce_pool(charge);
                        if (GlobalV::MY_RANK == 0)
                        {
                    #endif
                    printf("integrate rho(0): %f e\n", charge);
                    #ifdef __MPI
                        }
                    #endif

                    // rho(r) -> rho(G) -> rho(r)
                    this->rhopw->real2recip(rho0, rhog_i.data()); // -> rho(G)
                    this->rhopw->recip2real(rhog_i.data(), rho_i.data()); std::copy(rho_i.begin(), rho_i.end(), rho1.begin()); // save rho1(r)
                    // calculate error w.r.t. rho0
                    std::transform(rho_i.begin(), rho_i.end(), rho0, rho_temp.begin(), std::minus<double>());
                    err_r = std::accumulate(rho_temp.begin(), rho_temp.end(), 0.0, [](double sum, double val) { return sum + std::abs(val); });
                    std::transform(rhog_i.begin(), rhog_i.end(), rhog0, rhog_temp.begin(), std::minus<std::complex<double>>());
                    err_g = std::accumulate(rhog_temp.begin(), rhog_temp.end(), 0.0, [](double sum, std::complex<double> val) { return sum + std::abs(val); });
                    #ifdef __MPI
                        Parallel_Reduce::reduce_pool(err_r);
                        if (GlobalV::MY_RANK == 0)
                        {
                    #endif
                    printf("initialize rho(1) by rho(0) -> rho(G) -> rho(1).\n"
                           "sum(abs(err)) of rho(1) w.r.t. rho(0) = %16.10e\n", err_r);
                    printf("sum(abs(err)) of rhog(1) w.r.t. rhog(0) = %16.10e\n", err_g);
                    #ifdef __MPI
                        }
                    #endif
                    // more tests...
                    for (int itest = 0; itest < 8; itest++)
                    {
                        // rho(r) -> rho(G) -> rho(r)
                        this->rhopw->real2recip(rho_i.data(), rhog_i.data()); // -> rho(G)
                        this->rhopw->recip2real(rhog_i.data(), rho_i.data()); // -> rho(r)
                        charge = std::accumulate(rho_i.begin(), rho_i.end(), 0.0) * GlobalC::ucell.omega / this->rhopw->nxyz;
                        #ifdef __MPI
                            Parallel_Reduce::reduce_pool(charge);
                            if (GlobalV::MY_RANK == 0)
                            {
                        #endif
                        printf("Test #%d:\n", itest+1);
                        printf("integrate rho(%d): %f e\n", itest+2, charge);
                        #ifdef __MPI
                            }
                        #endif
                        // calculate error w.r.t. rho0
                        std::transform(rho_i.begin(), rho_i.end(), rho0, rho_temp.begin(), std::minus<double>());
                        err_r = std::accumulate(rho_temp.begin(), rho_temp.end(), 0.0, [](double sum, double val) { return sum + std::abs(val); });
                        std::transform(rhog_i.begin(), rhog_i.end(), rhog0, rhog_temp.begin(), std::minus<std::complex<double>>());
                        err_g = std::accumulate(rhog_temp.begin(), rhog_temp.end(), 0.0, [](double sum, std::complex<double> val) { return sum + std::abs(val); });
                        #ifdef __MPI
                            Parallel_Reduce::reduce_pool(err_r);
                            if (GlobalV::MY_RANK == 0)
                            {
                        #endif
                        printf("sum(abs(err)) of rho(%d) w.r.t. rho(0) & rho(1) = %16.10e, ", itest+2, err_r);
                        #ifdef __MPI
                            }
                        #endif
                        // then calculate error w.r.t. rho1
                        std::transform(rho_i.begin(), rho_i.end(), rho1.begin(), rho_temp.begin(), std::minus<double>());
                        err_r = std::accumulate(rho_temp.begin(), rho_temp.end(), 0.0, [](double sum, double val) { return sum + std::abs(val); });
                        #ifdef __MPI
                            Parallel_Reduce::reduce_pool(err_r);
                            if (GlobalV::MY_RANK == 0)
                            {
                        #endif
                        printf("%16.10e\n", err_r);
                        printf("sum(abs(err)) of rhog(%d) w.r.t. rhog(0) = %16.10e\n\n", itest+2, err_g);
                        #ifdef __MPI
                            }
                        #endif
                    }
                    // overwrite rho0 with rho1
                    // std::copy(rho1.begin(), rho1.end(), rho0);
                    //ModuleBase::WARNING_QUIT("Charge::init_rho", "Test recip2real & real2recip error");
                }
                else if (is > 0)
                {
                    if (prenspin == 1)
                    {
                        GlobalV::ofs_running << " Didn't read in the charge density but autoset it for spin " << is + 1
                                             << std::endl;
                        for (int ir = 0; ir < this->rhopw->nrxx; ir++)
                        {
                            this->rho[is][ir] = 0.0;
                        }
                    }
                    //
                    else if (prenspin == 2)
                    { // read up and down , then rearrange them.
                        if (is == 1)
                        {
                            std::cout << "Incomplete charge density file!" << std::endl;
                            read_error = true;
                            break;
                        }
                        else if (is == 2)
                        {
                            GlobalV::ofs_running << " Didn't read in the charge density but would rearrange it later. "
                                                 << std::endl;
                        }
                        else if (is == 3)
                        {
                            GlobalV::ofs_running << " rearrange charge density " << std::endl;
                            for (int ir = 0; ir < this->rhopw->nrxx; ir++)
                            {
                                this->rho[3][ir] = this->rho[0][ir] - this->rho[1][ir];
                                this->rho[0][ir] = this->rho[0][ir] + this->rho[1][ir];
                                this->rho[1][ir] = 0.0;
                                this->rho[2][ir] = 0.0;
                            }
                        }
                    }
                }
                else
                {
                    read_error = true;
                    break;
                }
            }

            if (XC_Functional::get_func_type() == 3 || XC_Functional::get_func_type() == 5)
            {
                for (int is = 0; is < GlobalV::NSPIN; is++)
                {
                    std::stringstream ssc;
                    ssc << GlobalV::global_readin_dir << "SPIN" << is + 1 << "_TAU.cube";
                    GlobalV::ofs_running << " try to read kinetic energy density from file : " << ssc.str()
                                         << std::endl;
                    // mohan update 2012-02-10, sunliang update 2023-03-09
                    if (ModuleIO::read_rho(
#ifdef __MPI
                            &(GlobalC::Pgrid),
#endif
                            GlobalV::MY_RANK,
                            GlobalV::ESOLVER_TYPE,
                            GlobalV::RANK_IN_STOGROUP,
                            is,
                            GlobalV::ofs_running,
                            GlobalV::NSPIN,
                            ssc.str(),
                            this->kin_r[is],
                            this->rhopw->nx,
                            this->rhopw->ny,
                            this->rhopw->nz,
                            eferm_iout.ef,
                            &(GlobalC::ucell),
                            this->prenspin))
                    {
                        GlobalV::ofs_running << " Read in the kinetic energy density: " << ssc.str() << std::endl;
                    }
                }
            }
        }
    }
    if (read_error)
    {
        const std::string warn_msg = " WARNING: \"init_chg\" is enabled but ABACUS failed to read charge density from file.\n"
                                     " Please check if there is SPINX_CHG.cube (X=1,...) or {suffix}-CHARGE-DENSITY.restart in the directory.\n";
        std::cout << std::endl << warn_msg;
        if (GlobalV::init_chg == "auto")
        {
            std::cout << " Charge::init_rho: use atomic initialization instead." << std::endl << std::endl;
        }
        else if (GlobalV::init_chg == "file")
        {
            ModuleBase::WARNING_QUIT("Charge::init_rho", "Failed to read in charge density from file.\nIf you want to use atomic charge initialization, \nplease set init_chg to atomic in INPUT.");
        }
    }

    if (GlobalV::init_chg == "atomic" || 
        (GlobalV::init_chg == "auto" && read_error)) // mohan add 2007-10-17
    {
        this->atomic_rho(GlobalV::NSPIN, GlobalC::ucell.omega, rho, strucFac, GlobalC::ucell);

        // liuyu 2023-06-29 : move here from atomic_rho(), which will be called several times in charge extrapolation
        // wenfei 2021-7-29 : initial tau = 3/5 rho^2/3, Thomas-Fermi
        if (XC_Functional::get_func_type() == 3 || XC_Functional::get_func_type() == 5)
        {
            const double pi = 3.141592653589790;
            const double fact = (3.0 / 5.0) * pow(3.0 * pi * pi, 2.0 / 3.0);
            for (int is = 0; is < GlobalV::NSPIN; ++is)
            {
                for (int ir = 0; ir < this->rhopw->nrxx; ++ir)
                {
                    kin_r[is][ir] = fact * pow(std::abs(rho[is][ir]) * GlobalV::NSPIN, 5.0 / 3.0) / GlobalV::NSPIN;
                }
            }
        }
    }

    // Peize Lin add 2020.04.04
    if (GlobalC::restart.info_load.load_charge && !GlobalC::restart.info_load.load_charge_finish)
    {
        for (int is = 0; is < GlobalV::NSPIN; ++is)
        {
            GlobalC::restart.load_disk("charge", is, this->nrxx, rho[is]);
        }
        GlobalC::restart.info_load.load_charge_finish = true;
    }
#ifdef __MPI
    this->init_chgmpi(nbz, bz);
#endif
}

//==========================================================
// computes the core charge on the real space 3D mesh.
//==========================================================
void Charge::set_rho_core(
    const ModuleBase::ComplexMatrix &structure_factor
)
{
    ModuleBase::TITLE("Charge","set_rho_core");
    ModuleBase::timer::tick("Charge","set_rho_core");

    // double eps = 1.e-10;
    //----------------------------------------------------------
    // LOCAL VARIABLES :
    // counter on mesh points
    // counter on atomic types
    // counter on g vectors
    //----------------------------------------------------------
    // int ir = 0;
    // int it = 0;
    // int ig = 0;

    bool bl = false;
    for (int it = 0; it<GlobalC::ucell.ntype; it++)
    {
        if (GlobalC::ucell.atoms[it].ncpp.nlcc)
        {
            bl = true;
            break;
        }
    }

    if (!bl)
    {
        ModuleBase::GlobalFunc::ZEROS( this->rho_core, this->rhopw->nrxx);
    	ModuleBase::timer::tick("Charge","set_rho_core");
        return;
    }

    double *rhocg = new double[this->rhopw->ngg];
    ModuleBase::GlobalFunc::ZEROS(rhocg, this->rhopw->ngg );

	// three dimension.
    std::complex<double> *vg = new std::complex<double>[this->rhopw->npw];	

    for (int it = 0; it < GlobalC::ucell.ntype;it++)
    {
        if (GlobalC::ucell.atoms[it].ncpp.nlcc)
        {
//----------------------------------------------------------
// EXPLAIN : drhoc compute the radial fourier transform for
// each shell of g vec
//----------------------------------------------------------
            this->non_linear_core_correction(
                GlobalC::ppcell.numeric,
                GlobalC::ucell.atoms[it].ncpp.msh,
                GlobalC::ucell.atoms[it].ncpp.r.data(),
                GlobalC::ucell.atoms[it].ncpp.rab.data(),
                GlobalC::ucell.atoms[it].ncpp.rho_atc.data(),
                rhocg);
//----------------------------------------------------------
// EXPLAIN : multiply by the structure factor and sum
//----------------------------------------------------------
            for (int ig = 0; ig < this->rhopw->npw ; ig++)
            {
                vg[ig] += structure_factor(it, ig) * rhocg[this->rhopw->ig2igg[ig]];
            }
        }
    }

	// for tmp use.
	for(int ig=0; ig< this->rhopw->npw; ig++)
	{
		this->rhog_core[ig] = vg[ig];
	}

    this->rhopw->recip2real(vg, this->rho_core);

    // test on the charge and computation of the core energy
    double rhoima = 0.0;
    double rhoneg = 0.0;
    for (int ir = 0; ir < this->rhopw->nrxx; ir++)
    {
        rhoneg += std::min(0.0, this->rhopw->ft.get_auxr_data<double>()[ir].real());
        rhoima += std::abs(this->rhopw->ft.get_auxr_data<double>()[ir].imag());
        // NOTE: Core charge is computed in reciprocal space and brought to real
        // space by FFT. For non smooth core charges (or insufficient cut-off)
        // this may result in negative values in some grid points.
        // Up to October 1999 the core charge was forced to be positive definite.
        // This induces an error in the force, and probably stress, calculation if
        // the number of grid points where the core charge would be otherwise neg
        // is large. The error disappears for sufficiently high cut-off, but may be
        // rather large and it is better to leave the core charge as it is.
        // If you insist to have it positive definite (with the possible problems
        // mentioned above) uncomment the following lines.  SdG, Oct 15 1999
    }

	// mohan fix bug 2011-04-03
    Parallel_Reduce::reduce_pool(rhoneg);
    Parallel_Reduce::reduce_pool(rhoima);

	// mohan changed 2010-2-2, make this same as in atomic_rho.
	// still lack something......
    rhoneg /= this->rhopw->nxyz * GlobalC::ucell.omega;
    rhoima /= this->rhopw->nxyz * GlobalC::ucell.omega;

    // calculate core_only exch-corr energy etxcc=E_xc[rho_core] if required
    // The term was present in previous versions of the code but it shouldn't
    delete [] rhocg;
    delete [] vg;
    ModuleBase::timer::tick("Charge","set_rho_core");
    return;
} // end subroutine set_rhoc

void Charge::set_rho_core_paw()
{
    ModuleBase::TITLE("Charge","set_rho_core_paw");
#ifdef USE_PAW
    double* tmp = new double[nrxx];
    GlobalC::paw_cell.get_vloc_ncoret(tmp,this->rho_core);
    delete[] tmp;

    this->rhopw->real2recip(this->rho_core,this->rhog_core);
#endif
}

void Charge::non_linear_core_correction
(
    const bool &numeric,
    const int mesh,
    const double *r,
    const double *rab,
    const double *rhoc,
    double *rhocg) const
{
    ModuleBase::TITLE("charge","drhoc");

	// use labmda instead of repeating codes
	const auto kernel = [&](int num_threads, int thread_id)
	{

	double gx = 0.0;
    double rhocg1 = 0.0;
    double *aux;

    // here we compute the fourier transform is the charge in numeric form
    if (numeric)
    {
        aux = new double [mesh];
        // G=0 term

        int igl0 = 0;
        if (this->rhopw->gg_uniq [0] < 1.0e-8)
        {
			// single thread term
			if (thread_id == 0)
			{
				for (int ir = 0;ir < mesh; ir++)
				{
					aux [ir] = r [ir] * r [ir] * rhoc [ir];
				}
				ModuleBase::Integral::Simpson_Integral(mesh, aux, rab, rhocg1);
				//rhocg [1] = fpi * rhocg1 / omega;
				rhocg [0] = ModuleBase::FOUR_PI * rhocg1 / GlobalC::ucell.omega;//mohan modify 2008-01-19
			}
            igl0 = 1;
        }

		int igl_beg, igl_end;
		// exclude igl0
		ModuleBase::TASK_DIST_1D(num_threads, thread_id, this->rhopw->ngg - igl0, igl_beg, igl_end);
		igl_beg += igl0;
		igl_end += igl_beg;

        // G <> 0 term
        for (int igl = igl_beg; igl < igl_end;igl++) 
        {
            gx = sqrt(this->rhopw->gg_uniq[igl] * GlobalC::ucell.tpiba2);
            ModuleBase::Sphbes::Spherical_Bessel(mesh, r, gx, 0, aux);
            for (int ir = 0;ir < mesh; ir++) 
            {
                aux [ir] = r[ir] * r[ir] * rhoc [ir] * aux [ir];
            } //  enddo
            ModuleBase::Integral::Simpson_Integral(mesh, aux, rab, rhocg1);
            rhocg [igl] = ModuleBase::FOUR_PI * rhocg1 / GlobalC::ucell.omega;
        } //  enddo
        delete [] aux;
    }
    else
    {
        // here the case where the charge is in analytic form,
        // check old version before 2008-12-9
    }

	}; // end kernel

	// do not use omp parallel when this function is already in parallel block
	// 
	// it is called in parallel block in Forces::cal_force_cc,
	// but not in other funtcion such as Stress_Func::stress_cc.
	ModuleBase::TRY_OMP_PARALLEL(kernel);

    return;
}