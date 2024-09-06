#include "module_io/input_conv.h"

#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_cell/module_symmetry/symmetry.h"
#include "module_cell/unitcell.h"
#include "module_elecstate/occupy.h"
#include "module_hamilt_general/module_surchem/surchem.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_io/berryphase.h"
#include "module_parameter/parameter.h"
#include "module_relax/relax_old/ions_move_basic.h"
#include "module_relax/relax_old/lattice_change_basic.h"

#include <algorithm>

#ifdef __EXX
#include "module_ri/exx_abfs-jle.h"
#endif

#ifdef __LCAO
#include "module_basis/module_ao/ORB_read.h"
#include "module_elecstate/potentials/H_TDDFT_pw.h"
#include "module_hamilt_lcao/hamilt_lcaodft/FORCE_STRESS.h"
#include "module_hamilt_lcao/module_dftu/dftu.h"
#include "module_hamilt_lcao/module_tddft/evolve_elec.h"
#include "module_hamilt_lcao/module_tddft/td_velocity.h"
#endif
#ifdef __PEXSI
#include "module_hsolver/module_pexsi/pexsi_solver.h"
#endif
#ifdef __MPI
#include "module_hsolver/diago_elpa.h"
#include "module_hsolver/diago_elpa_native.h"
#endif

#include "module_base/module_device/device.h"
#include "module_base/timer.h"
#include "module_elecstate/elecstate_lcao.h"
#include "module_elecstate/potentials/efield.h"
#include "module_elecstate/potentials/gatefield.h"
#include "module_hsolver/hsolver_lcao.h"
#include "module_hsolver/hsolver_pw.h"
#include "module_md/md_func.h"

#ifdef __LCAO
std::vector<double> Input_Conv::convert_units(std::string params, double c) {
    std::vector<double> params_ori;
    std::vector<double> params_out;
    parse_expression(params, params_ori);
    for (auto param: params_ori)
        params_out.emplace_back(param * c);

    return params_out;
}

void Input_Conv::read_td_efield()
{
    elecstate::H_TDDFT_pw::stype = PARAM.inp.td_stype;
    if (PARAM.inp.esolver_type == "tddft" && elecstate::H_TDDFT_pw::stype == 1)
    {
        TD_Velocity::tddft_velocity = true;
    } else {
        TD_Velocity::tddft_velocity = false;
    }
    if (PARAM.inp.out_mat_hs2 == 1)
    {
        TD_Velocity::out_mat_R = true;
    } else {
        TD_Velocity::out_mat_R = false;
    }
    parse_expression(PARAM.inp.td_ttype, elecstate::H_TDDFT_pw::ttype);

    elecstate::H_TDDFT_pw::tstart = PARAM.inp.td_tstart;
    elecstate::H_TDDFT_pw::tend = PARAM.inp.td_tend;

    elecstate::H_TDDFT_pw::dt = PARAM.mdp.md_dt / ModuleBase::AU_to_FS;
    elecstate::H_TDDFT_pw::dt_int = elecstate::H_TDDFT_pw::dt;

    // space domain parameters

    // length gauge
    elecstate::H_TDDFT_pw::lcut1 = PARAM.inp.td_lcut1;
    elecstate::H_TDDFT_pw::lcut2 = PARAM.inp.td_lcut2;

    // time domain parameters

    // Gauss
    elecstate::H_TDDFT_pw::gauss_omega = convert_units(PARAM.inp.td_gauss_freq,
                                                       2 * ModuleBase::PI * ModuleBase::AU_to_FS); // time(a.u.)^-1
    elecstate::H_TDDFT_pw::gauss_phase = convert_units(PARAM.inp.td_gauss_phase, 1.0);
    elecstate::H_TDDFT_pw::gauss_sigma = convert_units(PARAM.inp.td_gauss_sigma, 1 / ModuleBase::AU_to_FS);
    elecstate::H_TDDFT_pw::gauss_t0 = convert_units(PARAM.inp.td_gauss_t0, 1.0);
    elecstate::H_TDDFT_pw::gauss_amp = convert_units(PARAM.inp.td_gauss_amp,
                                                     ModuleBase::BOHR_TO_A / ModuleBase::Ry_to_eV); // Ry/bohr
    // init ncut for velocity gauge integral
    for (auto omega: elecstate::H_TDDFT_pw::gauss_omega) {
        int ncut
            = int(100.0 * omega * elecstate::H_TDDFT_pw::dt / ModuleBase::PI);
        if (ncut % 2 == 0) {
            ncut += 2;
        } else {
            ncut += 1;
        }
        if (elecstate::H_TDDFT_pw::stype == 0)
            ncut = 1;
        elecstate::H_TDDFT_pw::gauss_ncut.push_back(ncut);
    }
    // trapezoid
    elecstate::H_TDDFT_pw::trape_omega = convert_units(PARAM.inp.td_trape_freq,
                                                       2 * ModuleBase::PI * ModuleBase::AU_to_FS); // time(a.u.)^-1
    elecstate::H_TDDFT_pw::trape_phase = convert_units(PARAM.inp.td_trape_phase, 1.0);
    elecstate::H_TDDFT_pw::trape_t1 = convert_units(PARAM.inp.td_trape_t1, 1.0);
    elecstate::H_TDDFT_pw::trape_t2 = convert_units(PARAM.inp.td_trape_t2, 1.0);
    elecstate::H_TDDFT_pw::trape_t3 = convert_units(PARAM.inp.td_trape_t3, 1.0);
    elecstate::H_TDDFT_pw::trape_amp = convert_units(PARAM.inp.td_trape_amp,
                                                     ModuleBase::BOHR_TO_A / ModuleBase::Ry_to_eV); // Ry/bohr
    // init ncut for velocity gauge integral
    for (auto omega: elecstate::H_TDDFT_pw::trape_omega) {
        int ncut
            = int(100.0 * omega * elecstate::H_TDDFT_pw::dt / ModuleBase::PI);
        if (ncut % 2 == 0) {
            ncut += 2;
        } else {
            ncut += 1;
        }
        if (elecstate::H_TDDFT_pw::stype == 0)
            ncut = 1;
        elecstate::H_TDDFT_pw::trape_ncut.push_back(ncut);
    }
    // Trigonometric
    elecstate::H_TDDFT_pw::trigo_omega1 = convert_units(PARAM.inp.td_trigo_freq1,
                                                        2 * ModuleBase::PI * ModuleBase::AU_to_FS); // time(a.u.)^-1
    elecstate::H_TDDFT_pw::trigo_omega2 = convert_units(PARAM.inp.td_trigo_freq2,
                                                        2 * ModuleBase::PI * ModuleBase::AU_to_FS); // time(a.u.)^-1
    elecstate::H_TDDFT_pw::trigo_phase1 = convert_units(PARAM.inp.td_trigo_phase1, 1.0);
    elecstate::H_TDDFT_pw::trigo_phase2 = convert_units(PARAM.inp.td_trigo_phase2, 1.0);
    elecstate::H_TDDFT_pw::trigo_amp = convert_units(PARAM.inp.td_trigo_amp,
                                                     ModuleBase::BOHR_TO_A / ModuleBase::Ry_to_eV); // Ry/bohr
    // init ncut for velocity gauge integral
    for (auto omega: elecstate::H_TDDFT_pw::trigo_omega1) {
        int ncut
            = int(100.0 * omega * elecstate::H_TDDFT_pw::dt / ModuleBase::PI);
        if (ncut % 2 == 0) {
            ncut += 2;
        } else {
            ncut += 1;
        }
        if (elecstate::H_TDDFT_pw::stype == 0)
            ncut = 1;
        elecstate::H_TDDFT_pw::trigo_ncut.push_back(ncut);
    }
    // Heaviside
    elecstate::H_TDDFT_pw::heavi_t0 = convert_units(PARAM.inp.td_heavi_t0, 1.0);
    elecstate::H_TDDFT_pw::heavi_amp = convert_units(PARAM.inp.td_heavi_amp,
                                                     ModuleBase::BOHR_TO_A / ModuleBase::Ry_to_eV); // Ry/bohr

    return;
}
#endif

void Input_Conv::Convert()
{
    ModuleBase::TITLE("Input_Conv", "Convert");
    ModuleBase::timer::tick("Input_Conv", "Convert");
    GlobalV::CALCULATION = PARAM.globalv.global_calculation;
    GlobalV::double_grid = PARAM.globalv.double_grid;
    //-----------------------------------------------
    // set read_file_dir
    //-----------------------------------------------
    if (PARAM.inp.read_file_dir == "auto")
    {
        GlobalV::global_readin_dir = GlobalV::global_out_dir;
    }
    else
    {
        GlobalV::global_readin_dir = PARAM.inp.read_file_dir + '/';
    }
    //----------------------------------------------------------
    // main parameters / electrons / spin ( 10/16 )
    //----------------------------------------------------------
    //  suffix
    if (PARAM.inp.calculation == "md" && PARAM.mdp.md_restart) // md restart  liuyu add 2023-04-12
    {
        int istep = 0;
        double temperature = 0.0;
        MD_func::current_md_info(GlobalV::MY_RANK, GlobalV::global_readin_dir, istep, temperature);
        if (PARAM.inp.read_file_dir == "auto")
        {
            GlobalV::stru_file = GlobalV::global_stru_dir + "STRU_MD_" + std::to_string(istep);
        }
    } else if (PARAM.inp.stru_file != "") {
        GlobalV::stru_file = PARAM.inp.stru_file;
    }
    if (PARAM.inp.kpoint_file != "")
    {
        GlobalV::global_kpoint_card = PARAM.inp.kpoint_file;
    }
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "pseudo_dir", PARAM.inp.pseudo_dir);
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "orbital_dir", PARAM.inp.orbital_dir);
    // GlobalV::global_pseudo_type = PARAM.inp.pseudo_type;
    GlobalC::ucell.setup(PARAM.inp.latname,
                         PARAM.inp.ntype,
                         PARAM.inp.lmaxmax,
                         PARAM.inp.init_vel,
                         PARAM.inp.fixed_axes);

    if (PARAM.inp.calculation == "relax" || PARAM.inp.calculation == "cell-relax")
    {
        GlobalV::fixed_atoms = PARAM.inp.fixed_atoms;
    }


    GlobalV::NBANDS = PARAM.inp.nbands;

    GlobalV::device_flag = base_device::information::get_device_flag(PARAM.inp.device,
                                                                     PARAM.inp.ks_solver,
                                                                     PARAM.inp.basis_type);

    if (GlobalV::device_flag == "gpu" && PARAM.inp.basis_type == "pw")
    {
        GlobalV::KPAR = base_device::information::get_device_kpar(PARAM.inp.kpar);
    }
#ifdef __LCAO
    else if (PARAM.inp.basis_type == "lcao") {
        /// GlobalV::KPAR_LCAO is used in LCAO diagonalization only
        GlobalV::KPAR_LCAO = PARAM.inp.kpar;
        /// all other parts of the code use GlobalV::KPAR = 1
        GlobalV::KPAR = 1;
    }
#endif
    else
    {
        GlobalV::KPAR = PARAM.inp.kpar;
        GlobalV::NSTOGROUP = PARAM.inp.bndpar;
    }
    GlobalV::precision_flag = PARAM.inp.precision;
    if (GlobalV::device_flag == "cpu" and GlobalV::precision_flag == "single")
    {
// cpu single precision is not supported while float_fftw lib is not available
#ifndef __ENABLE_FLOAT_FFTW
        ModuleBase::WARNING_QUIT(
            "Input_Conv",
            "Single precision with cpu is not supported while float_fftw lib is not available; \
            \n Please recompile with cmake flag \"-DENABLE_FLOAT_FFTW=ON\".\n");
#endif // __ENABLE_FLOAT_FFTW
    }
    GlobalV::CALCULATION = PARAM.inp.calculation;
    GlobalV::ESOLVER_TYPE = PARAM.inp.esolver_type;

    GlobalV::PSEUDORCUT = PARAM.inp.pseudo_rcut;

    GlobalV::DFT_FUNCTIONAL = PARAM.inp.dft_functional;
    GlobalV::NSPIN = PARAM.inp.nspin;

    GlobalV::CAL_FORCE = PARAM.inp.cal_force;
    GlobalV::FORCE_THR = PARAM.inp.force_thr;

#ifdef __LCAO
    Force_Stress_LCAO<double>::force_invalid_threshold_ev = PARAM.inp.force_thr_ev2;
    Force_Stress_LCAO<std::complex<double>>::force_invalid_threshold_ev = PARAM.inp.force_thr_ev2;
#endif

    BFGS_Basic::relax_bfgs_w1 = PARAM.inp.relax_bfgs_w1;
    BFGS_Basic::relax_bfgs_w2 = PARAM.inp.relax_bfgs_w2;

    Ions_Move_Basic::relax_bfgs_rmax = PARAM.inp.relax_bfgs_rmax;
    Ions_Move_Basic::relax_bfgs_rmin = PARAM.inp.relax_bfgs_rmin;
    Ions_Move_Basic::relax_bfgs_init = PARAM.inp.relax_bfgs_init;
    Ions_Move_Basic::out_stru = PARAM.inp.out_stru; // mohan add 2012-03-23
    Lattice_Change_Basic::fixed_axes = PARAM.inp.fixed_axes;

    GlobalV::CAL_STRESS = PARAM.inp.cal_stress;


    GlobalV::RELAX_METHOD = PARAM.inp.relax_method;
    GlobalV::relax_new = PARAM.inp.relax_new;

    GlobalV::use_paw = PARAM.inp.use_paw;

    GlobalV::OUT_LEVEL = PARAM.inp.out_level;
    Ions_Move_CG::RELAX_CG_THR = PARAM.inp.relax_cg_thr; // pengfei add 2013-09-09

    ModuleSymmetry::Symmetry::symm_flag = std::stoi(PARAM.inp.symmetry);
    ModuleSymmetry::Symmetry::symm_autoclose = PARAM.inp.symmetry_autoclose;
    GlobalV::BASIS_TYPE = PARAM.inp.basis_type;
    GlobalV::KS_SOLVER = PARAM.inp.ks_solver;
    GlobalV::SEARCH_RADIUS = PARAM.inp.search_radius;

    //----------------------------------------------------------
    // planewave (8/8)
    //----------------------------------------------------------
    GlobalV::GAMMA_ONLY_LOCAL = PARAM.globalv.gamma_only_local;

    //----------------------------------------------------------
    // diagonalization  (5/5)
    //----------------------------------------------------------
    GlobalV::PW_DIAG_NMAX = PARAM.inp.pw_diag_nmax;
    GlobalV::PW_DIAG_NDIM = PARAM.inp.pw_diag_ndim;

    GlobalV::PW_DIAG_THR = PARAM.inp.pw_diag_thr;
    GlobalV::NB2D = PARAM.inp.nb2d;
    GlobalV::TEST_FORCE = PARAM.inp.test_force;
    GlobalV::TEST_STRESS = PARAM.inp.test_stress;

    //----------------------------------------------------------
    // iteration (1/3)
    //----------------------------------------------------------
    GlobalV::SCF_THR_TYPE = PARAM.inp.scf_thr_type;

#ifdef __LCAO
    if (PARAM.inp.dft_plus_u)
    {
        GlobalV::dft_plus_u = PARAM.inp.dft_plus_u;
        GlobalC::dftu.Yukawa = PARAM.inp.yukawa_potential;
        GlobalC::dftu.omc = PARAM.inp.omc;
        GlobalC::dftu.orbital_corr = PARAM.inp.orbital_corr;
        GlobalC::dftu.uramping = PARAM.globalv.uramping;
        GlobalC::dftu.mixing_dftu = PARAM.inp.mixing_dftu;
        GlobalC::dftu.U = PARAM.globalv.hubbard_u;
        GlobalC::dftu.U0 = PARAM.globalv.hubbard_u;
        if (PARAM.globalv.uramping > 0.01)
        {
            ModuleBase::GlobalFunc::ZEROS(GlobalC::dftu.U.data(), GlobalC::ucell.ntype);
        }
    }
#endif
    //--------------------------------------------
    // added by zhengdy-soc
    //--------------------------------------------
    if (PARAM.inp.noncolin || PARAM.inp.lspinorb)
    {
        GlobalV::NSPIN = 4;
    }

    if (GlobalV::NSPIN == 4)
    {
        GlobalV::NONCOLIN = PARAM.inp.noncolin;
        // wavefunctions are spinors with 2 components
        GlobalV::NPOL = 2;
        // set the domag variable to make a spin-orbit calculation with zero
        // magnetization
        GlobalV::DOMAG = false;
        GlobalV::DOMAG_Z = true;
        GlobalV::LSPINORB = PARAM.inp.lspinorb;
        if (PARAM.globalv.gamma_only_local)
        {
            ModuleBase::WARNING_QUIT("input_conv",
                                     "nspin=4(soc or noncollinear-spin) does "
                                     "not support gamma only calculation");
        }
    } else {
        GlobalV::LSPINORB = false;
        GlobalV::NONCOLIN = false;
        GlobalV::DOMAG = false;
        GlobalV::DOMAG_Z = false;
        GlobalV::NPOL = 1;
    }

    //----------------------------------------------------------
    // Yu Liu add 2022-05-18
    //----------------------------------------------------------
    GlobalV::EFIELD_FLAG = PARAM.inp.efield_flag;
    GlobalV::DIP_COR_FLAG = PARAM.inp.dip_cor_flag;
    elecstate::Efield::efield_dir = PARAM.inp.efield_dir;
    elecstate::Efield::efield_pos_max = PARAM.inp.efield_pos_max;
    elecstate::Efield::efield_pos_dec = PARAM.inp.efield_pos_dec;
    elecstate::Efield::efield_amp = PARAM.inp.efield_amp;

    //----------------------------------------------------------
    // Yu Liu add 2022-09-13
    //----------------------------------------------------------
    GlobalV::GATE_FLAG = PARAM.inp.gate_flag;
    GlobalV::nelec = PARAM.inp.nelec;
    if (PARAM.globalv.two_fermi)
    {
        GlobalV::TWO_EFERMI = true;
        GlobalV::nupdown = PARAM.inp.nupdown;
    }
    elecstate::Gatefield::zgate = PARAM.inp.zgate;
    elecstate::Gatefield::relax = PARAM.inp.relax;
    elecstate::Gatefield::block = PARAM.inp.block;
    elecstate::Gatefield::block_down = PARAM.inp.block_down;
    elecstate::Gatefield::block_up = PARAM.inp.block_up;
    elecstate::Gatefield::block_height = PARAM.inp.block_height;

//----------------------------------------------------------
// Fuxiang He add 2016-10-26
//----------------------------------------------------------
#ifdef __LCAO
    module_tddft::Evolve_elec::td_force_dt = PARAM.inp.td_force_dt;
    module_tddft::Evolve_elec::td_vext = PARAM.inp.td_vext;
    if (module_tddft::Evolve_elec::td_vext)
    {
        parse_expression(PARAM.inp.td_vext_dire, module_tddft::Evolve_elec::td_vext_dire_case);
    }
    module_tddft::Evolve_elec::out_dipole = PARAM.inp.out_dipole;
    module_tddft::Evolve_elec::out_efield = PARAM.inp.out_efield;
    module_tddft::Evolve_elec::td_print_eij = PARAM.inp.td_print_eij;
    module_tddft::Evolve_elec::td_edm = PARAM.inp.td_edm;
    TD_Velocity::out_current = PARAM.inp.out_current;
    TD_Velocity::out_current_k = PARAM.inp.out_current_k;
    TD_Velocity::out_vecpot = PARAM.inp.out_vecpot;
    TD_Velocity::init_vecpot_file = PARAM.inp.init_vecpot_file;
    read_td_efield();
#endif

   

    //----------------------------------------------------------
    // about restart, // Peize Lin add 2020-04-04
    //----------------------------------------------------------
    if (PARAM.inp.restart_save)
    {
        std::string dft_functional_lower = PARAM.inp.dft_functional;
        std::transform(PARAM.inp.dft_functional.begin(),
                       PARAM.inp.dft_functional.end(),
                       dft_functional_lower.begin(),
                       tolower);
        GlobalC::restart.folder = GlobalV::global_readin_dir + "restart/";
        ModuleBase::GlobalFunc::MAKE_DIR(GlobalC::restart.folder);
        if (dft_functional_lower == "hf" || dft_functional_lower == "pbe0"
            || dft_functional_lower == "hse"
            || dft_functional_lower == "opt_orb"
            || dft_functional_lower == "scan0") {
            GlobalC::restart.info_save.save_charge = true;
            GlobalC::restart.info_save.save_H = true;
        } else {
            GlobalC::restart.info_save.save_charge = true;
        }
    }
    if (PARAM.inp.restart_load)
    {
        std::string dft_functional_lower = PARAM.inp.dft_functional;
        std::transform(PARAM.inp.dft_functional.begin(),
                       PARAM.inp.dft_functional.end(),
                       dft_functional_lower.begin(),
                       tolower);
        GlobalC::restart.folder = GlobalV::global_readin_dir + "restart/";
        if (dft_functional_lower == "hf" || dft_functional_lower == "pbe0"
            || dft_functional_lower == "hse"
            || dft_functional_lower == "opt_orb"
            || dft_functional_lower == "scan0") {
            GlobalC::restart.info_load.load_charge = true;
            GlobalC::restart.info_load.load_H = true;
        } else {
            GlobalC::restart.info_load.load_charge = true;
        }
    }

//----------------------------------------------------------
// about exx, Peize Lin add 2018-06-20
//----------------------------------------------------------
#ifdef __EXX
#ifdef __LCAO

    std::string dft_functional_lower = PARAM.inp.dft_functional;
    std::transform(PARAM.inp.dft_functional.begin(),
                   PARAM.inp.dft_functional.end(),
                   dft_functional_lower.begin(),
                   tolower);
    if (dft_functional_lower == "hf" || dft_functional_lower == "pbe0" || dft_functional_lower == "scan0")
    {
        GlobalC::exx_info.info_global.cal_exx = true;
        GlobalC::exx_info.info_global.ccp_type
            = Conv_Coulomb_Pot_K::Ccp_Type::Hf;
    } else if (dft_functional_lower == "hse") {
        GlobalC::exx_info.info_global.cal_exx = true;
        GlobalC::exx_info.info_global.ccp_type
            = Conv_Coulomb_Pot_K::Ccp_Type::Hse;
    } else if (dft_functional_lower == "opt_orb") {
        GlobalC::exx_info.info_global.cal_exx = false;
        Exx_Abfs::Jle::generate_matrix = true;
    } else {
        GlobalC::exx_info.info_global.cal_exx = false;
    }

    if (GlobalC::exx_info.info_global.cal_exx || Exx_Abfs::Jle::generate_matrix || PARAM.inp.rpa)
    {
        // EXX case, convert all EXX related variables
        // GlobalC::exx_info.info_global.cal_exx = true;
        GlobalC::exx_info.info_global.hybrid_alpha = std::stod(PARAM.inp.exx_hybrid_alpha);
        XC_Functional::get_hybrid_alpha(std::stod(PARAM.inp.exx_hybrid_alpha));
        GlobalC::exx_info.info_global.hse_omega = PARAM.inp.exx_hse_omega;
        GlobalC::exx_info.info_global.separate_loop = PARAM.inp.exx_separate_loop;
        GlobalC::exx_info.info_global.hybrid_step = PARAM.inp.exx_hybrid_step;
        GlobalC::exx_info.info_global.mixing_beta_for_loop1 = PARAM.inp.exx_mixing_beta;
        GlobalC::exx_info.info_lip.lambda = PARAM.inp.exx_lambda;

        GlobalC::exx_info.info_ri.real_number = std::stoi(PARAM.inp.exx_real_number);
        GlobalC::exx_info.info_ri.pca_threshold = PARAM.inp.exx_pca_threshold;
        GlobalC::exx_info.info_ri.C_threshold = PARAM.inp.exx_c_threshold;
        GlobalC::exx_info.info_ri.V_threshold = PARAM.inp.exx_v_threshold;
        GlobalC::exx_info.info_ri.dm_threshold = PARAM.inp.exx_dm_threshold;
        GlobalC::exx_info.info_ri.cauchy_threshold = PARAM.inp.exx_cauchy_threshold;
        GlobalC::exx_info.info_ri.C_grad_threshold = PARAM.inp.exx_c_grad_threshold;
        GlobalC::exx_info.info_ri.V_grad_threshold = PARAM.inp.exx_v_grad_threshold;
        GlobalC::exx_info.info_ri.C_grad_R_threshold = PARAM.inp.exx_c_grad_r_threshold;
        GlobalC::exx_info.info_ri.V_grad_R_threshold = PARAM.inp.exx_v_grad_r_threshold;
        GlobalC::exx_info.info_ri.cauchy_force_threshold = PARAM.inp.exx_cauchy_force_threshold;
        GlobalC::exx_info.info_ri.cauchy_stress_threshold = PARAM.inp.exx_cauchy_stress_threshold;
        GlobalC::exx_info.info_ri.ccp_rmesh_times = std::stod(PARAM.inp.exx_ccp_rmesh_times);

        Exx_Abfs::Jle::Lmax = PARAM.inp.exx_opt_orb_lmax;
        Exx_Abfs::Jle::Ecut_exx = PARAM.inp.exx_opt_orb_ecut;
        Exx_Abfs::Jle::tolerence = PARAM.inp.exx_opt_orb_tolerence;

        // EXX does not support symmetry=1
        if (PARAM.inp.calculation != "nscf" && PARAM.inp.symmetry == "1")
        {
            ModuleSymmetry::Symmetry::symm_flag = 0;
        }
    }
#endif                                                   // __LCAO
#endif                                                   // __EXX
    GlobalC::ppcell.cell_factor = PARAM.inp.cell_factor; // LiuXh add 20180619

    //----------------------------------------------------------
    // reset symmetry flag to avoid error
    //----------------------------------------------------------
    // In these case, symmetry should be reset to 0
    // efield does not support symmetry=1
    if (PARAM.inp.efield_flag && ModuleSymmetry::Symmetry::symm_flag == 1)
    {
        ModuleSymmetry::Symmetry::symm_flag = 0;
    }
    // In these case, inversion symmetry is also not allowed, symmetry should be
    // reset to -1
    if (GlobalV::LSPINORB)
    {
        ModuleSymmetry::Symmetry::symm_flag = -1;
    }
    // end of symmetry reset

    //----------------------------------------------------------
    // main parameters / electrons / spin ( 2/16 )
    //----------------------------------------------------------
    //	electrons::nelup = PARAM.inp.nelup;
    //	electrons::neldw = PARAM.inp.neldw;

    //----------------------------------------------------------
    // occupation (3/3)
    //----------------------------------------------------------
    std::string occupations = "smearing";
    Occupy::decision(occupations, PARAM.inp.smearing_method, PARAM.inp.smearing_sigma);

    //----------------------------------------------------------
    // iteration
    //----------------------------------------------------------
    GlobalV::SCF_NMAX = PARAM.inp.scf_nmax;

    //----------------------------------------------------------
    // wavefunction / charge / potential / (2/4)
    //----------------------------------------------------------
    GlobalV::init_chg = PARAM.inp.init_chg;
    GlobalV::init_wfc = PARAM.inp.init_wfc;
    GlobalV::psi_initializer = PARAM.inp.psi_initializer;
    GlobalV::chg_extrap = PARAM.inp.chg_extrap; // xiaohui modify 2015-02-01
    GlobalV::nelec = PARAM.inp.nelec;
    GlobalV::out_pot = PARAM.inp.out_pot;
    GlobalV::out_app_flag = PARAM.inp.out_app_flag;

#ifdef __LCAO
    hsolver::HSolverLCAO<double>::out_mat_hs = PARAM.inp.out_mat_hs;
    hsolver::HSolverLCAO<double>::out_mat_hsR = PARAM.inp.out_mat_hs2; // LiuXh add 2019-07-16
    hsolver::HSolverLCAO<double>::out_mat_t = PARAM.inp.out_mat_t;
    hsolver::HSolverLCAO<double>::out_mat_dh = PARAM.inp.out_mat_dh;
    hsolver::HSolverLCAO<std::complex<double>>::out_mat_hs = PARAM.inp.out_mat_hs;
    hsolver::HSolverLCAO<std::complex<double>>::out_mat_hsR = PARAM.inp.out_mat_hs2; // LiuXh add 2019-07-16
    hsolver::HSolverLCAO<std::complex<double>>::out_mat_t = PARAM.inp.out_mat_t;
    hsolver::HSolverLCAO<std::complex<double>>::out_mat_dh = PARAM.inp.out_mat_dh;
    if (GlobalV::GAMMA_ONLY_LOCAL)
    {
        elecstate::ElecStateLCAO<double>::out_wfc_lcao = PARAM.inp.out_wfc_lcao;
    }
    else if (!GlobalV::GAMMA_ONLY_LOCAL)
    {
        elecstate::ElecStateLCAO<std::complex<double>>::out_wfc_lcao = PARAM.inp.out_wfc_lcao;
    }
    if (PARAM.inp.calculation == "nscf" && !PARAM.inp.towannier90 && !PARAM.inp.berry_phase)
    {
        if (GlobalV::GAMMA_ONLY_LOCAL)
        {
            elecstate::ElecStateLCAO<double>::need_psi_grid = false;
        } else if (!GlobalV::GAMMA_ONLY_LOCAL) {
            elecstate::ElecStateLCAO<std::complex<double>>::need_psi_grid
                = false;
        }
    }
    if (PARAM.inp.calculation == "test_neighbour" && GlobalV::NPROC > 1)
    {
        ModuleBase::WARNING_QUIT("Input_conv", "test_neighbour must be done with 1 processor");
    }
#endif

    //----------------------------------------------------------
    // About LCAO
    //----------------------------------------------------------
    // mohan add 2021-04-16
    //	ORB.ecutwfc = PARAM.inp.lcao_ecut;
    //	ORB.dk = PARAM.inp.lcao_dk;
    //	ORB.dR = PARAM.inp.lcao_dr;
    //	ORB.Rmax = PARAM.inp.lcao_rmax;

    // mohan add 2021-02-16
    berryphase::berry_phase_flag = PARAM.inp.berry_phase;

//-----------------------------------------------
// caoyu add for DeePKS
//-----------------------------------------------
#ifdef __DEEPKS
    GlobalV::deepks_scf = PARAM.inp.deepks_scf;
    GlobalV::deepks_bandgap = PARAM.inp.deepks_bandgap; // QO added for bandgap label 2021-12-15
    GlobalV::deepks_v_delta = PARAM.inp.deepks_v_delta;
    GlobalV::deepks_out_labels = PARAM.inp.deepks_out_labels;
    GlobalV::deepks_equiv = PARAM.inp.deepks_equiv;

    if (GlobalV::deepks_equiv && GlobalV::deepks_bandgap) {
        ModuleBase::WARNING_QUIT(
            "Input_conv",
            "deepks_equiv and deepks_bandgap cannot be used together");
    }
    if (PARAM.inp.deepks_out_unittest)
    {
        GlobalV::deepks_out_labels = true;
        GlobalV::deepks_scf = true;
        if (GlobalV::NPROC > 1)
        {
            ModuleBase::WARNING_QUIT("Input_conv", "generate deepks unittest with only 1 processor");
        }
        if (GlobalV::CAL_FORCE != 1)
        {
            ModuleBase::WARNING_QUIT("Input_conv", "force is required in generating deepks unittest");
        }
        if (GlobalV::CAL_STRESS != 1)
        {
            ModuleBase::WARNING_QUIT("Input_conv", "stress is required in generating deepks unittest");
        }
    }
    if (GlobalV::deepks_scf || GlobalV::deepks_out_labels)
    {
        GlobalV::deepks_setorb = true;
    }
#else
    if (PARAM.inp.deepks_scf || PARAM.inp.deepks_out_labels || PARAM.inp.deepks_bandgap || PARAM.inp.deepks_v_delta)
    {
        ModuleBase::WARNING_QUIT("Input_conv", "please compile with DeePKS");
    }
#endif
    //-----------------------------------------------
    // sunml add for implicit solvation model
    //-----------------------------------------------
    GlobalV::imp_sol = PARAM.inp.imp_sol;
    GlobalV::eb_k = PARAM.inp.eb_k;

    //-----------------------------------------------
    // Deltaspin related parameters
    //-----------------------------------------------
    GlobalV::sc_thr = PARAM.inp.sc_thr;

    // mixing parameters
    GlobalV::MIXING_MODE = PARAM.inp.mixing_mode;
    GlobalV::MIXING_BETA = PARAM.inp.mixing_beta;
    GlobalV::MIXING_NDIM = PARAM.inp.mixing_ndim;
    GlobalV::MIXING_RESTART = PARAM.inp.mixing_restart;
    GlobalV::MIXING_GG0 = PARAM.inp.mixing_gg0;
    GlobalV::MIXING_BETA_MAG = PARAM.inp.mixing_beta_mag;
    GlobalV::MIXING_GG0_MAG = PARAM.inp.mixing_gg0_mag;
    GlobalV::MIXING_GG0_MIN = PARAM.inp.mixing_gg0_min;
    GlobalV::MIXING_ANGLE = PARAM.inp.mixing_angle;
    GlobalV::MIXING_TAU = PARAM.inp.mixing_tau;
    GlobalV::MIXING_DMR = PARAM.inp.mixing_dmr;

    //-----------------------------------------------
    // Quasiatomic Orbital analysis
    //-----------------------------------------------
    GlobalV::qo_thr = PARAM.inp.qo_thr;
    GlobalV::qo_screening_coeff = PARAM.inp.qo_screening_coeff;

    //-----------------------------------------------
    // PEXSI related parameters
    //-----------------------------------------------
#ifdef __PEXSI
    pexsi::PEXSI_Solver::pexsi_npole = PARAM.inp.pexsi_npole;
    pexsi::PEXSI_Solver::pexsi_inertia = PARAM.inp.pexsi_inertia;
    pexsi::PEXSI_Solver::pexsi_nmax = PARAM.inp.pexsi_nmax;
    // pexsi::PEXSI_Solver::pexsi_symbolic = PARAM.inp.pexsi_symbolic;
    pexsi::PEXSI_Solver::pexsi_comm = PARAM.inp.pexsi_comm;
    pexsi::PEXSI_Solver::pexsi_storage = PARAM.inp.pexsi_storage;
    pexsi::PEXSI_Solver::pexsi_ordering = PARAM.inp.pexsi_ordering;
    pexsi::PEXSI_Solver::pexsi_row_ordering = PARAM.inp.pexsi_row_ordering;
    pexsi::PEXSI_Solver::pexsi_nproc = PARAM.inp.pexsi_nproc;
    pexsi::PEXSI_Solver::pexsi_symm = PARAM.inp.pexsi_symm;
    pexsi::PEXSI_Solver::pexsi_trans = PARAM.inp.pexsi_trans;
    pexsi::PEXSI_Solver::pexsi_method = PARAM.inp.pexsi_method;
    pexsi::PEXSI_Solver::pexsi_nproc_pole = PARAM.inp.pexsi_nproc_pole;
    // pexsi::PEXSI_Solver::pexsi_spin = PARAM.inp.pexsi_spin;
    pexsi::PEXSI_Solver::pexsi_temp = PARAM.inp.pexsi_temp;
    pexsi::PEXSI_Solver::pexsi_gap = PARAM.inp.pexsi_gap;
    pexsi::PEXSI_Solver::pexsi_delta_e = PARAM.inp.pexsi_delta_e;
    pexsi::PEXSI_Solver::pexsi_mu_lower = PARAM.inp.pexsi_mu_lower;
    pexsi::PEXSI_Solver::pexsi_mu_upper = PARAM.inp.pexsi_mu_upper;
    pexsi::PEXSI_Solver::pexsi_mu = PARAM.inp.pexsi_mu;
    pexsi::PEXSI_Solver::pexsi_mu_thr = PARAM.inp.pexsi_mu_thr;
    pexsi::PEXSI_Solver::pexsi_mu_expand = PARAM.inp.pexsi_mu_expand;
    pexsi::PEXSI_Solver::pexsi_mu_guard = PARAM.inp.pexsi_mu_guard;
    pexsi::PEXSI_Solver::pexsi_elec_thr = PARAM.inp.pexsi_elec_thr;
    pexsi::PEXSI_Solver::pexsi_zero_thr = PARAM.inp.pexsi_zero_thr;
#endif

    // elpa related
#ifdef __MPI
    hsolver::DiagoElpa<std::complex<double>>::elpa_num_thread = PARAM.inp.elpa_num_thread;
    ;
    hsolver::DiagoElpa<double>::elpa_num_thread = PARAM.inp.elpa_num_thread;
    ;
    hsolver::DiagoElpaNative<std::complex<double>>::elpa_num_thread = PARAM.inp.elpa_num_thread;
    ;
    hsolver::DiagoElpaNative<double>::elpa_num_thread = PARAM.inp.elpa_num_thread;
    ;
#endif
    ModuleBase::timer::tick("Input_Conv", "Convert");
    return;
}
