#ifndef K_VECTORS_H
#define K_VECTORS_H

/*
    On the refactor to remove GlobalV and GlobalC from class implementation
    2024/03/30 Kirk0830

    What IS kpoint?
    kpoint remarks translational symmetry of a crystal wavefunction. Each kpoint represents a translational
    symmetry that the wavefunction can only strictly recover itself after a specific translation operation.
    A more modern understanding is the translational operator and Hamiltonian operator is commutable, thus
    the eigenstates of the translational operator can be used to linearly combined to form the eigenstates of
    Hamiltonian operator.

    For programming, what can kpoint class have?
    0. Basic
    kpoint class is for storing kpoint related information, certainly it should have functionalities to read
    and even write ABACUS KPT files. Therefore it is a class, if compactly impelemented, it should hold the
    kpoint coordinates data, and have methods for I/O the kpoint coordinates.
    Comparatively the old implementation let kpoint know about nspin, the number of spin channels. However
    it is, totally unaccpectable, because the translational symmetry has nothing related to the spin symmetry.

    1. Symmetry: kpoint reduction, irreducible Brillouin zone
    More specifically, the import of kpoints from external files, should be accompied with a kpoint processing,
    such as the reduction of kpoints by symmetry operations. 
    Then it comes to the question "How to determine and reduce the number of kpoints according to symmetry".
    In symmetry module, there are symmetrical operations instantiated according to the symmetry detected by
    program. Imposing these operations on kpoints, if any two kpoints are equivalent, then they are the same,
    and one of them should be removed.

    2. Parallelization
    At least for planewave, the parallelization on kpoints is natural. Prsent parallelization strategy is to
    distribute kpoints by simple % and // way.
    
    Regulations on implementing MPI related functions:
    1. no matter if it is really MPI environment, always keep the workflow strictly unchanged, and if it is
    non-MPI, then some functions are left empty but they are still be called in workflow. Which means it is
    bad to implement functions like:
    #ifdef __MPI
    void mpi_k()
    {
        // do something
    }
    #endif
    , it is recommended to implement like:
    void mpi_k()
    {
        #ifdef __MPI
        // do something
        #endif
    }

    Future demands: HFX q-grid, phonon q-grid
    Thus it needs a new name, might be "BrillouinZoneSamplingGenerator" or something.
    bz_sampl.h

    namespace bz_sampl{
        static int nprocs;
        static int iproc;

        std::tuple<std::vector<container::Tensor>, std::vector<double>> generate(const int& nspin,
                                                                                 const std::string& fkpt,
                                                                                 const int& irank,
                                                                                 const int& nrank);
    }

    // in esolver_ks.h
    std::vector<container::Tensor> kpoints;
    // in esolver_ks.cpp
    // let what to be static?

    result = bz_sampl::generate(nspin, fkpt, irank, nrank);
    kpoints = std::get<0>(result);
    kweights = std::get<1>(result);
*/

#include <vector>

#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_base/matrix3.h"
#include "module_cell/unitcell.h"

class K_Vectors
{
public:

    std::vector<ModuleBase::Vector3<double>> kvec_c;		// Cartesian coordinates of k points
    std::vector<ModuleBase::Vector3<double>> kvec_d;		// Direct coordinates of k points
    std::vector<ModuleBase::Vector3<double>> kvec_d_ibz;	// ibz Direct coordinates of k points

    std::vector<double> wk;						// wk, weight of k points
    std::vector<double> wk_ibz;					// ibz kpoint wk ,weight of k points

    std::vector<int> ngk;						// ngk, number of plane waves for each k point
    std::vector<int> isk;						// distinguish spin up and down k points
    std::vector<int> ibz2bz;					// mohan added 2009-05-18

    int nks;						// number of k points in this pool(processor, up+dw)
    int nkstot;						/// total number of k points, equal to nkstot_ibz after reducing k points
    int nkstot_ibz;             /// number of k points in IBZ
    int nkstot_full;    /// number of k points in full k mesh

    int nmp[3];						// Number of Monhorst-Pack
    std::vector<int> kl_segids;	// index of kline segment

    K_Vectors();
    ~K_Vectors();

    void set(const ModuleSymmetry::Symmetry &symm,
             const std::string &k_file_name,
             const int& nspin,
             const ModuleBase::Matrix3 &reciprocal_vec,
             const ModuleBase::Matrix3 &latvec);

    void ibz_kpoint(const ModuleSymmetry::Symmetry &symm, 
                    bool use_symm,std::string& skpt, 
                    const UnitCell& ucell, 
                    bool& match);
    //LiuXh add 20180515
    void set_after_vc(const ModuleSymmetry::Symmetry &symm,
                     const std::string &k_file_name,
                     const int& nspin,
                     const ModuleBase::Matrix3 &reciprocal_vec,
                     const ModuleBase::Matrix3 &latvec);
    //get global index for ik
    inline int getik_global(const int& ik) const;

// REFACTORING K_Vectors functions below >>>
private:
    // newly built functions on refactor to remove GlobalC::Pkpoints and parallel_kpoints class
    K_Vectors(const std::string& fkpt,
              const double& lat0,
              const std::vector<std::vector<double>>& avecs,
              const std::vector<std::vector<double>>& bvecs,
              const int& npools, 
              const int& nprocs_inpool);

    // partial alternative to set() function
    bool build_kpt(const std::string& fkpt);

    // I/O functions
    // read ABACUS KPT file
    void read_abacus_kpt(const std::string& fkpt);
    // from special kpoints specified in KPT file, generate a kpoint-path
    static void interpolate_knodes(const std::vector<std::vector<double>>& knodes, //< [in] special kpoints direct coordinates
                                   const std::vector<int>& nks,                    //< [in] number of kpoints in each segment
                                   std::vector<ModuleBase::Vector3<double>>& kvec, //< [out] kpoints direct coordinates
                                   std::vector<int>& kseg_ids,                     //< [out] segment ids of kpoints
                                   int& nkstot,                                    //< [out] total number of kpoints
                                   std::vector<double>& wk);                       //< [out] weight of kpoints
    // convert kspacing to exact Monkhorst-Pack mesh dimensions
    static std::vector<int> kspacing_tompmesh(const std::vector<std::vector<double>>& bvecs,    //< [in] reciprocal lattice vectors
                                       const double& lat0,                               //< [in] lattice_constant
                                       const std::vector<double>& kspacing);             //< [in] kspacing
    // overloaded version for ModuleBase::Matrix3
    static std::vector<int> kspacing_tompmesh(const ModuleBase::Matrix3& bmat,                  //< [in] reciprocal lattice matrix
                                       const double& lat0,                               //< [in] lattice_constant
                                       const std::vector<double>& kspacing);             //< [in] kspacing
    // write Monkhorst-Pack kpoints to KPT file
    static std::string write_abacus_mpkmesh(const std::string& center,         //< [in] can be "Gamma" or "Monkhorst-Pack"
                                     const std::vector<int>& nmp,       //< [in] number of Monkhorst-Pack kpoints
                                     const std::vector<int>& shifts);   //< [in] shifts of Monkhorst-Pack kpoints
    // write Line mode kpoints to KPT file
    static std::string write_abacus_kline(const std::string& scale,                             //< [in] can be "Direct" or "Cartesian"
                                          const std::vector<std::vector<double>>& kvec,         //< [in] kpoints coordinates
                                          const std::vector<int>& nks);                         //< [in] number of kpoints in each segment
    // overloaded version for ModuleBase::Vector3, which is not safe and has partially been deprecated by new compilers like icpx
    static std::string write_abacus_kline(const std::string& scale,                             //< [in] can be "Direct" or "Cartesian"
                                          const std::vector<ModuleBase::Vector3<double>>& kvec, //< [in] kpoints coordinates
                                          const std::vector<int>& nks);                         //< [in] number of kpoints in each segment

    // Synchronize
    // synchronize between alpha and beta spin. Because kpoint does not distinguish spin physically, this way
    // is just a computational implementation convention, rather than physically strict/correct way.
    void sync_kvec_betweenspin(const int& nspin); // because in principle there is forever 1 set of kpoints,
                                                  // the special treatment like treating different spin
                                                  // kpoints separately is not physically correct. Thus the
                                                  // nspin is set as one parameter instead of a member variable.
    // synchronize between kvec_c and kvec_d
    void sync_kvec_betweencd(const bool& direct,            //< [in] true is from kvec_d to kvec_c, false is from kvec_c to kvec_d
                             const ModuleBase::Matrix3& t); //< [in] transformation matrix
    // synchronize between MPI processes
    void sync_kvec_betweenproc();

    std::vector<std::vector<double>> kvec(const std::vector<int>& iks,
                                          const bool& direct = true,
                                          const bool& irreducible = true) const
    {
        std::vector<std::vector<double>> kvecs(iks.size());
        std::transform(iks.begin(), iks.end(), kvecs.begin(), [&](const int& ik)
        {
            return kvec(ik, direct, irreducible);
        });
        return kvecs;
    }
    std::vector<double> kvec(const int& ik,
                             const bool& direct = true,
                             const bool& irreducible = true) const
    {
        if((direct)&&(irreducible)) return kvec_d_[ikibz2ik_[ik]];
        if((direct)&&(!irreducible)) return kvec_d_[ik];
        if((!direct)&&(irreducible)) return kvec_c_[ikibz2ik_[ik]];
        if((!direct)&&(!irreducible)) return kvec_c_[ik];
        return std::vector<double>();
    }
    std::vector<std::vector<double>> kvec_c_;
    std::vector<std::vector<double>> kvec_d_;
    std::vector<int> ik2ikibz_; // mapping from kpoint index to irreducible kpoint index
    std::vector<int> ikibz2ik_; // mapping from irreducible kpoint index to kpoint index
// <<< REFACTORING K_Vectors functions above

private:
    int nspin;
    bool kc_done;
    bool kd_done;
    double koffset[3];     			// used only in automatic k-points.
    std::string k_kword; //LiuXh add 20180619
    int k_nkstot; //LiuXh add 20180619
    bool is_mp = false; //Monkhorst-Pack
    
    void renew( const int &kpoint_number );

    // step 1 : generate kpoints
    bool read_kpoints(const std::string &fn); // return 0: something wrong.
    void monkhorst_pack(const int *nmp_in, 
                        const double *koffset_in,
                        const int tipo);
    double Monkhorst_Pack_formula(const int &k_type, 
                                  const double &offset,
                                  const int& n, 
                                  const int &dim);

    // step 2 : set both kvec and kved; normalize weight
    void update_use_ibz( void );
    void set_both_kvec(const ModuleBase::Matrix3 &G,const ModuleBase::Matrix3 &R, std::string& skpt);
    void normalize_wk( const int &degspin );

    // step 3 : mpi kpoints information.
    void mpi_k();

    // step 4 : *2 or *4 kpoints.
    // *2 for LSDA
    // *4 for non-collinear
    void set_kup_and_kdw();

    // step 5
    // print k lists.
    void print_klists(std::ofstream &fn);
    //bool read_kpoints_after_vc(const std::string &fn); //LiuXh add 20180515
    //void Monkhorst_Pack_after_vc(const int *nmp_in,const double *koffset_in,const int tipo); //LiuXh add 20180515
    void mpi_k_after_vc(); //LiuXh add 20180515 //Useless now, it should be removed after several versions' testing.
    void set_both_kvec_after_vc(const ModuleBase::Matrix3 &G,const ModuleBase::Matrix3 &R); //Useless now, it should be removed after several versions' testing.
    void set_kup_and_kdw_after_vc();
};

inline int K_Vectors:: getik_global(const int& ik) const
{
    int nkp = this->nkstot / GlobalV::KPAR;
    int rem = this->nkstot % GlobalV::KPAR;
    if(GlobalV::MY_POOL < rem)
    {
        return GlobalV::MY_POOL*nkp + GlobalV::MY_POOL + ik;
    }
    else
    {
        return GlobalV::MY_POOL*nkp + rem + ik;       
    }
}

#endif // KVECT_H
