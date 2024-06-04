#ifndef K_VECTORS_H
#define K_VECTORS_H
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

    int nmp[3];						// Number of Monhorst-Pack
    std::vector<int> kl_segids;	// index of kline segment

    K_Vectors();
    ~K_Vectors();

    void set(const ModuleSymmetry::Symmetry &symm,
             const std::string &k_file_name,
             const int& nspin,
             const ModuleBase::Matrix3 &reciprocal_vec,
             const ModuleBase::Matrix3 &latvec,
             std::ofstream&);

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
        { return kvec(ik, direct, irreducible); });
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

public:
    // Write the simple getter function into one line would be helpful for future maintainment
    int get_nks() const { return this->nks; }
    int get_nkstot() const { return this->nkstot; }
    int get_nkstot_ibz() const { return this->nkstot_ibz; }
    int get_nkstot_full() const { return this->nkstot_full; }

    // No I cannot understand why are setter functions here. Why need setter functions? who will need them?
    void set_nks(int value) { this->nks = value; }
    void set_nkstot(int value) { this->nkstot = value; }
    void set_nkstot_ibz(int value) { this->nkstot_ibz = value; }
    void set_nkstot_full(int value) { this->nkstot_full = value; }

private:
    int nks;						// number of k points in this pool(processor, up+dw)
    int nkstot;						/// total number of k points, equal to nkstot_ibz after reducing k points
    int nkstot_ibz;             /// number of k points in IBZ
    int nkstot_full;    /// number of k points in full k mesh
    
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
