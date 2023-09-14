#include "psi_initializer.h"
#include "module_base/memory.h"


psi_initializer::psi_initializer(Structure_Factor* sf_in, ModulePW::PW_Basis_K* pw_wfc_in): sf(sf_in), pw_wfc(pw_wfc_in)
{
    this->ixy2is = new int[this->pw_wfc->fftnxy];
    this->pw_wfc->getfftixy2is(this->ixy2is);
}


psi_initializer::~psi_initializer()
{
    delete[] this->ixy2is;
    if (this->psig != nullptr) delete this->psig;
}


int psi_initializer::get_starting_nw() const
{
    ModuleBase::timer::tick("psi_initializer", "get_starting_nw");
    if (GlobalV::init_wfc == "file")
    {
        throw std::runtime_error("wavefunc::get_starting_nw. init_wfc file is not implemented yet! "
                                 + ModuleBase::GlobalFunc::TO_STRING(__FILE__) + " line "
                                 + ModuleBase::GlobalFunc::TO_STRING(__LINE__)); // Peize Lin change 2019-05-01
        //**********************************************************************
        // ... read the wavefunction into memory (if it is not done in c_bands)
        //**********************************************************************
    }
    else if (GlobalV::init_wfc.substr(0,6) == "atomic")
    {
        if (GlobalC::ucell.natomwfc >= GlobalV::NBANDS)
        {
            if(GlobalV::test_wf) GlobalV::ofs_running << " Start wave functions are all pseudo atomic wave functions." << std::endl;
        }
        else
        {
            if(GlobalV::test_wf) GlobalV::ofs_running<<" Start wave functions are atomic + "
                                                     <<GlobalV::NBANDS - GlobalC::ucell.natomwfc
                                                     <<" random wavefunctions."<< std::endl;
        }
        return std::max(GlobalC::ucell.natomwfc,  GlobalV::NBANDS);
    }
    else if (GlobalV::init_wfc == "random")
    {
        if(GlobalV::test_wf) GlobalV::ofs_running << " Start wave functions are all random." << std::endl;
        return GlobalV::NBANDS;
    }
    else
    {
		throw std::runtime_error("wavefunc::get_starting_nw. Bad init_wfc parameter set in input file! "
                                + ModuleBase::GlobalFunc::TO_STRING(__FILE__) + " line "
                                + ModuleBase::GlobalFunc::TO_STRING(__LINE__)); // Peize Lin change 2019-05-01
    }
    ModuleBase::timer::tick("psi_initializer", "get_starting_nw");
}


psi::Psi<std::complex<double>>* psi_initializer::allocate()
{
    ModuleBase::timer::tick("psi_initializer", "allocate");
	int prefactor = 1;
    /*
    Here the method for determining memory space (number of band, the first dimension) of psi is:
    1. if nbands > nlocal and natomwfc, will use nbands
    2. if nbands < nlocal or natomwfc, will use nlocal or natomwfc
    */
    int nbands_files = (GlobalC::ucell.natomwfc >= GlobalV::NLOCAL)? 
                            GlobalC::ucell.natomwfc : GlobalV::NLOCAL;
    int nbands_actual = 0;
    if(GlobalV::init_wfc == "random") 
    {
        nbands_actual = GlobalV::NBANDS;
        this->nbands_complem = 0;
    }
    else
    {
        if(nbands_files >= GlobalV::NBANDS)
        {
            nbands_actual = nbands_files;
            this->nbands_complem = nbands_files - GlobalV::NBANDS;
        }
        else
        {
            nbands_actual = GlobalV::NBANDS;
            this->nbands_complem = 0;
        }
    }
	int nkpts_actual = (GlobalV::CALCULATION == "nscf" && this->mem_saver == 1)? 
                            1 : this->pw_wfc->nks;
    int nbasis_actual = this->pw_wfc->npwk_max * GlobalV::NPOL;
    psi::Psi<std::complex<double>>* psi_out = nullptr;
    psi_out = new psi::Psi<std::complex<double>>(
        nkpts_actual, 
            GlobalV::NBANDS, // because no matter what, the wavefunction finally needed has GlobalV::NBANDS bands
                nbasis_actual, 
                    this->pw_wfc->npwk);
    this->psig = new psi::Psi<std::complex<double>>(
        nkpts_actual, 
            nbands_actual, 
                nbasis_actual, 
                    this->pw_wfc->npwk);
    std::cout << "TEST: psi allocated memory for: " 
              << nkpts_actual << "*" 
              << nbands_actual << "*" 
              << nbasis_actual << std::endl;
    const size_t memory_cost = 
        nkpts_actual*
            nkpts_actual*
                this->pw_wfc->npwk_max * GlobalV::NPOL*
                    sizeof(std::complex<double>);
	std::cout << " MEMORY FOR PSI (MB)  : " << double(memory_cost)/1024.0/1024.0 << std::endl;
	ModuleBase::Memory::record("Psi_PW", memory_cost);
    ModuleBase::timer::tick("psi_initializer", "allocate");
    return psi_out;
}


void psi_initializer::print_status(psi::Psi<std::complex<double>>& psi) const
{
    std::cout << "Current method: " << this->method << std::endl;
    std::cout << "Psi status:" << std::endl;
    std::cout << "  number of kpoints: " << psi.get_nk() << std::endl;
    std::cout << "  number of bands: " << psi.get_nbands() << std::endl;
    std::cout << "  number of planewaves: " << psi.get_nbasis() << std::endl;
}

#ifdef __MPI
void psi_initializer::stick_to_pool(double* stick, const int& ir, double* out, const ModulePW::PW_Basis_K* wfc_basis) const
{	
    ModuleBase::timer::tick("psi_initializer", "stick_to_pool");
	MPI_Status ierror;
    const int is = this->ixy2is[ir];
	const int ip = wfc_basis->fftixy2ip[ir];
    const int nz = wfc_basis->nz;

	if(ip == 0 && GlobalV::RANK_IN_POOL ==0)
	{
		for(int iz=0; iz<nz; iz++)
		{
			out[is*nz+iz] = stick[iz];
		}
	}
	else if(ip == GlobalV::RANK_IN_POOL )
	{
		MPI_Recv(stick, nz, MPI_DOUBLE, 0, ir, POOL_WORLD,&ierror);
		for(int iz=0; iz<nz; iz++)
		{
			out[is*nz+iz] = stick[iz];
		}
	}
	else if(GlobalV::RANK_IN_POOL==0)
	{
		MPI_Send(stick, nz, MPI_DOUBLE, ip, ir, POOL_WORLD);
	}

	return;	
    ModuleBase::timer::tick("psi_initializer", "stick_to_pool");
}
#endif
