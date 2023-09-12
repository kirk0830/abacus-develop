#include "psi_initializer.h"
#include "module_base/memory.h"

int psi_initializer::get_starting_nw() const
{
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
}

psi::Psi<std::complex<double>>* psi_initializer::allocate() {

	int prefactor = 1;
    int nbands_actual = (GlobalV::init_wfc.substr(0, 3) == "nao")? 
                            GlobalV::NLOCAL : GlobalV::NBANDS;
	int nkpts_actual = (GlobalV::CALCULATION == "nscf" && this->mem_saver == 1)? 
                            1 : this->pw_wfc->nks;
    psi::Psi<std::complex<double>>* psi_out = nullptr;
    psi_out = new psi::Psi<std::complex<double>>(
        nkpts_actual, 
            nbands_actual, 
                this->pw_wfc->npwk_max * GlobalV::NPOL, 
                    this->pw_wfc->npwk);
    const size_t memory_cost = 
        nkpts_actual*
            nkpts_actual*
                this->pw_wfc->npwk_max * GlobalV::NPOL*
                    sizeof(std::complex<double>);
	std::cout << " MEMORY FOR PSI (MB)  : " << double(memory_cost)/1024.0/1024.0 << std::endl;
	ModuleBase::Memory::record("Psi_PW", memory_cost);
    return psi_out;
}
