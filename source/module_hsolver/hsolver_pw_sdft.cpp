#include "hsolver_pw_sdft.h"
#include "module_base/timer.h"
#include "module_base/global_function.h"
#include "src_pw/symmetry_rho.h"
#include "module_base/timer.h"
#include "module_base/tool_title.h"

//temporary
#include "src_pw/global.h"
namespace hsolver
{
    void HSolverPW_SDFT::solve(hamilt::Hamilt* pHamilt, 
                           psi::Psi<std::complex<double>>& psi, 
                           elecstate::ElecState* pes, 
                           Stochastic_WF& stowf,
                           const int iter,
                           const std::string method_in, 
                           const bool skip_charge)
    {
        ModuleBase::TITLE(this->classname, "solve");
        ModuleBase::timer::tick(this->classname, "solve");
        const int npwx = psi.get_nbasis();
        const int nbands = psi.get_nbands();
        const int nks = psi.get_nk();

        // prepare for the precondition of diagonalization
        this->precondition.resize(psi.get_nbasis());

        // select the method of diagonalization
        this->method = method_in;
        this->initpdiagh();

        // part of KSDFT to get KS orbitals
	    for (int ik = 0; ik < nks; ++ik)
	    {
		    if(nbands > 0 && GlobalV::MY_STOGROUP == 0)
		    {
                pHamilt->updateHk(ik);
                psi.fix_k(ik);

                // template add precondition calculating here
                update_precondition(precondition, ik, this->wfc_basis->npwk[ik]);
		    	/// solve eigenvector and eigenvalue for H(k)
                double* p_eigenvalues = &(pes->ekb(ik, 0));
                this->hamiltSolvePsiK(pHamilt, psi, p_eigenvalues);
		    }
		    else
		    {
		    	pHamilt->updateHk(ik);
                psi.fix_k(ik);
		    	//In fact, hm.hpw.init_k has been done in wf.wfcinit();
		    }

            // DiagoCG would keep 9*nbasis memory in cache during loop-k
            // it should be deleted before calculating charge
            if(this->method == "cg")
            {
                delete pdiagh;
                pdiagh = nullptr;
            }
            
		    stoiter.stohchi.current_ik = ik;
		
#ifdef __MPI
			if(nbands > 0)
			{
				MPI_Bcast(psi.get_pointer(), npwx*nbands*2, MPI_DOUBLE , 0, PARAPW_WORLD);
				MPI_Bcast(&(pes->ekb(ik, 0)), nbands, MPI_DOUBLE, 0, PARAPW_WORLD);
			}
#endif
			stoiter.orthog(ik,psi,stowf);
			stoiter.checkemm(ik,iter, stowf);	//check and reset emax & emin
		}

		for (int ik = 0;ik < nks;ik++)
		{
			//init k
			if(nks > 1) pHamilt->updateHk(ik);
			stoiter.stohchi.current_ik = ik;
			stoiter.calPn(ik, stowf);
		}

		stoiter.itermu(iter,pes);
		//(5) calculate new charge density 
		// calculate KS rho.
		if(nbands > 0)
		{
			pes->psiToRho(psi);
#ifdef __MPI
			MPI_Bcast(&pes->eband,1, MPI_DOUBLE, 0,PARAPW_WORLD);
#endif
		}
		else
		{
			for(int is=0; is < GlobalV::NSPIN; is++)
			{
				ModuleBase::GlobalFunc::ZEROS(pes->charge->rho[is], GlobalC::rhopw->nrxx);
			}
		}
		// calculate stochastic rho
		stoiter.sum_stoband(stowf,pes);
		

		//(6) calculate the delta_harris energy 
		// according to new charge density.
		// mohan add 2009-01-23
		//en.calculate_harris(2);

		if(GlobalV::MY_STOGROUP==0)
		{
			Symmetry_rho srho;
			for(int is=0; is < GlobalV::NSPIN; is++)
			{
				srho.begin(is, GlobalC::CHR,GlobalC::rhopw, GlobalC::Pgrid, GlobalC::symm);
			}
		}
		else
		{
#ifdef __MPI
			if(ModuleSymmetry::Symmetry::symm_flag)	MPI_Barrier(MPI_COMM_WORLD);
#endif
		}

		if(GlobalV::MY_STOGROUP == 0)
		{
        	GlobalC::en.deband = GlobalC::en.delta_e();
		}
        ModuleBase::timer::tick(this->classname, "solve");
        return;
    }
    
}