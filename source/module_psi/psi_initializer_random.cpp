#include "psi_initializer_random.h"
#include <mpi.h>
#include "module_base/parallel_global.h"
#include "module_cell/parallel_kpoints.h"

psi_initializer_random::psi_initializer_random(Structure_Factor* sf_in, ModulePW::PW_Basis_K* pw_wfc_in) : psi_initializer(sf_in, pw_wfc_in)
{
    this->ixy2is = new int[this->pw_wfc->fftnxy];
    this->pw_wfc->getfftixy2is(this->ixy2is);
}

psi_initializer_random::~psi_initializer_random()
{
    delete[] this->ixy2is;
}

#ifdef __MPI
void psi_initializer_random::stick_to_pool(float* stick, const int& ir, float* out, const ModulePW::PW_Basis_K* wfc_basis) const
{	
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
		MPI_Recv(stick, nz, MPI_FLOAT, 0, ir, POOL_WORLD, &ierror);
		for(int iz=0; iz<nz; iz++)
		{
			out[is*nz+iz] = stick[iz];
		}
	}
	else if(GlobalV::RANK_IN_POOL==0)
	{
		MPI_Send(stick, nz, MPI_FLOAT, ip, ir, POOL_WORLD);
	}

	return;	
}
void psi_initializer_random::stick_to_pool(double* stick, const int& ir, double* out, const ModulePW::PW_Basis_K* wfc_basis) const
{	
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
}
#endif

void psi_initializer_random::random(std::complex<double>* psi,
                       const int iw_start,
                       const int iw_end,
                       const int ik,
                       const ModulePW::PW_Basis_K* wfc_basis)
{
    this->random_t(psi, iw_start, iw_end, ik, wfc_basis);
}

void psi_initializer_random::random(std::complex<float>* psi,
                       const int iw_start,
                       const int iw_end,
                       const int ik,
                       const ModulePW::PW_Basis_K* wfc_basis)
{
    this->random_t(psi, iw_start, iw_end, ik, wfc_basis);
}

template <typename FPTYPE>
void psi_initializer_random::random_t(std::complex<FPTYPE>* psi,
                         const int iw_start,
                         const int iw_end,
                         const int ik,
                         const ModulePW::PW_Basis_K* wfc_basis)
{
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
}

void psi_initializer_random::initialize(psi::Psi<std::complex<double>>& psi, int ik)
{
    this->random(psi.get_pointer(), 0, psi.get_nbands(), ik, this->pw_wfc);
    // we still need to diagonalize the obtained psi from hsolver::DiagoIterAssist::diagH_subspace
    // will do it in HSolver function...
}