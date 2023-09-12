#ifndef WAVEFUNC_IN_PW_H
#define WAVEFUNC_IN_PW_H

#include "module_base/complexmatrix.h"
#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_base/realarray.h"
#include "module_base/vector3.h"
#include "module_basis/module_pw/pw_basis_k.h"
#include "module_hamilt_pw/hamilt_pwdft/structure_factor.h"
#include "module_base/spherical_bessel_transformer.h"
#include "module_psi/psi.h"
//---------------------------------------------------
// FUNCTION: expand the local basis sets into plane
// wave basis sets
//---------------------------------------------------
namespace Wavefunc_in_pw
{

	void cal_ovlp_flzjlq(
		std::vector<std::string> &orbital_files, 
		ModuleBase::realArray &ovlp_flzjlq);

	void integral(
		const int meshr, // number of mesh points 
		const double *psir,
		const double *r,
		const double *rab, 
		const int &l, 
		double* table);
	
	//mohan add 2010-04-20
	double cosine_interpolation(
		const double &energy_x,
		const double &ecut,
		const double &beta);

    void init_psi_lcao(const int& ik,
                     const ModulePW::PW_Basis_K* wfc_basis,
                     const Structure_Factor& sf,
                     psi::Psi<std::complex<double>>& psi,
                     const ModuleBase::realArray& ovlp_flzjlq);
	ModuleBase::SphericalBesselTransformer sbt;
}
#endif
