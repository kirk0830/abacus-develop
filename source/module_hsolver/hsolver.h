#ifndef HSOLVER_H
#define HSOLVER_H

#include "diagh.h"
#include "module_elecstate/elecstate.h"
#include "module_hamilt/hamilt.h"
#include "module_psi/psi.h"
#include "src_pw/sto_wf.h"

#include <complex>

namespace hsolver
{

class HSolver
{
  public:
    HSolver(){};
    virtual ~HSolver(){
        delete pdiagh;
    };
    /*//initialization, used in construct function or restruct a new HSolver
    virtual void init(
        const Basis* pbas //We need Basis class here, use global class for this initialization first
        //const Input &in, //We need new Input class here, use global variable for this initialization first
        //elecstate::ElecState *pes
        )=0;
    //initialization, only be called for change some parameters only
    virtual void update(
        Input &in )=0;*/

    // solve Hamiltonian to electronic density in ElecState
    virtual void solve
    (
        hamilt::Hamilt* phm, 
        psi::Psi<std::complex<double>>& ppsi, 
        elecstate::ElecState* pes, 
        const std::string method, 
        const bool skip_charge=false
    )
    {
        return;
    }
    virtual void solve
    (
        hamilt::Hamilt* phm, 
        psi::Psi<double>& ppsi, 
        elecstate::ElecState* pes, 
        const std::string method, 
        const bool skip_charge=false
    )
    {
        return;
    }

    virtual void solve
    (
        hamilt::Hamilt* phm, 
        psi::Psi<std::complex<double>>& ppsi, 
        elecstate::ElecState* pes, 
        Stochastic_WF& stowf,
        const int iter,
        const std::string method, 
        const bool skip_charge=false
    )
    {
        return;
    }

    std::string classname = "none";
    // choose method of DiagH for solve Hamiltonian matrix
    // cg, dav, elpa, scalapack, hpseps, cusolver
    std::string method = "none";

  protected:
    DiagH* pdiagh = nullptr; // for single Hamiltonian matrix diagonal solver

};

} // namespace hsolver
#endif