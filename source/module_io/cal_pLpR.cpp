#include <cmath>
#include <vector>
#include <map>
#include <tuple>
#include <complex>
#include <fstream>
#include <memory>
#include "module_cell/unitcell.h"
#include "module_base/spherical_bessel_transformer.h"
#include "module_basis/module_nao/two_center_integrator.h"
#include "module_cell/module_neighbor/sltk_grid_driver.h"
#include "module_cell/module_neighbor/sltk_atom_arrange.h"
#include "module_parameter/parameter.h"

void _radial_indexing(const UnitCell& ucell,
                      std::vector<std::tuple<int, int, int, int, int>>& lin2comp,
                      std::map<std::tuple<int, int, int, int, int>, int>& comp2lin)
{
    int i = 0;
    for (int it = 0; it < ucell.ntype; it++)
    {
        const Atom& atom = ucell.atoms[it];
        const int lmax = atom.nwl;
        for (int ia = 0; ia < atom.na; ia++)
        {
            for (int l = 0; l < lmax; l++)
            {
                const int nzeta = atom.l_nchi[l];
                for (int iz = 0; iz < nzeta; iz++)
                {
                    for (int m = -l; m <= l; m++)
                    {
                        lin2comp.push_back(std::make_tuple(it, ia, l, iz, m));
                        comp2lin[std::make_tuple(it, ia, l, iz, m)] = i;
                        i++;
                    }
                }
            }
        }
    }
}

/**
 * 
 * FIXME: the following part will be transfered to TwoCenterIntegrator soon
 * 
 */

// L+|l, m> = sqrt((l-m)(l+m+1))|l, m+1>, return the sqrt((l-m)(l+m+1))
double _lplus_on_ylm(const int l, const int m)
{
    return std::sqrt((l - m) * (l + m + 1));
}

// L-|l, m> = sqrt((l+m)(l-m+1))|l, m-1>, return the sqrt((l+m)(l-m+1))
double _lminus_on_ylm(const int l, const int m)
{
    return std::sqrt((l + m) * (l - m + 1));
}

void _cal_pLzpR(const std::unique_ptr<TwoCenterIntegrator>& calculator,
                const int it, const int ia, const int il, const int iz, const int mi,
                const int jt, const int ja, const int jl, const int jz, const int mj,
                const ModuleBase::Vector3<double>& vR,
                std::complex<double>& val)
{
    double val_ = 0;
    calculator->calculate(it, il, iz, mi, jt, jl, jz, mj, vR, &val_);
    val = std::complex<double>(mi) * val_;
}

void _cal_pLypR(const std::unique_ptr<TwoCenterIntegrator>& calculator,
                const int it, const int ia, const int il, const int iz, const int im,
                const int jt, const int ja, const int jl, const int jz, const int jm,
                const ModuleBase::Vector3<double>& vR,
                std::complex<double>& val)
{
    // Ly = -i/2 * (L+ - L-)
    const double plus_ = _lplus_on_ylm(jl, jm);
    const double minus_ = _lminus_on_ylm(jl, jm);
    double val_plus = 0, val_minus = 0;
    if (plus_ != 0)
    {
        calculator->calculate(it, il, iz, im, jt, jl, jz, jm + 1, vR, &val_plus);
        val_plus *= plus_;
    }
    if (minus_ != 0)
    {
        calculator->calculate(it, il, iz, im, jt, jl, jz, jm - 1, vR, &val_minus);
        val_minus *= minus_;
    }
    val = std::complex<double>(0, -0.5) * (val_plus - val_minus);
}

void _cal_pLxpR(const std::unique_ptr<TwoCenterIntegrator>& calculator,
                const int it, const int ia, const int il, const int iz, const int im,
                const int jt, const int ja, const int jl, const int jz, const int jm,
                const ModuleBase::Vector3<double>& vR,
                std::complex<double>& val)
{   
    // Lx = 1/2 * (L+ + L-)
    const double plus_ = _lplus_on_ylm(jl, jm);
    const double minus_ = _lminus_on_ylm(jl, jm);
    double val_plus = 0, val_minus = 0;
    if (plus_ != 0)
    {
        calculator->calculate(it, il, iz, im, jt, jl, jz, jm + 1, vR, &val_plus);
        val_plus *= plus_;
    }
    if (minus_ != 0)
    {
        calculator->calculate(it, il, iz, im, jt, jl, jz, jm - 1, vR, &val_minus);
        val_minus *= minus_;
    }
    val = std::complex<double>(0.5) * (val_plus + val_minus);
}

void _cal_pLpR(const std::unique_ptr<TwoCenterIntegrator>& calculator,
               const int it, const int ia, const int il, const int iz, const int im,
               const int jt, const int ja, const int jl, const int jz, const int jm,
               const ModuleBase::Vector3<double>& vR,
               std::complex<double>& val,
               const char dir = 'z')
{
    switch (dir)
    {
    case 'z':
        _cal_pLzpR(calculator, it, ia, il, iz, im, jt, ja, jl, jz, jm, vR, val);
        break;
    case 'y':
        _cal_pLypR(calculator, it, ia, il, iz, im, jt, ja, jl, jz, jm, vR, val);
        break;
    case 'x':
        _cal_pLxpR(calculator, it, ia, il, iz, im, jt, ja, jl, jz, jm, vR, val);
        break;
    default:
        break;
    }
}

// calculate the matrix of <phi_i|Lx/Ly/Lz|phi_j>
void calculate(const std::unique_ptr<TwoCenterIntegrator>& calculator,
               const UnitCell& ucell,
               const double rcut_max,
               std::vector<double>& out,
               std::ofstream* ptrlog = nullptr)
{
    // I follow the void ESolver_LJ::runner to find the neighboring atoms
    const int t_atom = PARAM.inp.test_atom_input;
    const int t_deconstructor = PARAM.inp.test_deconstructor;
    const int t_grid = PARAM.inp.test_grid;
    Grid_Driver neighbor_searcher(t_deconstructor, t_grid);
    atom_arrange::search(PARAM.inp.search_pbc, 
                         *(ptrlog), 
                         neighbor_searcher, 
                         ucell, 
                         rcut_max,
                         t_atom);
    
    // calculate the matrix elements
    
}

/**
 * @brief Write a 2D matrix to a plain text file.
 * 
 * @tparam T the type of the matrix elements
 * @param matrix the matrix that will be written
 * @param ncols the number of columns of the matrix
 * @param ofs the output file stream
 * @param oflog the log file stream
 * @param rank the identifier to distinguish different MPI processes
 */
template <typename T>
void _write_plain_matrix(const std::vector<T>& matrix, 
                         const int ncols,
                         std::ofstream& ofs,
                         std::ofstream& oflog,
                         const int precision = 8,
                         const int rank = 0)
{
    if (!ofs.good()) {
        return;
    }

    if (rank == 0) {
        for (int i = 0; i < matrix.size(); i++) {
            ofs << std::setw(precision + 6) << std::right << std::scientific << matrix[i];
            if ((i + 1) % ncols == 0) {
                ofs << std::endl;
            } else {
                ofs << " ";
            }
        }
    }
}