/**
 * calculate the <phi_i|Lx/Ly/Lz|phi_j> matrix elements
 * 
 * Formulation
 * -----------
 * 
 * Calculate the <phi_i|Lx/Ly/Lz|phi_j> with ladder operator L+ and L-.
 * 
 * The relation between Lx, Ly and L+, L- are:
 * 
 *                                  Lx = (L+ + L-) / 2
 *                                  Ly = (L+ - L-) / 2i
 * 
 * With L+, the spherical harmonic function Ylm (denoted as |l, m> in the following)
 * can be raised:
 * 
 *                          L+|l, m> = sqrt((l-m)(l+m+1))|l, m+1>
 * 
 * Likely, with L-, the spherical harmonic function Ylm can be lowered:
 * 
 *                          L-|l, m> = sqrt((l+m)(l-m+1))|l, m-1>
 * 
 * Therefore the Lx matrix element can be calculated as:
 * 
 *                <l, m|Lx|l, m'> =   sqrt((l-m)(l+m+1)) * delta(m, m'+1) / 2
 *                                  + sqrt((l+m)(l-m+1)) * delta(m, m'-1) / 2
 * 
 * The Ly matrix element can be calculated as:
 * 
 *                <l, m|Ly|l, m'> =   sqrt((l-m)(l+m+1)) * delta(m, m'+1) / 2i
 *                                  - sqrt((l+m)(l-m+1)) * delta(m, m'-1) / 2i
 * 
 * The Lz matrix element can be calculated as:
 * 
 *                          <l, m|Lz|l, m'> = m * delta(m, m')
 * 
 * However, things will change when there are more than one centers.
 * 
 * Technical Details
 * -----------------
 * 
 * 0. The calculation of matrix elements involves the two-center-integral calculation,
 *    this is supported by ABACUS built-in class TwoCenterIntegrator.
 *    see: source/module_basis/module_nao/two_center_integrator.h.
 *  
 * 1. The interface of it is RadialCollection, which is a collection of radial functions. 
 *    see: source/module_basis/module_nao/radial_collection.h
 * 
 * 2. The radial functions are stored in AtomicRadials class,
 *    see: source/module_basis/module_nao/atomic_radials.h
 * 
 * 3. The construction of AtomicRadials involves the filename of orbital, it is stored
 *    in the UnitCell instance
 * 
 * 
 */
#include <cmath>
#include <vector>
#include <map>
#include <tuple>
#include <complex>
#include "module_cell/unitcell.h"
#include "module_basis/module_nao/two_center_integrator.h"
namespace ModuleIO
{

} // namespace ModuleIO