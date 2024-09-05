#include <sstream>
#include <cassert>
#include <cmath>
#include <numeric>
#include <map>
#include "module_hamilt_pw/hamilt_pwdft/radial_proj.h"
#include "module_base/constants.h"
#include "module_base/matrix.h"
#include "module_base/math_ylmreal.h"
#include "module_base/spherical_bessel_transformer.h"
#include "module_basis/module_pw/pw_basis.h"

void _backward_cutoff(const int nr,
                      const double* x_in,
                      const double* y_in,
                      double* x_out,
                      double* y_out,
                      const double cut_thr = 1e-8);

void _backward_cutoff(const std::vector<double>& x_in,
                      const std::vector<double>& y_in,
                      std::vector<double>& x_out,
                      std::vector<double>& y_out,
                      const double cthr = 1e-8);

void mask(const std::vector<double>& r,
          const std::vector<double>& w, // the radial function
          const std::vector<double>& m, // the mask function
          std::vector<double>& wm,
          const double cut_thr = 1e-8);

void mask(const int nr,
          const double* r,
          const double* w,
          const int nrm,
          const double* m,
          double* wm,
          const double cut_thr = 1e-8);

void RadialProjection::RadialProjector::_build_backward_map(const std::vector<std::vector<int>>& it2iproj,
                                                            const std::vector<int>& iproj2l,
                                                            std::vector<int>& irow2it,
                                                            std::vector<int>& irow2iproj,
                                                            std::vector<int>& irow2m)
{
    const int ntype = it2iproj.size(); // the ntype here only count the valid, that is, with the projector.

    int nproj_tot = 0; // count the total number of projectors (atom-position-irrelevant)
    for(int it = 0; it < ntype; it++) // for all types with projector...
    {
        for(auto& iproj: it2iproj[it]) // for each projector, the projector is indexed with global iproj
        {
            const int l = iproj2l[iproj]; // easily get the angular momentum of this projector
            nproj_tot += (2*l + 1); // add 2l+1 to the total number of projectors
        }
    }
    // resize/allcoate the memory for the output
    irow2it.resize(nproj_tot);
    irow2iproj.resize(nproj_tot);
    irow2m.resize(nproj_tot);

    int irow = 0;
    for(int it = 0; it < ntype; it++)
    {
        const int nproj = it2iproj[it].size();
        for(int iproj = 0; iproj < nproj; iproj++)
        {
            const int l = iproj2l[it2iproj[it][iproj]];
            for(int m = -l; m <= l; m++)
            {
                irow2it[irow] = it;
                irow2iproj[irow] = iproj;
                irow2m[irow] = m;
                irow++;
            }
        }
    }
}

void RadialProjection::RadialProjector::_build_forward_map(const std::vector<std::vector<int>>& it2ia,
                                                           const std::vector<std::vector<int>>& it2iproj,
                                                           const std::vector<int>& iproj2l,
                                                           std::map<std::tuple<int, int, int, int>, int>& itiaiprojm2irow)
{
    const int ntype = it2ia.size();
    int irow = 0;
    for(int it = 0; it < ntype; it++)
    {
        const int nproj = it2iproj[it].size();
        for(auto& ia: it2ia[it]) // the index is from UnitCell, so it allows the case ia is not continuous
        {
            for(int iproj = 0; iproj < nproj; iproj++)
            {
                const int l = iproj2l[it2iproj[it][iproj]]; // what is the iproj of the i-th projector of it atomtype?
                for(int m = -l; m <= l; m++)
                {
                    itiaiprojm2irow[std::make_tuple(it, ia, iproj, m)] = irow;
                    irow++;
                }
            }
        }
    }
}

void RadialProjection::RadialProjector::build_sbt_tab(const int nr,
                                                      const double* r,
                                                      const std::vector<double*>& radials,
                                                      const std::vector<int>& l,
                                                      const int nq,
                                                      const double& dq,
                                                      const bool do_mask,
                                                      const double& eta,
                                                      const double& cut_thr)
{
    l_ = l;
    const int nrad = radials.size();
    assert(nrad == l.size());
    std::vector<double> qgrid(nq);
    std::iota(qgrid.begin(), qgrid.end(), 0);
    std::transform(qgrid.begin(), qgrid.end(), qgrid.begin(), [dq](const double& q){return q*dq;});

    if(cubspl_.get()) { cubspl_.reset(); } // release the old one if it is not the first time
    cubspl_ = std::unique_ptr<ModuleBase::CubicSpline>(new ModuleBase::CubicSpline(nq,              // int
                                                                                   qgrid.data()));  // double*
    cubspl_->reserve(nrad);
    ModuleBase::SphericalBesselTransformer sbt_(true); // bool: enable cache

    std::vector<double> _temp(nq);
    // the SphericalBesselTransformer's result is multiplied by one extra factor sqrt(2/pi), should remove it
    // see module_base/spherical_bessel_transformer.h and module_base/spherical_bessel_transformer.cpp:328
    const double pref = std::sqrt(2.0/std::acos(-1.0));

    std::vector<double> mask_(201, 1.0);
    if (do_mask) { maskgen(eta, mask_); }
    for(int i = 0; i < nrad; i++)
    {
        std::vector<double> wm(nr);
        mask(nr, r, radials[i], 201, mask_.data(), wm.data(), cut_thr);
        sbt_.direct(l[i], nr, r, wm.data(), nq, qgrid.data(), _temp.data());
        std::for_each(_temp.begin(), _temp.end(), [pref](double& x){x = x/pref;});
        cubspl_->add(_temp.data());
    }
}

void RadialProjection::RadialProjector::build_sbt_tab(const std::vector<double>& r,
                                                      const std::vector<std::vector<double>>& radials,
                                                      const std::vector<int>& l,
                                                      const int nq,
                                                      const double& dq,
                                                      const bool do_mask,
                                                      const double& eta,
                                                      const double& cut_thr)
{
    const int nr = r.size();
    const int nrad = radials.size();
    for(int i = 0; i < nrad; i++) { assert(radials[i].size() == nr); }
    std::vector<double*> radptrs(radials.size());
    for(int i = 0; i < radials.size(); i++) { radptrs[i] = const_cast<double*>(radials[i].data()); }
    build_sbt_tab(nr, r.data(), radptrs, l, nq, dq, do_mask, eta, cut_thr);
}

void RadialProjection::RadialProjector::sbtft(const std::vector<ModuleBase::Vector3<double>>& qs,
                                              std::vector<std::complex<double>>& out,
                                              const char type,
                                              const double& omega,
                                              const double& tpiba)
{
    // first cache the Ylm values
    const int lmax_ = *std::max_element(l_.begin(), l_.end());
    const int total_lm = std::pow(lmax_+1, 2);
    const int npw = qs.size();
    ModuleBase::matrix ylm_(total_lm, npw);
    ModuleBase::YlmReal::Ylm_Real(total_lm, npw, qs.data(), ylm_);

    const int nrad = l_.size();
    int nchannel = 0;
    for(auto l: l_) { nchannel += 2*l+1; }
    out.resize(nchannel*npw);

    std::vector<double> qnorm(npw);
    std::transform(qs.begin(), qs.end(), qnorm.begin(), [tpiba](const ModuleBase::Vector3<double>& q){return tpiba*q.norm();});
    
    std::vector<double> Jlfq(npw);
    int iproj = 0;
    for(int i = 0; i < nrad; i++)
    {
        const int l = l_[i];
        std::complex<double> pref = (type == 'r')? std::pow(ModuleBase::IMAG_UNIT, l) : std::pow(ModuleBase::NEG_IMAG_UNIT, l);
        pref = pref * ModuleBase::FOUR_PI/std::sqrt(omega);
        cubspl_->eval(npw, qnorm.data(), Jlfq.data(), nullptr, nullptr, i);
        for(int m = -l; m <= l; m++)
        {
            for(int iq = 0; iq < npw; iq++)
            {
                out[iproj*npw+iq] = pref * Jlfq[iq] * ylm_(l*l + l + m, iq);
            }
            iproj++;
        }
    }
    assert(iproj == nchannel); // should write to inflate each radial to 2l+1 channels
}

void RadialProjection::maskgen(const double& eta, 
                               std::vector<double>& mask)
{
    /* mask function is hard coded here, eta = 15 */
    mask.resize(201);
    std::string src = "0.10000000E+01";
    src += " 0.10000000E+01 0.99948662E+00 0.99863154E+00 0.99743557E+00";
    src += " 0.99589985E+00 0.99402586E+00 0.99181538E+00 0.98927052E+00";
    src += " 0.98639370E+00 0.98318766E+00 0.97965544E+00 0.97580040E+00";
    src += " 0.97162618E+00 0.96713671E+00 0.96233623E+00 0.95722924E+00";
    src += " 0.95182053E+00 0.94611516E+00 0.94011842E+00 0.93383589E+00";
    src += " 0.92727338E+00 0.92043693E+00 0.91333282E+00 0.90596753E+00";
    src += " 0.89834777E+00 0.89048044E+00 0.88237263E+00 0.87403161E+00";
    src += " 0.86546483E+00 0.85667987E+00 0.84768450E+00 0.83848659E+00";
    src += " 0.82909416E+00 0.81951535E+00 0.80975838E+00 0.79983160E+00";
    src += " 0.78974340E+00 0.77950227E+00 0.76911677E+00 0.75859548E+00";
    src += " 0.74794703E+00 0.73718009E+00 0.72630334E+00 0.71532544E+00";
    src += " 0.70425508E+00 0.69310092E+00 0.68187158E+00 0.67057566E+00";
    src += " 0.65922170E+00 0.64781819E+00 0.63637355E+00 0.62489612E+00";
    src += " 0.61339415E+00 0.60187581E+00 0.59034914E+00 0.57882208E+00";
    src += " 0.56730245E+00 0.55579794E+00 0.54431609E+00 0.53286431E+00";
    src += " 0.52144984E+00 0.51007978E+00 0.49876105E+00 0.48750040E+00";
    src += " 0.47630440E+00 0.46517945E+00 0.45413176E+00 0.44316732E+00";
    src += " 0.43229196E+00 0.42151128E+00 0.41083069E+00 0.40025539E+00";
    src += " 0.38979038E+00 0.37944042E+00 0.36921008E+00 0.35910371E+00";
    src += " 0.34912542E+00 0.33927912E+00 0.32956851E+00 0.31999705E+00";
    src += " 0.31056799E+00 0.30128436E+00 0.29214897E+00 0.28316441E+00";
    src += " 0.27433307E+00 0.26565709E+00 0.25713844E+00 0.24877886E+00";
    src += " 0.24057988E+00 0.23254283E+00 0.22466884E+00 0.21695884E+00";
    src += " 0.20941357E+00 0.20203357E+00 0.19481920E+00 0.18777065E+00";
    src += " 0.18088790E+00 0.17417080E+00 0.16761900E+00 0.16123200E+00";
    src += " 0.15500913E+00 0.14894959E+00 0.14305240E+00 0.13731647E+00";
    src += " 0.13174055E+00 0.12632327E+00 0.12106315E+00 0.11595855E+00";
    src += " 0.11100775E+00 0.10620891E+00 0.10156010E+00 0.97059268E-01";
    src += " 0.92704295E-01 0.88492966E-01 0.84422989E-01 0.80492001E-01";
    src += " 0.76697569E-01 0.73037197E-01 0.69508335E-01 0.66108380E-01";
    src += " 0.62834685E-01 0.59684561E-01 0.56655284E-01 0.53744102E-01";
    src += " 0.50948236E-01 0.48264886E-01 0.45691239E-01 0.43224469E-01";
    src += " 0.40861744E-01 0.38600231E-01 0.36437098E-01 0.34369520E-01";
    src += " 0.32394681E-01 0.30509780E-01 0.28712032E-01 0.26998673E-01";
    src += " 0.25366964E-01 0.23814193E-01 0.22337676E-01 0.20934765E-01";
    src += " 0.19602844E-01 0.18339338E-01 0.17141711E-01 0.16007467E-01";
    src += " 0.14934157E-01 0.13919377E-01 0.12960772E-01 0.12056034E-01";
    src += " 0.11202905E-01 0.10399183E-01 0.96427132E-02 0.89313983E-02";
    src += " 0.82631938E-02 0.76361106E-02 0.70482151E-02 0.64976294E-02";
    src += " 0.59825322E-02 0.55011581E-02 0.50517982E-02 0.46327998E-02";
    src += " 0.42425662E-02 0.38795566E-02 0.35422853E-02 0.32293218E-02";
    src += " 0.29392897E-02 0.26708663E-02 0.24227820E-02 0.21938194E-02";
    src += " 0.19828122E-02 0.17886449E-02 0.16102512E-02 0.14466132E-02";
    src += " 0.12967606E-02 0.11597692E-02 0.10347601E-02 0.92089812E-03";
    src += " 0.81739110E-03 0.72348823E-03 0.63847906E-03 0.56169212E-03";
    src += " 0.49249371E-03 0.43028657E-03 0.37450862E-03 0.32463165E-03";
    src += " 0.28016004E-03 0.24062948E-03 0.20560566E-03 0.17468305E-03";
    src += " 0.14748362E-03 0.12365560E-03 0.10287226E-03 0.84830727E-04";
    src += " 0.69250769E-04 0.55873673E-04 0.44461100E-04 0.34793983E-04";
    src += " 0.26671449E-04 0.19909778E-04 0.14341381E-04 0.98138215E-05";
    std::stringstream ss(src);
    for(int i = 0; i < mask.size(); i++)
    {
        ss >> mask[i];
    }
}

// do cutoff for the radial functions, discard all the values below the cut_thr
void _backward_cutoff(const int nr,
                      const double* x_in,
                      const double* y_in,
                      double* x_out,
                      double* y_out,
                      const double cut_thr)
{
    // this function works in a delete-new manner, which is susceptible to memory leak
    // do not use this function if possible
    delete[] x_out;
    delete[] y_out;
    assert(false); // not implemented
}
// a std::vector overload
void _backward_cutoff(const std::vector<double>& x_in,
                      const std::vector<double>& y_in,
                      std::vector<double>& x_out,
                      std::vector<double>& y_out,
                      const double cthr)
{
    assert(x_in.size() == y_in.size()); // assert the x and y to be homogeneous
    const int n = x_in.size();
    int n_out = n;
    while (n_out > 0 && y_in[n_out-1] <= cthr) { n_out--; }
    x_out.resize(n_out);
    y_out.resize(n_out);
    std::copy(x_in.begin(), x_in.begin()+n_out, x_out.begin());
    std::copy(y_in.begin(), y_in.begin()+n_out, y_out.begin());
}

// wm(r) = w(r)/m(r), in which w(r) is the function to "mask" and m(r) is the so-called mask function
// the mask function is rescaled so that r ranges from 0 to 1
void mask(const int nr,
          const double* r,
          const double* w,
          const int nrm,
          const double* m,
          double* wm,
          const double cut_thr)
{
    std::vector<double> r_(nr), w_(nr);
    std::copy(r, r+nr, r_.begin());
    std::copy(w, w+nr, w_.begin());
    std::vector<double> m_(nrm);
    std::copy(m, m+nrm, m_.begin());
    std::vector<double> wm_(nr);
    mask(r_, w_, m_, wm_, cut_thr);
    std::copy(wm_.begin(), wm_.end(), wm);
}

void mask(const std::vector<double>& r,
          const std::vector<double>& w, // the radial function
          const std::vector<double>& m, // the mask function
          std::vector<double>& wm,
          const double cut_thr)
{
    assert(r.size() == w.size()); // x and y should have the same size
    // for mask function, only y is needed to provided because x is default to be in range [0, 1].
    // then the cutoff of mask function is always the 1.5*rmax, in which rmax is the maximum of r of
    // function to mask.

    // first do the backward cutoff for the radial function to reduce unnecessary computation
    std::vector<double> r_, w_;
    _backward_cutoff(r, w, r_, w_, cut_thr);

    // then build the cubic spline interpolation table with the newly generated r_
    const double rmax = *std::max_element(r_.begin(), r_.end());
    std::vector<double> rmask(m.size(), 0.0);
    std::iota(rmask.begin(), rmask.end(), 0.0);
    // remember the realspace cutoff for mask function is about 1.5*rmax
    std::transform(rmask.begin(), rmask.end(), rmask.begin(), [rmax](double x){return 1.5*x/rmax;});
    // build the cubic spline for the mask function
    ModuleBase::CubicSpline maskspl_(rmask.size(), rmask.data(), m.data());

    // do mask, first evaluate the mask function on the new r_
    std::vector<double> m_(r_.size());
    maskspl_.eval(r_.size(), r_.data(), m_.data());
    wm.resize(r_.size());
    // then do the mask
    std::transform(w_.begin(), w_.end(), m_.begin(), wm.begin(), [](double x, double y){return x/y;});
    // then zero-padding the wm to the original size
    wm.resize(r.size(), 0.0);
}

// Small box FFT then requires another spherical FFT grid that the |q|<2qmax - qc, in which the qmax is
// the norm of wavevector of the maximal kinetic energy cutoff, usually corresponds to that of ecutrho. 
// qc always corresponds to ecutwfc.



// the qmax and qc may change if the volumn of cell changes, otherwise all the wm(q) can be cached 
// somewhere

// transform wm(r) to wm(q). Then do sample and integration on real space grid. The sample of realspace
// grid points can be done by either primitive sampling or FFT by PW_Basis::recip2real.




















/**
 * Additional-bidirectional mapping for the projector. 
 * 
 * These two methods are commented out because of minimal-implementation consideration.
 */

// void build_itiprojm_map(const std::vector<std::vector<int>>& it2iproj,
//                         const std::vector<int>& iproj2l,
//                         std::vector<int>& irow2it,
//                         std::vector<int>& irow2iproj,
//                         std::vector<int>& irow2m,
//                         std::map<std::tuple<int, int, int>, int>& itiprojm2irow)
// {
//     const int ntype = it2iproj.size();

//     int nproj_tot = 0;
//     for(int it = 0; it < ntype; it++)
//     {
//         for(auto& iproj: it2iproj[it])
//         {
//             const int l = iproj2l[iproj];
//             nproj_tot += (2*l + 1);
//         }
//     }
//     irow2it.resize(nproj_tot);
//     irow2iproj.resize(nproj_tot);
//     irow2m.resize(nproj_tot);

//     int irow = 0;
//     for(int it = 0; it < ntype; it++)
//     {
//         const int nproj = it2iproj[it].size();
//         for(int iproj = 0; iproj < nproj; iproj++)
//         {
//             const int l = iproj2l[it2iproj[it][iproj]];
//             for(int m = -l; m <= l; m++)
//             {
//                 irow2it[irow] = it;
//                 irow2iproj[irow] = iproj;
//                 irow2m[irow] = m;
//                 itiprojm2irow[std::make_tuple(it, iproj, m)] = irow;
//                 irow++;
//             }
//         }
//     }
// }

// void build_itiaiprojm_map(const std::vector<std::vector<int>>& it2ia,
//                           const std::vector<std::vector<int>>& it2iproj,
//                           const std::vector<int>& iproj2l,
//                           std::vector<int>& irow2it,
//                           std::vector<int>& irow2ia,
//                           std::vector<int>& irow2iproj,
//                           std::vector<int>& irow2m,
//                           std::map<std::tuple<int, int, int, int>, int>& itiaiprojm2irow)
// {
//     const int ntype = it2ia.size();
//     int nproj_tot = 0;
//     for(int it = 0; it < ntype; it++)
//     {
//         for(auto& ia: it2ia[it])
//         {
//             for(auto& iproj: it2iproj[it])
//             {
//                 const int l = iproj2l[iproj];
//                 nproj_tot += (2*l + 1);
//             }
//         }
//     }
//     irow2it.resize(nproj_tot);
//     irow2ia.resize(nproj_tot);
//     irow2iproj.resize(nproj_tot);
//     irow2m.resize(nproj_tot);

//     int irow = 0;
//     for(int it = 0; it < ntype; it++)
//     {
//         const int nproj = it2iproj[it].size();
//         for(auto& ia: it2ia[it])
//         {
//             for(int iproj = 0; iproj < nproj; iproj++)
//             {
//                 const int l = iproj2l[it2iproj[it][iproj]];
//                 for(int m = -l; m <= l; m++)
//                 {
//                     irow2it[irow] = it;
//                     irow2ia[irow] = ia;
//                     irow2iproj[irow] = iproj;
//                     irow2m[irow] = m;
//                     itiaiprojm2irow[std::make_tuple(it, ia, iproj, m)] = irow;
//                     irow++;
//                 }
//             }
//         }
//     }
// }