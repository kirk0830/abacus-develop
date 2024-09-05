#include "snap_psibeta_half_tddft.h"

#include "module_base/constants.h"
#include "module_base/math_integral.h"
#include "module_base/timer.h"
#include "module_base/ylm.h"

namespace module_tddft
{

// nlm[0] : <phi|exp^{-iAr}|beta>
// nlm[1, 2, 3,] : <phi|r_a * exp^{-iAr}|beta>, which a = x, y, z.
void snap_psibeta_half_tddft(const LCAO_Orbitals& orb,
                             const InfoNonlocal& infoNL_,
                             std::vector<std::vector<std::complex<double>>>& nlm,
                             const ModuleBase::Vector3<double>& R1,
                             const int& T1,
                             const int& L1,
                             const int& m1,
                             const int& N1,
                             const ModuleBase::Vector3<double>& R0, // The projector.
                             const int& T0,
                             const ModuleBase::Vector3<double>& A,
                             const bool& calc_r)
{
    ModuleBase::timer::tick("module_tddft", "snap_psibeta_half_tddft");

    // find number of projectors on atom R0
    const int nproj = infoNL_.nproj[T0];
    if (nproj == 0)
    {
        if (calc_r)
        {
            nlm.resize(4);
        }
        else
        {
            nlm.resize(1);
        }
        return;
    }

    std::vector<bool> calproj;
    calproj.resize(nproj);
    std::vector<int> rmesh1;
    rmesh1.resize(nproj);

    if (calc_r)
    {
        nlm.resize(4);
    }
    else
    {
        nlm.resize(1);
    }

    // Count number of projectors (l,m)
    int natomwfc = 0;
    for (int ip = 0; ip < nproj; ip++)
    {
        //============================
        // Use pseudo-atomic orbitals
        //============================

        const int L0 = infoNL_.Beta[T0].Proj[ip].getL(); // mohan add 2021-05-07
        natomwfc += 2 * L0 + 1;
    }

    for (int dim = 0; dim < nlm.size(); dim++)
    {
        nlm[dim].resize(natomwfc);
        for (auto& x: nlm[dim])
        {
            x = 0.0;
        }
    }

    // rcut of orbtials and projectors
    // in our calculation, we always put orbital phi at the left side of <phi|beta>
    // because <phi|beta> = <beta|phi>
    const double Rcut1 = orb.Phi[T1].getRcut();
    const ModuleBase::Vector3<double> dRa = R0 - R1;

    double distance10 = dRa.norm();

    bool all_out = true;
    for (int ip = 0; ip < nproj; ip++)
    {
        const double Rcut0 = infoNL_.Beta[T0].Proj[ip].getRcut();
        if (distance10 > (Rcut1 + Rcut0))
        {
            calproj[ip] = false;
        }
        else
        {
            all_out = false;
            calproj[ip] = true;
        }
    }

    if (all_out)
    {
        ModuleBase::timer::tick("module_tddft", "snap_psibeta_half_tddft");
        return;
    }

    auto Polynomial_Interpolation
        = [](const int& mesh_r, const double* psi_r, const double* r_radial, const double& x) -> double {
        int left = 0;
        int right = mesh_r - 1;
        while (left <= right)
        {
            int mid = left + (right - left) / 2;
            if (r_radial[mid] == x)
            {
                left = mid;
                right = mid;
                return psi_r[mid];
            }
            else if (r_radial[mid] < x)
            {
                left = mid + 1;
            }
            else
            {
                right = mid - 1;
            }
        }

        double y = 0.0;
        if (right > mesh_r - 4) {
            return y;
        }

        double x0 = r_radial[right];
        double x1 = r_radial[right + 1];
        double x2 = r_radial[right + 2];
        double x3 = r_radial[right + 3];

        double y0 = psi_r[right];
        double y1 = psi_r[right + 1];
        double y2 = psi_r[right + 2];
        double y3 = psi_r[right + 3];

        y = (x - x1) * (x - x2) * (x - x3) / (x0 - x1) / (x0 - x2) / (x0 - x3) * y0
            + (x - x0) * (x - x2) * (x - x3) / (x1 - x0) / (x1 - x2) / (x1 - x3) * y1
            + (x - x0) * (x - x1) * (x - x3) / (x2 - x0) / (x2 - x1) / (x2 - x3) * y2
            + (x - x0) * (x - x1) * (x - x2) / (x3 - x0) / (x3 - x1) / (x3 - x2) * y3;

        return y;
    };

    const int mesh_r1 = orb.Phi[T1].PhiLN(L1, N1).getNr();
    const double* psi_1 = orb.Phi[T1].PhiLN(L1, N1).getPsi();
    const double* radial1 = orb.Phi[T1].PhiLN(L1, N1).getRadial();

    int ridial_grid_num = 140;
    int angular_grid_num = 110;
    double* r_ridial = new double[ridial_grid_num];
    double* weights_ridial = new double[ridial_grid_num];

    int index = 0;
    for (int nb = 0; nb < nproj; nb++)
    {
        const int L0 = infoNL_.Beta[T0].Proj[nb].getL();
        if (!calproj[nb])
        {
            index += 2 * L0 + 1;
            continue;
        }

        const int mesh_r0 = infoNL_.Beta[T0].Proj[nb].getNr();
        const double* beta_r = infoNL_.Beta[T0].Proj[nb].getBeta_r();
        const double* radial0 = infoNL_.Beta[T0].Proj[nb].getRadial();

        double Rcut0 = infoNL_.Beta[T0].Proj[nb].getRcut();
        ModuleBase::Integral::Gauss_Legendre_grid_and_weight(radial0[0],
                                                             radial0[mesh_r0 - 1],
                                                             ridial_grid_num,
                                                             r_ridial,
                                                             weights_ridial);

        double A_phase = A * R0;
        std::complex<double> exp_iAR0 = std::exp(ModuleBase::IMAG_UNIT * A_phase);

        for (int ir = 0; ir < ridial_grid_num; ir++)
        {
            std::vector<std::complex<double>> result_angular(2 * L0 + 1, 0.0);
            std::vector<std::complex<double>> result_angular_r_commu_x;
            std::vector<std::complex<double>> result_angular_r_commu_y;
            std::vector<std::complex<double>> result_angular_r_commu_z;
            if (calc_r)
            {
                result_angular_r_commu_x.resize(2 * L0 + 1, 0.0);
                result_angular_r_commu_y.resize(2 * L0 + 1, 0.0);
                result_angular_r_commu_z.resize(2 * L0 + 1, 0.0);
            }

            for (int ian = 0; ian < angular_grid_num; ian++)
            {
                double x = ModuleBase::Integral::LebedevLaikovGrid110_x[ian];
                double y = ModuleBase::Integral::LebedevLaikovGrid110_y[ian];
                double z = ModuleBase::Integral::LebedevLaikovGrid110_z[ian];
                double weights_angular = ModuleBase::Integral::LebedevLaikovGrid110_w[ian];
                ModuleBase::Vector3<double> r_angular_tmp(x, y, z);

                ModuleBase::Vector3<double> r_coor = r_ridial[ir] * r_angular_tmp;
                ModuleBase::Vector3<double> tmp_r_coor = r_coor + dRa;
                double tmp_r_coor_norm = tmp_r_coor.norm();
                ModuleBase::Vector3<double> tmp_r_unit;
                if (tmp_r_coor_norm > 1e-10)
                {
                    tmp_r_unit = tmp_r_coor / tmp_r_coor_norm;
                }

                if (tmp_r_coor_norm > Rcut1) {
                    continue;
                }

                std::vector<double> rly0;
                ModuleBase::Ylm::rl_sph_harm(L0, x, y, z, rly0);

                std::vector<double> rly1;
                ModuleBase::Ylm::rl_sph_harm(L1, tmp_r_unit.x, tmp_r_unit.y, tmp_r_unit.z, rly1);

                double phase = A * r_coor;
                std::complex<double> exp_iAr = std::exp(ModuleBase::IMAG_UNIT * phase);

                ModuleBase::Vector3<double> tmp_r_coor_r_commu = r_coor + R0;

                for (int m0 = 0; m0 < 2 * L0 + 1; m0++)
                {
                    double temp_interpolation_value = Polynomial_Interpolation(mesh_r1, psi_1, radial1, tmp_r_coor_norm);

                    result_angular[m0] += exp_iAr * rly0[L0 * L0 + m0] * rly1[L1 * L1 + m1]
                                          * temp_interpolation_value
                                          * weights_angular;

                    if (calc_r)
                    {
                        result_angular_r_commu_x[m0] += exp_iAr * tmp_r_coor_r_commu.x * rly0[L0*L0+m0] * rly1[L1*L1+m1]
                                                        * temp_interpolation_value
                                                        * weights_angular;
                        
                        result_angular_r_commu_y[m0] += exp_iAr * tmp_r_coor_r_commu.y * rly0[L0*L0+m0] * rly1[L1*L1+m1]
                                                        * temp_interpolation_value
                                                        * weights_angular;
                        
                        result_angular_r_commu_z[m0] += exp_iAr * tmp_r_coor_r_commu.z * rly0[L0*L0+m0] * rly1[L1*L1+m1]
                                                        * temp_interpolation_value
                                                        * weights_angular;
                    }
                }
            }

            int index_tmp = index;
            double temp = Polynomial_Interpolation(mesh_r0, beta_r, radial0, r_ridial[ir]) * r_ridial[ir] * weights_ridial[ir];
            if (!calc_r)
            {
                for (int m0 = 0; m0 < 2 * L0 + 1; m0++)
                {
                    nlm[0][index_tmp] += temp * result_angular[m0] * exp_iAR0;
                    index_tmp++;
                }
            }
            else
            {
                for (int m0 = 0; m0 < 2 * L0 + 1; m0++)
                {
                    nlm[0][index_tmp] += temp * result_angular[m0] * exp_iAR0;
                    nlm[1][index_tmp] += temp * result_angular_r_commu_x[m0] * exp_iAR0;
                    nlm[2][index_tmp] += temp * result_angular_r_commu_y[m0] * exp_iAR0;
                    nlm[3][index_tmp] += temp * result_angular_r_commu_z[m0] * exp_iAR0;
                    index_tmp++;
                }
            }
        }

        index += 2 * L0 + 1;
    }

    for(int dim = 0; dim < nlm.size(); dim++)
    {
        for (auto &x : nlm[dim])
        {
            // nlm[0] is <phi|exp^{-iAr}|beta>
            // nlm[1 or 2 or 3] is <phi|r_a * exp^{-iAr}|beta>, a = x, y, z
            x = std::conj(x); 
        }
    }

    delete[] r_ridial;
    delete[] weights_ridial;
    assert(index == natomwfc);
    ModuleBase::timer::tick("module_tddft", "snap_psibeta_half_tddft");

    return;
}

} // namespace module_tddft
