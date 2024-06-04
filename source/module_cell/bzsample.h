/**
 * @file bzsample.cpp
 * @author kirk0830
 * @brief for Brillouin zone sampling
 * @version 0.1
 * @date 2024-06-04
 * 
 * @copyright Copyright (c) 2024
 * 
 */

// in: file with param: kspacing, gamma_only
// out: vecs (only provide the one in direct coords), weights, segids
#ifndef BZSAMPLE_H
#define BZSAMPLE_H

#include <vector>
#include <string>
#include <fstream>
#include "module_base/vector3.h"
#include "mpi.h"
struct BillouinZonePointSet
{
    // direct coordinates of kpoints, of dim nk*3
    std::vector<ModuleBase::Vector3<double>> vecs;
    // weights of kpoints, of dim nk
    std::vector<double> weights;
    // segment ids of kpoints, of dim nk
    std::vector<int> segids;
    // index of the pool of which the kpoint belongs to
    int ipool;
};

class BrillouinZoneSampler
{
public:
    // from KPT file
    explicit BrillouinZoneSampler(const std::string& fkpt);
    // from INPUT the kspacing
    BrillouinZoneSampler(const ModuleBase::Vector3<double>& kspacing,
                         const ModuleBase::Vector3<double>& a,
                         const ModuleBase::Vector3<double>& b,
                         const ModuleBase::Vector3<double>& c);
    // gamma point only or default
    BrillouinZoneSampler();

    // split kpoints into several groups/pools according to the number of pools
    int cal_pool_color(const int& npools, const int& iproc, const int& nproc)
    {
        const int nproc_per_pool = nproc / npools;
        int color = -1;
        if (iproc < nproc_per_pool * npools)
        {
            color = iproc / nproc_per_pool;
        }
        else
        {
            color = nproc_per_pool + (iproc - nproc_per_pool * npools) / (nproc_per_pool + 1);
        }
        return color;
    }
    // build function, returns set of kpoints (packed-up in BillouinZonePointSet struct) according to the input
    // considerations are mainly on kpoint pools, or say the kpar
    BillouinZonePointSet build(const int& iproc);

};
#endif