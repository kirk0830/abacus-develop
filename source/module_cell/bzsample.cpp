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
class BrillouinZoneSampler
{
    explicit BrillouinZoneSampler(const std::string& fkpt);
    BrillouinZoneSampler(const ModuleBase::Vector3<double>& kspacing,
                         const ModuleBase::Vector3<double>& a,
                         const ModuleBase::Vector3<double>& b,
                         const ModuleBase::Vector3<double>& c);
    BrillouinZoneSampler();
    void build(std::vector<ModuleBase::Vector3<double>>& vecs, 
               std::vector<double>& weights, 
               std::vector<int>& segids);
    
};
#endif