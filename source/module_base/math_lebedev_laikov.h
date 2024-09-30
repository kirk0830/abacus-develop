#ifndef MATH_LEBEDEV_LAIKOV_H
#define MATH_LEBEDEV_LAIKOV_H

#include "module_base/vector3.h"
#include <set>
#include <string>

namespace ModuleBase
{
/**
 * Lebedev-Laikov Spherical Grid Generator
 * 
 * Generate Lebedev-Laikov spherical grid points and weights for numerical integration
 * on the unit sphere. 
 * 
 * 
 */
class Lebedev_laikov_grid
{
public:
    /**
     * @brief Build a new Lebedev-Laikov grid
     * 
     * @param degree the degree of the grid, can only take the following values:
     * 6, 14, 26, 38, 50, 74, 86, 110, 146, 170, 194, 230, 266, 302, 350, 434, 590, 770, 974,
     * 1202, 1454, 1730, 2030, 2354, 2702, 3074, 3470, 3890, 4334, 4802, 5294, 5810
     */
    Lebedev_laikov_grid(const int degree);
    ~Lebedev_laikov_grid();

    /**
     * @brief Generate the grid points
     * 
     */
    void generate_grid_points();

    /**
     * @brief Get all the grid coordinates in ModuleBase::Vector3<double>
     * 
     * @return const ModuleBase::Vector3<double>* 
     */
    const ModuleBase::Vector3<double>* get_grid_coor() const
    {
        return grid_coor;
    };

    /**
     * @brief Get the weight object, for numerical integration
     * 
     * @return const double* weights of the grid points
     */
    const double* get_weight() const
    {
        return weight;
    };

    void print_grid_and_weight(std::string filename);

    // degree: can only take the following values
    // degree = { 6, 14, 26, 38, 50, 74, 86, 110, 146, 170, 194, 230, 266, 302, 350, 434, 590, 770, 974, 
    // 1202, 1454, 1730, 2030, 2354, 2702, 3074, 3470, 3890, 4334, 4802, 5294, 5810};
    int degree = 6;

private:
    int getLebedevReccurencePoints(int type, int start, double a, double b, double v);

    std::set<int> allowed_degree = { 
        6, 14, 26, 38, 50, 74, 86, 110, 146, 170, 194, 230, 266, 
        302, 350, 434, 590, 770, 974, 1202, 1454, 1730, 2030, 
        2354, 2702, 3074, 3470, 3890, 4334, 4802, 5294, 5810
    };

    ModuleBase::Vector3<double> *grid_coor = nullptr;
    double* weight = nullptr;
};

}

#endif