#ifndef XYZ_IO_H
#define XYZ_IO_H
#include <vector>
#include <string>
#include <utility>

namespace ModuleIO
{
    /**
     * ABACUS (ext-)xyz parser
     * see documents of ASE for more details and example: https://wiki.fysik.dtu.dk/ase/ase/io/formatoptions.html#extxyz
     */
    class XYZParser
    {
        public:
            enum PropertyDatatype {
                Integer,
                Real,
                String,
                Logical
            };
            static void write(const std::string& fn,
                              const int natom,
                              const std::string* elem,
                              const double* pos,
                              const std::string& comment = "",
                              const int nprop = 0,
                              const std::string* propname = nullptr,
                              const int* propdim = nullptr,
                              const PropertyDatatype* propdtype = nullptr,
                              const std::string* prop = nullptr);
            static void write(const std::string& fn,
                              const std::vector<std::string>& elem,
                              const std::vector<double>& pos,
                              const std::string& comment = "",
                              const std::vector<std::string>& propname = {}, 
                              const std::vector<int>& propdim = {}, 
                              const std::vector<PropertyDatatype>& propdtype = {},
                              const std::vector<std::vector<std::string>>& prop = {});
            static void read(const std::string& fn,
                             int& natom,
                             std::string*& elem,
                             double*& pos,
                             std::string& comment,
                             int& nprop,
                             std::string*& propname,
                             int*& propdim,
                             PropertyDatatype*& propdtype,
                             std::string*& prop);
            static void read(const std::string& fn,
                             std::vector<std::string>& elem,
                             std::vector<double>& pos,
                             std::string& comment,
                             std::vector<std::string>& propname,
                             std::vector<int>& propdim,
                             std::vector<PropertyDatatype>& propdtype,
                             std::vector<std::vector<std::string>>& prop);

    }; // class XYZParser
} // namespace ModuleIO
#endif