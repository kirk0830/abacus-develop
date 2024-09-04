#ifndef XYZ_IO_H
#define XYZ_IO_H
#include <vector>
#include <string>
#include <utility>
#include <memory>
#include "module_base/formatter.h"

namespace ModuleIO
{
    /**
     * ABACUS (ext-)xyz parser for reading and writing (extended) xyz files.
     * A typical extxyz file looks like:
     * 8
     * Lattice="5.44 0.0 0.0 0.0 5.44 0.0 0.0 0.0 5.44" Properties=species:S:1:pos:R:3 Time=0.0
     * Si        0.00000000      0.00000000      0.00000000
     * Si        1.36000000      1.36000000      1.36000000
     * ...
     * With more properties, the Properties=... part will be longer, for example:
     * Properties="species:S:1:pos:R:3:vel:R:3:select:I:1"
     * then the rest of extxyz file can be:
     * Si        4.08000000      4.08000000      1.36000000   0.00000000      0.00000000      0.00000000       1
     * 
     * So generally the extxyz is actually a table, in which the source data may have different
     * datatypes, such as int, double, std::string and bool.
     * 
     * A online document for extxyz format can be found at:
     * https://wiki.fysik.dtu.dk/ase/ase/io/formatoptions.html#extxyz
     */
    class XYZParser
    {
        public:
            // there are only limited types of properties in xyz file
            enum PropertyDatatype {
                Integer,
                Real,
                String,
                Logical
            };

            XYZParser() = delete;
            XYZParser(const std::string& fn); // read
            XYZParser(const std::string& fn,  // write 
                      const int natom,
                      const std::vector<std::string>& items,
                      const std::vector<std::string>& fmts,
                      const std::vector<int> ncols);

        private:
            std::unique_ptr<FmtTable> table_;

    }; // class XYZParser
} // namespace ModuleIO
#endif