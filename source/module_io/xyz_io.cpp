#include "module_io/xyz_io.h"
#include "module_base/formatter.h"

// write string Lattice="5.44 0.0 0.0 0.0 5.44 0.0 0.0 0.0 5.44"
std::string _write_lattice(const double* eij,
                           const int ndigits = 2)
{
    std::string fmtstr = "%." + std::to_string(ndigits) + "f";
    std::string out = "Lattice=";
    for (int i = 0; i < 9; i++)
    {
        out += FmtCore::format(fmtstr.c_str(), eij[i]);
        if (i < 8)
        {
            out += " ";
        }
    }
    return out;
}

// write string Properties="species:S:1:pos:R:3:vel:R:3:select:I:1"
std::string _write_item_header(const int n,
                               const std::string* name,
                               const ModuleIO::XYZParser::PropertyDatatype* dtype,
                               const int* ncol)
{
    // predefined datatype id
    const std::vector<std::string> dtype_id = {"I", "R", "S", "L"};
    std::vector<std::string> frags(3*n);
    for (int i = 0; i < n; i++)
    {
        frags[3*i] = name[i];
        frags[3*i+1] = dtype_id[static_cast<int>(dtype[i])];
        frags[3*i+2] = std::to_string(ncol[i]);
    }
    const std::string out = "Properties=\"" + FmtCore::join(":", frags) + "\"";
    return out;
}

std::string _wash_table(const std::string table)
{
    std::vector<std::string> rows = FmtCore::split("\n", table);
    // what should be removed is the title, and three lines of dashes, which are the first line,
    // the line below title and the last line
    for (int i = 0; i < 4; i++)
    {
        rows.erase(rows.begin());
    }
    // remove the last line
    rows.pop_back();
    return FmtCore::join("\n", rows);
}

